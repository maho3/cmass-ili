#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <climits>

#include <vector>
#include <functional>
#include <map>
#include <utility>

#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
using namespace pybind11::literals;
namespace pyb = pybind11;

#include "cuboid.h"

namespace cmangle {
    extern "C" {
        #include "mangle.h"
    }
}

#include "healpix_base.h"

#ifndef M_PIl
#define M_PIl 3.14159265358979323846264338327950288419716939937510L
#endif

typedef uint64_t galid_t;  // Matt: I had to increase this from 32->64
typedef uint64_t hloidx_t;

typedef std::function<pyb::tuple(int, pyb::array_t<hloidx_t>, pyb::array_t<double>)> HOD_fct_t;

static constexpr auto pyb_arr_style = pyb::array::c_style | pyb::array::forcecast;

namespace Geometry
{
    int remaps[][9] =
                      { // 1.4142 1.0000 0.7071
                        { 1, 1, 0,
                          0, 0, 1,
                          1, 0, 0, },
                        // 1.7321 0.8165 0.7071
                        { 1, 1, 1,
                          1, 0, 0,
                          0, 1, 0, },
                        // 1.0000 1.0000 1.0000
                        // trivial case, only for debugging
                        { 1, 0, 0,
                          0, 1, 0,
                          0, 0, 1 },
                      };

    // get from the quadrant ra=[-90,90], dec=[0,90] to the NGC footprint
    // we only need a rotation around the y-axis I believe
    const double alpha = 97.0 * M_PI / 180.0; // rotation around y-axis
    const double beta = 6.0; // rotation around z-axis, in degrees

    // in units of L1, L2, L3
    const double origin[] = { 0.5, -0.058, 0.0 };
}

namespace Numbers
{
    static const double lightspeed = 299792.458; // km/s    

    // both figures from Chang
    // we are dealing with pretty small angle differences so better do things in long double
    static const long double angscale = 0.01722L * M_PIl / 180.0L; // in rad
    static const double collrate = 0.6;

    // before initial downsampling, we increase the fiber collision rate by this
    // to make sure that after fiber collisions we still have enough galaxies in each redshift bin
    static const double fibcoll_rate_correction = 0.05;

    // how many interpolation stencils we use to get from chi to z
    static const int N_interp = 1024;

    // gsl integration workspace size for chi(z) integration
    static const size_t ws_size = 8192;

    // maximum integration error in chi(z), Mpc/h
    static const double chi_epsabs = 1e-3;
}

struct Mask
{
    static const int Nveto = 6;
    // we order them by size so the cheap ones go first
    // it is important that the mask comes first, because the veto masks are inverted
    const char fnames[Nveto+1][128] = {
          "mask_DR12v5_CMASS_North.ply",
          "bright_object_mask_rykoff_pix.ply", 
          "centerpost_mask_dr12.ply", 
          "collision_priority_mask_dr12.ply",
          "badfield_mask_postprocess_pixs8.ply", 
          "allsky_bright_star_mask_pix.ply",
          "badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply",
    };

    std::vector<cmangle::MangleMask *> masks;

    Mask (const char *boss_dir, bool veto=true)
    {
        for (int ii=0; ii<1+veto*Nveto; ++ii)
            masks.push_back(cmangle::mangle_new());

        // pymangle uses the great convention that status=0 means failure
        int status = 1;
        
        #pragma omp parallel for reduction(*:status)
        for (int ii=0; ii<1+veto*Nveto; ++ii)
        {
            int status_;
            char fname[512];
            std::snprintf(fname, 512, "%s/%s", boss_dir, fnames[ii]);
            status_ = cmangle::mangle_read(masks[ii], fname);
            status *= status_;
            if (!status_) continue;
            status_ = cmangle::set_pixel_map(masks[ii]);
            status *= status_;
        }

        if (!status) throw std::runtime_error("Mask loading failed");
    }
    
    ~Mask ()
    {
        for (auto m: masks) cmangle::mangle_free(m);
    }

    bool masked (const cmangle::Point &pt, int &status) const
    {
        status = 1;
        int64_t poly_id; long double weight;
        for (size_t ii=0; ii<masks.size(); ++ii)
        {
            status *= cmangle::mangle_polyid_and_weight_pix(masks[ii], &pt, &poly_id, &weight);
            if ( (weight==0.0L) == (ii==0) )
                return true;
        }
        return false;
    }
};

struct Lightcone
{
    const char *boss_dir;
    const Mask &mask;
    const double BoxSize, Omega_m, zmin, zmax;
    const int remap_case;
    const bool do_downsample, verbose;
    const unsigned augment;
    const unsigned long seed;
    const std::vector<double> snap_times;
    size_t Nsnaps;
    std::vector<double> snap_redshifts, snap_chis, redshift_bounds, chi_bounds;
    gsl_spline *z_chi_interp; std::vector<gsl_interp_accel *> z_chi_interp_acc;
    const cuboid::Cuboid C; const double Li[3]; // the sidelengths in decreasing order
    gsl_histogram *boss_z_hist;
    std::vector<int> snap_indices_done;

    // halo positions after remapping and redshift range cut, as well as IDs
    // the CHI contains the velocity extrapolation correction
    std::vector<double> tmp_hlo_CHI, tmp_hlo_V, tmp_hlo_Z;
    std::vector<hloidx_t> tmp_hlo_ID;

    // galaxy positions and IDs
    std::vector<double> RA, DEC, Z;
    std::vector<galid_t> GALID;

    // the callback into python
    // const HOD_fct_t hod_fct;
    HOD_fct_t hod_fct;

    Lightcone (const Mask &mask_,
               // HOD_fct_t hod_fct_,
               double Omega_m_, double zmin_, double zmax_,
               const std::vector<double> &snap_times_,
               const char *boss_dir_=nullptr, 
               double BoxSize_=3e3, int remap_case_=0,
               bool verbose_=false, unsigned augment_=0,
               unsigned long seed_=137UL) :
        mask{mask_},
        // hod_fct{hod_fct_},
        Omega_m{Omega_m_}, zmin{zmin_}, zmax{zmax_},
        snap_times{snap_times_}, Nsnaps{snap_times_.size()},
        boss_dir{boss_dir_}, do_downsample(boss_dir_),
        BoxSize{BoxSize_}, remap_case{remap_case_},
        verbose{verbose_}, augment{augment_},
        seed{seed_},
        C{cuboid::Cuboid(Geometry::remaps[remap_case_])}, Li{ C.L1, C.L2, C.L3}
    {
        if (Nsnaps>=256) // can't fit in our data type
            throw std::runtime_error("too many snapshots");

        if (verbose) std::printf("process_times\n");
        process_times();

        z_chi_interp = gsl_spline_alloc(gsl_interp_cspline, Numbers::N_interp);
        for (int ii=0; ii<omp_get_max_threads(); ++ii)
            z_chi_interp_acc.push_back(gsl_interp_accel_alloc());
        if (verbose) std::printf("interpolate_chi_z\n");
        interpolate_chi_z();

        // the target redshift distribution
        if (do_downsample)
        {
            if (verbose) std::printf("read_boss_nz\n");
            read_boss_nz();
        }
    }

    void set_hod (HOD_fct_t hod_fct_) { hod_fct = hod_fct_; }

    ~Lightcone ()
    {
        gsl_spline_free(z_chi_interp);
        for (auto x: z_chi_interp_acc) gsl_interp_accel_free(x);
        if (boss_z_hist) gsl_histogram_free(boss_z_hist);
    }

    void _add_snap (int snap_idx, size_t Nhlo,
                    std::vector<double> &xhlo, std::vector<double> &vhlo)
    {
        std::vector<double> delta_xgal, delta_vgal; std::vector<hloidx_t> hlo_idx;
        size_t Ngal;

        if (std::find(snap_indices_done.begin(), snap_indices_done.end(), snap_idx)
            != snap_indices_done.end())
            throw std::runtime_error("adding the same snapshot twice");

        if (verbose) std::printf("\tremap_snapshot\n");
        remap_snapshot(Nhlo, xhlo, vhlo);

        if (verbose) std::printf("\tchoose halos\n");
        choose_halos(snap_idx, Nhlo, xhlo, vhlo);

        // update the number of halos
        Nhlo = tmp_hlo_ID.size();
        if (!Nhlo) goto snap_done;

        // callback into python to obtain the galaxy delta_x, delta_v, hlo_idx
        if (verbose) std::printf("\tcallback into python HOD\n");
        call_python_hod(snap_idx, Nhlo, hlo_idx, delta_xgal, delta_vgal);

        Ngal = hlo_idx.size();
        assert(3UL*Ngal==delta_xgal.size() && 3UL*Ngal==delta_vgal.size());
        if (!Ngal) goto snap_done;

        if ( Ngal > (std::numeric_limits<galid_t>::max()/256) )
            throw std::runtime_error("too many galaxies");

        if (verbose) std::printf("\tchoose_galaxies\n");
        choose_galaxies(snap_idx, Ngal, hlo_idx, delta_xgal, delta_vgal);

        snap_done:
            // clean up the temporary halo storage
            tmp_hlo_CHI.clear(); tmp_hlo_V.clear(); tmp_hlo_Z.clear(); tmp_hlo_ID.clear();
            if (verbose) std::printf("Done with snap index %d\n", snap_idx);
            snap_indices_done.push_back(snap_idx);
    }

    void add_snap (int snap_idx,
                   const pyb::array_t<double, pyb_arr_style> &xhlo_numpy,
                   const pyb::array_t<double, pyb_arr_style> &vhlo_numpy)
    {
        assert(xhlo_numpy.ndim()==2 && vhlo_numpy.ndim()==2);
        size_t Nhlo = xhlo_numpy.shape(0);
        assert(xhlo_numpy.shape(0)==Nhlo && vhlo_numpy.shape(0)==Nhlo);
        assert(xhlo_numpy.shape(1)==3 && vhlo_numpy.shape(1)==3);

        // copy into vectors which we will modify
        std::vector<double> xhlo(3UL*Nhlo), vhlo(3UL*Nhlo);
        std::memcpy(xhlo.data(), xhlo_numpy.data(), 3UL*Nhlo*sizeof(double));
        std::memcpy(vhlo.data(), vhlo_numpy.data(), 3UL*Nhlo*sizeof(double));

        _add_snap(snap_idx, Nhlo, xhlo, vhlo);
    }

    void _finalize ()
    {
        if (snap_indices_done.size() != Nsnaps)
            throw std::runtime_error("not all snapshots have been added");

        if (do_downsample)
        {
            // get an idea of the fiber collision rate in our sample,
            // so the subsequent downsampling is to the correct level.
            // This is justified as all this stuff is not super expensive.
            double fibcoll_rate = fibcoll</*only_measure=*/true>();

            // first downsampling before fiber collisions are applied
            if (verbose) std::printf("downsample\n");
            downsample(fibcoll_rate+Numbers::fibcoll_rate_correction);
        }

        // apply fiber collisions
        if (verbose) std::printf("fibcoll\n");
        fibcoll</*only_measure=*/false>(); 

        // now downsample to our final density
        if (do_downsample)
        {
            if (verbose) std::printf("downsample\n");
            downsample(0.0);
        }
    }

    pyb::tuple finalize ()
    {
        _finalize();

        // copy RA, DEC, Z into the numpy arrays to be used from python
        pyb::array_t<double> RAnumpy = pyb::array(RA.size(), RA.data()),
                             DECnumpy = pyb::array(DEC.size(), DEC.data()),
                             Znumpy = pyb::array(Z.size(), Z.data());
        pyb::array_t<galid_t> GALIDnumpy = pyb::array(GALID.size(), GALID.data());

        return pyb::make_tuple(RAnumpy, DECnumpy, Znumpy, GALIDnumpy);
    }

    void process_times ();
    void interpolate_chi_z ();
    void read_boss_nz ();
    void downsample (double);
    template<bool> double fibcoll ();
    void remap_snapshot (size_t,
                         std::vector<double> &,
                         std::vector<double> &) const;
    void choose_halos (int, size_t,
                       const std::vector<double> &,
                       const std::vector<double> &);
    void call_python_hod (int, size_t,
                          std::vector<hloidx_t> &,
                          std::vector<double> &,
                          std::vector<double> &) const;
    void choose_galaxies (int, size_t,
                          const std::vector<hloidx_t> &,
                          const std::vector<double> &,
                          const std::vector<double> &);
};

PYBIND11_MODULE(lc, m)
{
    pyb::class_<Mask> (m, "Mask")
        .def(pyb::init<const char *, bool>(), "boss_dir"_a, pyb::kw_only(), "veto"_a=true);

    pyb::class_<Lightcone> (m, "Lightcone")
        .def(pyb::init<const Mask&, double, double, double, const std::vector<double>&,
                       const char *,
                       double, int,
                       bool, unsigned,
                       unsigned long>(),
            "mask"_a, "Omega_m"_a, "zmin"_a, "zmax"_a, "snap_times"_a,
            pyb::kw_only(),
            "boss_dir"_a=nullptr,
            "BoxSize"_a=3e3, "remap_case"_a=0,
            "verbose"_a=false, "augment"_a=0,
            "seed"_a=137UL
        )
        .def("set_hod", &Lightcone::set_hod, "hod_fct"_a)
        .def("add_snap", &Lightcone::add_snap, "snap_idx"_a, "xhlo"_a, "vhlo"_a)
        .def("finalize", &Lightcone::finalize);
}

void Lightcone::read_boss_nz (void)
{

    std::vector<double> nz;

    char fname[512];
    // text file, each line number of objects in a redshift bin
    // we assume bins are equally spaced in redshift between zmin and zmax,
    // number of bins is inferred from the file
    std::snprintf(fname, 512, "%s/nz_DR12v5_CMASS_North_zmin%.4f_zmax%.4f.dat",
                  boss_dir, zmin, zmax);
    auto fp = std::fopen(fname, "r");
    char line[64];
    while (std::fgets(line, 64, fp))
        nz.push_back(std::atof(line));

    boss_z_hist = gsl_histogram_alloc(nz.size());
    gsl_histogram_set_ranges_uniform(boss_z_hist, zmin, zmax);
    std::memcpy(boss_z_hist->bin, nz.data(), nz.size() * sizeof(double));
}

static inline double Hz (double z, double Omega_m)
{
    // h km/s/Mpc
    return 100.0 * std::sqrt( Omega_m*gsl_pow_3(1.0+z) + (1.0-Omega_m) );
}

static double comoving_integrand (double z, void *p)
{
    double Omega_m = *(double *)p;
    return Numbers::lightspeed / Hz(z, Omega_m);
}

static void comoving (int N, double *z, double *chi, double Omega_m)
// populate chi with chi(z) in Mpc/h
{
    gsl_function F { .function=comoving_integrand, .params=&Omega_m };

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(Numbers::ws_size);

    double err;
    for (int ii=0; ii<N; ++ii)
        gsl_integration_qag(&F, 0.0, z[ii], Numbers::chi_epsabs, 0.0,
                            Numbers::ws_size, 6, ws, chi+ii, &err);

    gsl_integration_workspace_free(ws);
}

void Lightcone::process_times (void)
{
    // make sure redshifts are increasing
    for (int ii=0; ii<Nsnaps-1; ++ii) assert(snap_times[ii]>snap_times[ii+1]);

    for (int ii=0; ii<Nsnaps; ++ii) snap_redshifts.push_back(1.0/snap_times[ii] - 1.0);

    snap_chis.resize(Nsnaps);
    comoving(Nsnaps, snap_redshifts.data(), snap_chis.data(), Omega_m);

    // define the redshift boundaries
    // TODO we can maybe do something more sophisticated here
    redshift_bounds.push_back(zmin);
    for (int ii=1; ii<Nsnaps; ++ii)
        redshift_bounds.push_back(0.5*(snap_redshifts[ii-1]+snap_redshifts[ii]));
    redshift_bounds.push_back(zmax);
    
    // convert to comoving distances
    chi_bounds.resize(Nsnaps+1);
    comoving(Nsnaps+1, redshift_bounds.data(), chi_bounds.data(), Omega_m);
}

void Lightcone::interpolate_chi_z (void)
{
    double z_interp_min = zmin-0.01, z_interp_max=zmax+0.01;
    double *z_interp = (double *)std::malloc(Numbers::N_interp * sizeof(double));
    double *chi_interp = (double *)std::malloc(Numbers::N_interp * sizeof(double));
    for (int ii=0; ii<Numbers::N_interp; ++ii)
        z_interp[ii] = z_interp_min + (z_interp_max-z_interp_min)
                                      * (double)ii / (double)(Numbers::N_interp-1);
    comoving(Numbers::N_interp, z_interp, chi_interp, Omega_m);
    gsl_spline_init(z_chi_interp, chi_interp, z_interp, Numbers::N_interp);
    std::free(z_interp); std::free(chi_interp);
}

template<typename T>
static inline double per_unit (T x, double BoxSize)
{
    double x1 = (double)x / BoxSize;
    x1 = std::fmod(x1, 1.0);
    return (x1<0.0) ? x1+1.0 : x1;
}

static inline void reflect (unsigned r, double *x, double *v,
                            double BoxSize)
{
    for (unsigned ii=0; ii<3; ++ii)
        if (r & (1<<ii))
        {
            x[ii] = BoxSize - x[ii];
            v[ii] *= -1.0;
        }
}

static inline void transpose (unsigned t, double *x)
{
    if      (t==0) return;
    else if (t==1) std::swap(x[0], x[1]);
    else if (t & 1)
    {
        std::swap(x[0], x[2]);
        if (t & 4) std::swap(x[0], x[1]);
    }
    else
    {
        std::swap(x[1], x[2]);
        if (t & 4) std::swap(x[0], x[1]);
    }
}

void Lightcone::remap_snapshot (size_t Nhlo,
                                std::vector<double> &xhlo,
                                std::vector<double> &vhlo) const
{
    std::vector<double> tmp_xhlo(3UL*Nhlo), tmp_vhlo(3UL*Nhlo);

    unsigned r = augment / 6; // labels the 8 reflection cases
    unsigned t = augment % 6; // labels the 6 transposition cases

    #pragma omp parallel for schedule (dynamic, 1024)
    for (size_t ii=0; ii<Nhlo; ++ii)
    {
        double x[3], v[3];
        for (int kk=0; kk<3; ++kk)
        {
            x[kk] = xhlo[3*ii+kk];
            v[kk] = vhlo[3*ii+kk];
        }

        reflect(r, x, v, BoxSize);

        transpose(t, x);
        C.Transform(per_unit(x[0], BoxSize), per_unit(x[1], BoxSize), per_unit(x[2], BoxSize),
                    tmp_xhlo[3*ii+0], tmp_xhlo[3*ii+1], tmp_xhlo[3*ii+2]);
        for (int kk=0; kk<3; ++kk) tmp_xhlo[3*ii+kk] *= BoxSize;

        transpose(t, v);
        C.VelocityTransform(v[0], v[1], v[2],
                            tmp_vhlo[3*ii+0], tmp_vhlo[3*ii+1], tmp_vhlo[3*ii+2]);
    }

    xhlo = std::move(tmp_xhlo);
    vhlo = std::move(tmp_vhlo);
}

static inline galid_t make_galid (galid_t snap_idx, galid_t gal_idx)
{
    return ( snap_idx << ((sizeof(galid_t)-1)*CHAR_BIT) ) + gal_idx;
}

void Lightcone::choose_halos (int snap_idx, size_t Nhlo,
                              const std::vector<double> &xhlo,
                              const std::vector<double> &vhlo)
{
    assert(tmp_hlo_CHI.empty() && tmp_hlo_V.empty() && tmp_hlo_Z.empty() && tmp_hlo_ID.empty());

    // h km/s/Mpc
    const double H_ = Hz(snap_redshifts[snap_idx], Omega_m);

    #pragma omp parallel
    {
        // the accelerator is non-const upon evaluation
        auto acc = z_chi_interp_acc[omp_get_thread_num()];
        double los[3];

        #pragma omp for schedule (dynamic, 1024)
        for (size_t jj=0; jj<Nhlo; ++jj)
        {
            for (int kk=0; kk<3; ++kk) los[kk] = xhlo[3*jj+kk] - Geometry::origin[kk]*BoxSize*Li[kk];
            double chi = std::hypot(los[0], los[1], los[2]);

            double vhloproj = (los[0]*vhlo[3*jj+0]+los[1]*vhlo[3*jj+1]+los[2]*vhlo[3*jj+2])
                              / chi;

            // map onto the lightcone -- we assume that this is a relatively small correction
            //                           so we just do it to first order
            double delta_z = (chi - snap_chis[snap_idx])
                             / (Numbers::lightspeed + vhloproj) * H_;

            // correct the position accordingly
            for (int kk=0; kk<3; ++kk) los[kk] -= delta_z * vhlo[3*jj+kk] / H_;

            // only a small correction I assume (and it is borne out by experiment!)
            chi = std::hypot(los[0], los[1], los[2]);

            if (chi>chi_bounds[snap_idx] && chi<chi_bounds[snap_idx+1])
            // we are in the comoving shell that's coming from this snapshot
            // since we're doing halos here, no testing for mask
            {
                double z = gsl_spline_eval(z_chi_interp, chi, acc);

                #pragma omp critical (HLO_CHOOSE_APPEND)
                {
                    for (int kk=0; kk<3; ++kk) tmp_hlo_CHI.push_back(los[kk]);
                    for (int kk=0; kk<3; ++kk) tmp_hlo_V.push_back(vhlo[3*jj+kk]);
                    tmp_hlo_Z.push_back(z);
                    tmp_hlo_ID.push_back(jj);
                }
            }
        }
    }
}

void Lightcone::choose_galaxies (int snap_idx, size_t Ngal,
                                 const std::vector<hloidx_t> &hlo_idx,
                                 const std::vector<double> &delta_xgal,
                                 const std::vector<double> &delta_vgal)
{
    // h km/s/Mpc
    const double H_ = Hz(snap_redshifts[snap_idx], Omega_m);

    // (Mpc/h)/(km/s)
    const double rsd_factor = (1.0+snap_redshifts[snap_idx]) / H_;

    int mangle_status = 1;

    #pragma omp parallel reduction(*:mangle_status)
    {
        // the accelerator is non-const upon evaluation
        auto acc = z_chi_interp_acc[omp_get_thread_num()];
        double los[3], vgal[3];
        int mangle_status_ = 1;

        #pragma omp for schedule (dynamic, 1024)
        for (size_t jj=0; jj<Ngal; ++jj)
        {
            if (!mangle_status_) [[unlikely]] continue;

            assert(hlo_idx[jj]<tmp_hlo_ID.size());

            for (int kk=0; kk<3; ++kk) los[kk] = tmp_hlo_CHI[3*hlo_idx[jj]+kk] + delta_xgal[3*jj+kk];
            double chi = std::hypot(los[0], los[1], los[2]);

            for (int kk=0; kk<3; ++kk) vgal[kk] = tmp_hlo_V[3*hlo_idx[jj]+kk] + delta_vgal[3*jj+kk];

            // the radial velocity
            double vproj = (los[0]*vgal[0]+los[1]*vgal[1]+los[2]*vgal[2]) / chi;

            // now add RSD (the order is important here, as RSD can be a pretty severe effect)
            for (int kk=0; kk<3; ++kk) los[kk] += rsd_factor * vproj * los[kk] / chi;

            // I'm lazy
            chi = std::hypot(los[0], los[1], los[2]);

            // discard galaxies falling out at the edges
            if (chi>chi_bounds[0] && chi<chi_bounds[Nsnaps]) [[likely]]
            {
                // rotate the line of sight into the NGC footprint and transpose the axes into
                // canonical order
                double x1, x2, x3;
                x1 = std::cos(Geometry::alpha) * los[2] - std::sin(Geometry::alpha) * los[1];
                x2 = los[0];
                x3 = std::sin(Geometry::alpha) * los[2] + std::cos(Geometry::alpha) * los[1];

                double z = gsl_spline_eval(z_chi_interp, chi, acc);
                double theta = std::acos(x3/chi);
                double phi = std::atan2(x2, x1);

                double dec = 90.0-theta/M_PI*180.0;
                double ra = phi/M_PI*180.0;
                if (ra<0.0) ra += 360.0;
                ra += Geometry::beta;

                galid_t galid;

                // for the angular mask
                cmangle::Point pt;
                cmangle::point_set_from_radec(&pt, ra, dec);
                bool m = mask.masked(pt, mangle_status_);
                mangle_status *= mangle_status_;
                if (!mangle_status_) [[unlikely]] continue;
                if (m) goto not_chosen;

                // compute the galaxy ID
                galid = make_galid(snap_idx, jj);

                // this is executed in parallel, modifying global variables
                #pragma omp critical (CHOOSE_APPEND)
                {
                    RA.push_back(ra);
                    DEC.push_back(dec);
                    Z.push_back(z);
                    GALID.push_back(galid);
                }
                continue;

                not_chosen : /* do nothing */;
            }
        }
    }

    if (!mangle_status) throw std::runtime_error("mangle failed");
}

void Lightcone::call_python_hod (int snap_idx, size_t Nhlo,
                                 std::vector<hloidx_t> &hlo_idx,
                                 std::vector<double> &delta_xgal,
                                 std::vector<double> &delta_vgal) const
{
    assert( hlo_idx.empty() && delta_xgal.empty() && delta_vgal.empty() );

    pyb::array_t<double> hlo_Z_numpy = pyb::array(Nhlo, tmp_hlo_Z.data());
    pyb::array_t<hloidx_t> hlo_ID_numpy = pyb::array(Nhlo, tmp_hlo_ID.data());

    // here we give control back to python
    pyb::tuple res = hod_fct(snap_idx, hlo_ID_numpy, hlo_Z_numpy);
    assert(res.size() == 3);

    pyb::array_t<hloidx_t> hlo_idx_numpy = res[0].cast<pyb::array_t<hloidx_t, pyb_arr_style>>();
    pyb::array_t<double> delta_xgal_numpy = res[1].cast<pyb::array_t<double, pyb_arr_style>>();
    pyb::array_t<double> delta_vgal_numpy = res[2].cast<pyb::array_t<double, pyb_arr_style>>();

    // make sure what we got makes sense
    assert(hlo_idx_numpy.ndim()==1 && delta_xgal_numpy.ndim()==2 && delta_vgal_numpy.ndim()==2);
    size_t Ngal = hlo_idx_numpy.shape(0);
    assert(delta_xgal_numpy.shape(0)==Ngal && delta_xgal_numpy.shape(1)==3);
    assert(delta_vgal_numpy.shape(0)==Ngal && delta_vgal_numpy.shape(1)==3);

    // copy into our nice internal arrays
    hlo_idx.resize(Ngal); delta_xgal.resize(3UL*Ngal); delta_vgal.resize(3UL*Ngal);
    std::memcpy(hlo_idx.data(), hlo_idx_numpy.data(), Ngal*sizeof(hloidx_t));
    std::memcpy(delta_xgal.data(), delta_xgal_numpy.data(), 3UL*Ngal*sizeof(double));
    std::memcpy(delta_vgal.data(), delta_vgal_numpy.data(), 3UL*Ngal*sizeof(double));
}

void Lightcone::downsample (double plus_factor)
// downsamples to the boss_z_hist, up to plus_factor
{
    const gsl_histogram *target_hist = boss_z_hist;

    gsl_histogram *sim_z_hist = gsl_histogram_clone(target_hist); // get the same ranges
    gsl_histogram_reset(sim_z_hist);
    for (auto z : Z)
        gsl_histogram_increment(sim_z_hist, z);
    double keep_fraction[target_hist->n];
    for (size_t ii=0; ii<target_hist->n; ++ii)
    {
        // note the GSL bin counts are doubles already
        keep_fraction[ii] = (1.0+plus_factor)
                            * gsl_histogram_get(target_hist, ii)
                            / gsl_histogram_get(sim_z_hist, ii);

        // give warning if something goes wrong
        if (keep_fraction[ii] > 1.02 && plus_factor==0.0)
            std::fprintf(stderr, "[WARNING] at z=%.4f, do not have enough galaxies, keep_fraction=%.4f\n",
                                 0.5*(target_hist->range[ii]+target_hist->range[ii+1]),
                                 keep_fraction[ii]);
    }

    // for debugging
    if (verbose)
    {
        std::printf("keep_fraction:\n");
        for (int ii=0; ii<target_hist->n; ++ii)
            std::printf("%.2f ", keep_fraction[ii]);
        std::printf("\n");
    }
    
    std::vector<double> ra_tmp, dec_tmp, z_tmp;
    std::vector<galid_t> galid_tmp;

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    for (size_t ii=0; ii<RA.size(); ++ii)
    {
        size_t idx;
        gsl_histogram_find(sim_z_hist, Z[ii], &idx);
        if (gsl_rng_uniform(rng) < keep_fraction[idx])
        {
            ra_tmp.push_back(RA[ii]);
            dec_tmp.push_back(DEC[ii]);
            z_tmp.push_back(Z[ii]);
            galid_tmp.push_back(GALID[ii]);
        }
    }

    // clean up
    gsl_histogram_free(sim_z_hist);
    gsl_rng_free(rng);

    // assign output
    RA = std::move(ra_tmp);
    DEC = std::move(dec_tmp);
    Z = std::move(z_tmp);
    GALID = std::move(galid_tmp);
}

struct GalHelper
{
    // we store the coordinates so we have more efficient lookup hopefully
    // (at the cost of more memory but should be fine)
    double ra, dec, z;
    galid_t galid;
    int64_t hp_idx;
    unsigned long id;
    pointing ang;

    GalHelper ( ) = default;

    GalHelper (unsigned long id_, double ra_, double dec_, double z_, galid_t galid_,
               const T_Healpix_Base<int64_t> &hp_base) :
        ra {ra_}, dec {dec_}, z{z_}, id {id_}, galid {galid_}
    {
        double theta = (90.0-dec) * M_PI/180.0;
        double phi = ra * M_PI/180.0;
        ang = pointing(theta, phi);
        hp_idx = hp_base.ang2pix(ang);
    };
};

static inline long double hav (long double theta)
{
    auto x = std::sin(0.5*theta);
    return x * x;
}

static inline long double haversine (const pointing &a1, const pointing &a2)
{
    long double t1=a1.theta, t2=a2.theta, p1=a1.phi, p2=a2.phi;
    return hav(t1-t2) + hav(p1-p2) * ( - hav(t1-t2) + hav(t1+t2) );
}

template<bool only_measure>
double Lightcone::fibcoll ()
{
    // for sampling from overlapping regions according to collision rate
    // and assignment of random galaxy IDs
    std::vector<gsl_rng *> rngs;
    for (int ii=0; ii<omp_get_max_threads(); ++ii)
    {
        rngs.push_back(gsl_rng_alloc(gsl_rng_default));
        gsl_rng_set(rngs.back(), seed+(unsigned long)ii);
    }

    // for HOPEFULLY efficient nearest neighbour search
    static const int64_t nside = 128; // can be tuned, with 128 pixarea/discarea~225
    const T_Healpix_Base hp_base (nside, RING, SET_NSIDE);

    // we assign random 64bit IDs to the galaxies
    // This is better than using their index in the input arrays, since these arrays
    // have some ordering (low redshifts go first), and this would slightly bias
    // the fiber collision removal implmented later
    // The collision rate is tiny for 64bit random numbers and thus has no impact.
    std::vector<GalHelper> all_vec { RA.size() };
    #pragma omp parallel
    {
        auto rng = rngs[omp_get_thread_num()];
        
        #pragma omp for
        for (size_t ii=0; ii<RA.size(); ++ii)
            all_vec[ii] = GalHelper(gsl_rng_get(rng), RA[ii], DEC[ii], Z[ii], GALID[ii], hp_base);
    }

    if constexpr (!only_measure)
    {
        // empty the output arrays
        RA.clear();
        DEC.clear();
        Z.clear();
        GALID.clear();
    }

    // sort these for efficient access with increasing healpix index
    std::sort(all_vec.begin(), all_vec.end(),
              [](const GalHelper &a, const GalHelper &b){ return a.hp_idx<b.hp_idx; });

    // again, for efficient access we compute the ranges
    // As we only cover a small sky area, a map is probably more efficient than an array
    std::map<int64_t, std::pair<size_t, size_t>> ranges;
    int64_t current_hp_idx = all_vec[0].hp_idx;
    size_t range_start = 0;
    for (size_t ii=0; ii<all_vec.size(); ++ii)
    {
        const auto &g = all_vec[ii];
        if (g.hp_idx != current_hp_idx)
        {
            ranges[current_hp_idx] = std::pair<size_t, size_t>(range_start, ii);
            range_start = ii;
            current_hp_idx = g.hp_idx;
        }
    }

    // we assume the inputs have sufficient entropy that always removing the first galaxy is fine

    size_t Nkept = 0;

    #pragma omp parallel reduction(+:Nkept)
    {
        // allocate for efficiency, I don't know how rangeset works internally
        // so we do this the old-fashioned way (not private etc)
        rangeset<int64_t> query_result;
        std::vector<int64_t> query_vector;

        // generator used by this thread
        auto rng = rngs[omp_get_thread_num()];

        #pragma omp for schedule (dynamic, 1024)
        for (size_t ii=0; ii<all_vec.size(); ++ii)
        {
            const auto &g = all_vec[ii];

            // one can gain performance here by playing with "fact"
            hp_base.query_disc_inclusive(g.ang, Numbers::angscale, query_result, /*fact=*/4);
            query_result.toVector(query_vector);

            for (auto hp_idx : query_vector)
            {
                const auto this_range_ptr = ranges.find(hp_idx);
                if (this_range_ptr == ranges.end())
                    // not found
                    continue;
                const auto this_range = this_range_ptr->second;
                for (size_t ii=this_range.first; ii<this_range.second; ++ii)
                    if (g.id > all_vec[ii].id
                        && haversine(g.ang, all_vec[ii].ang) < hav(Numbers::angscale)) [[unlikely]]
                    // use the fact that haversine is monotonic to avoid inverse operation
                    // by using the greater-than check, we ensure to remove only one member
                    // of each pair
                        goto collided;
            }
            goto not_collided;

            // TODO it could also make sense to call the rng every time we have a collision
            //      in the above loop instead. Maybe not super important though.
            collided :
            if (gsl_rng_uniform(rng)<Numbers::collrate) continue;

            not_collided :
            ++Nkept;
            if constexpr (!only_measure)
                // no collision, let's keep this galaxy. We're writing into global variables so need
                // to be careful
                #pragma omp critical (FIBCOLL_APPEND)
                {
                    RA.push_back(g.ra);
                    DEC.push_back(g.dec);
                    Z.push_back(g.z);
                    GALID.push_back(g.galid);
                }
        }
    }

    // clean up
    for (int ii=0; ii<rngs.size(); ++ii)
        gsl_rng_free(rngs[ii]);
    
    if constexpr (!only_measure)
        assert(Nkept == RA.size());

    double fibcoll_rate = (double)(all_vec.size()-Nkept)/(double)(all_vec.size());

    return fibcoll_rate;
}
