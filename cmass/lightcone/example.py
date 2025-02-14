import numpy as np
from matplotlib import pyplot as plt
import lc

def split_galid (gid) :
    return np.divmod(gid, 2**((gid.itemsize-1)*8))

# how many galaxies we randomly generate per volume
N = 300000

# where mask .ply files and redshift histogram are stored
boss_dir = './testdata'

# scale factors, need to be monotonically decreasing
snap_times = [0.7, 0.65, 0.6, ]

# some randomness
rng = np.random.default_rng(137)

# make the mask (this is a bit expensive, so we do it outside the lightcone so it can be re-used)
print('Starting mask loading...')
m = lc.Mask(boss_dir=boss_dir, veto=False)
print('...finished mask loading')

def hod_fct (
        snap_idx: int,
        hlo_idx: np.ndarray[np.uint64],
        z: np.ndarray[np.float64]) -> tuple :
    """ This is the callback function for the HOD

    The arguments are
    - snap_idx: an integer going into the snapshots list
    - hlo_idx: an array of integer indices going into the halo arrays
    - z: an array of doubles corresponding to the redshifts at which halos intersect lightcone

    This function is expected to return the following:
    - hlo_idx: a [N] array of integer indices going into the arrays *passed as inputs to this function*
               So these *do not* correspond to the halo indices that came as inputs!
    - delta_x: a [N,3] array of galaxy position offsets from their host halo centers,
               comoving Mpc/h
    - delta_v: a [N,3] array of galaxy velocity offsets from their host halo velocities,
               physical km/s

    NOTE that this can be any callable as usual,
         so you can use any object as long as it implements the __call__ method
    """

    # as an example, only halos with certain redshifts get galaxies. Each gets a pair of galaxies.
    # I set the delta_x and delta_v very small so it should be possible to see these pairs in the output
    hlo_idx_out = np.arange(0, len(hlo_idx), dtype=np.uint64)
    select = np.fabs(np.sin(100.0*z))<0.2
    hlo_idx_out = np.repeat(hlo_idx_out[select], 2)
    delta_x = (rng.random((len(hlo_idx_out), 3))-0.5) * 0.01
    delta_v = (rng.random((len(hlo_idx_out), 3))-0.5) * 1.0

    return hlo_idx_out, delta_x, delta_v


# run the constructor
l = lc.Lightcone(
    mask=m,
    Omega_m=0.3, zmin=0.40, zmax=0.65, snap_times=snap_times,
    boss_dir=None, # NOTE setting to None here disables n(z) downsampling now
    verbose=True
)

# set the HOD function
l.set_hod(hod_fct)

# add some snapshots
for snap_idx, a in enumerate(snap_times):
    xhlo = rng.random((N, 3)) * 3e3
    vhlo = (rng.random((N, 3))-0.5) * 300
    l.add_snap(snap_idx, xhlo, vhlo)

ra, dec, z, galid = l.finalize()

galsnap, galidx = split_galid(galid)
print(f'{galsnap.min()} <= galsnap <= {galsnap.max()}')
print(f'{galidx.min()} <= galidx <= {galidx.max()}')

# check if mask is working correctly
fig, ax = plt.subplots(figsize=(20, 10), ncols=2)
ax_mask = ax[0]
ax_snap = ax[1]
choose = (z > 0.5) * (z < 0.6)

ax_mask.plot(ra[choose], dec[choose], linestyle='none', marker='o', markersize=1)
ax_mask.set_xlabel('RA [deg]')
ax_mask.set_ylabel('DEC [deg]')

ax_snap.plot(galsnap, z, linestyle='none', marker='o', markersize=1)
ax_snap.set_xlabel('snapshot index')
ax_snap.set_ylabel('redshift')
plt.show()
