Code adapted from [here](https://github.com/leanderthiele/nuvoid_production).

The Makefile should hopefully work out of the box.
External dependences are GSL and pybind11 (I have copied Martin White's cuboid remapping code,
pymangle, and a minimal healpix into this repository).

Basic usage example please see `example.py`.
The heavy sections of the code are multi-threaded, control via `OMP_NUM_THREADS`.

For efficiency, I have factored out the mask into a separate class.
This way, an instance of this class only needs to be constructed once at the beginning
and can then be passed to the lightcone generator repeatedly.

## mask constructor

**Mandatory**
* `boss_dir`: a directory containing all `.ply` files for the angular mask.
  If `veto=False`, only `mask_DR12v5_CMASS_North.ply` is required.

**Optional**
* `veto = True`: whether the veto masks are to be applied, I see no reason
  why not to.


## lightcone Constructor

**Mandatory**
* `mask`: instance of the `Mask` class, constructed as described above
* `Omega_m`
* `zmin`, `zmax`: the redshift boundaries
* `snap_times`: a monotonically *decreasing* list of scale factors
  (i.e., increasing redshift)

**Optional**
* `boss_dir`: a directory containing a text file named
  `nz_DR12v5_CMASS_North_zmin%.4f_zmax%.4f.dat`,
  where the placeholders are filled by the `zmin` and `zmax` arguments.
  This file defines the desired redshift distribution.
  It is assumed that the histogramming is in uniform redshift bins
  between `zmin` and `zmax`.
  Each line of the file should contain a single integer for the number
  of galaxies in the corresponding bin.
  Set this to `None` to skip downsampling to any redshift distribution.
* `BoxSize = 3e3`: Mpc/h
* `remap_case = 0`: either 0 or 1
* `verbose = False`
* `augment = 0`: between 0 and 47. There are 96 possible augmentations implemented,
  namely the product of 2 remaps (above argument), 8 reflections, 6 transpositions.
  Please see the function `Lightcone::remap_snapshot` for the exact definition if it
  is important.
* `seed = 137`: randomness is introduced in downsampling to n(z) and fiber collisions.
  I believe determinism will only be possible when running on a single thread.

## setting the HOD

It is somehow hard to do this in the constructor, some pybind11 issue.
Hacky solution is just a separate method, `set_hod` which just takes a python callable.
I have some example in `example.py`, here I copy the documentation:
```python
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
```



## adding galaxies

Using the `add_snap` method. This can be called in any order of snapshots.

**Arguments**
* `snap_idx`: integer index into the `snap_times` array
* `xgal`: Nx3 numpy array, galaxy positions in Mpc/h
* `vgal`: Nx3 numpy array, physical km/s.
* `vhlo`: Nx3 numpy array, physical km/s. In principle this can be empty if `correct=False`,
          but I haven't tried this.


## getting result

Using the `finalize` method. Call after calling `add_snap` for all snapshots.
Returns a tuple RA, DEC, z, galid where the angles are in degrees.
These are already rotated into the NGC footprint.

The galid return value is a unique ID for each galaxy. This is an unsigned integer (currently 32bit),
whose leading byte is the snapshot index this galaxy came from, and the remaining lower bytes are the
index into the original galaxy list passed to `add_snap`.
The `example.py` file contains the following function to recover this information:
```python
def split_galid (gid) :
    # returns snapshot index, galaxy index
   return np.divmod(gid, 2**((gid.itemsize-1)*8))
```
