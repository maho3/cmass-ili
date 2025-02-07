import os
from os.path import join
import numpy as np
import h5py
from astropy.stats import scott_bin_width
from scipy.interpolate import InterpolatedUnivariateSpline

from ..utils import timing_decorator
from ..survey.tools import gen_randoms


def get_nofz(z, fsky, cosmo=None):
    ''' calculate nbar(z) given redshift values and f_sky (sky coverage
    fraction)
    Parameters
    ----------
    z : array like
        array of redshift values
    fsky : float
        sky coverage fraction
    cosmo : cosmology object
        cosmology to calculate comoving volume of redshift bins
    Returns
    -------
    number density at input redshifts: nbar(z)
    Notes
    -----
    * based on nbodykit implementation
    * Deprecated with pypower implementation
    '''
    # calculate nbar(z) for each galaxy
    _, edges = scott_bin_width(z, return_bins=True)

    dig = np.searchsorted(edges, z, "right")
    N = np.bincount(dig, minlength=len(edges)+1)[1:-1]

    R_hi = cosmo.comoving_distance(edges[1:])  # Mpc/h
    R_lo = cosmo.comoving_distance(edges[:-1])  # Mpc/h

    dV = (4./3.) * np.pi * (R_hi**3 - R_lo**3) * fsky

    nofz = InterpolatedUnivariateSpline(
        0.5*(edges[1:] + edges[:-1]), N/dV, ext='const')

    return nofz


@timing_decorator
def load_lightcone(
        source_dir, hod_seed=0, aug_seed=0,
        filter_name=None, is_North=True):
    pfx = 'ngc' if is_North else 'sgc'
    if filter_name is None:
        infile = join(
            source_dir, f'{pfx}_lightcone',
            f'hod{hod_seed:05}_aug{aug_seed:05}.h5')
    else:
        infile = join(
            source_dir, f'{pfx}_filtered',
            f'hod{hod_seed:05}_aug{aug_seed:05}_{filter_name}.h5')

    with h5py.File(infile, 'r') as f:
        ra = f['ra'][...]
        dec = f['dec'][...]
        z = f['z'][...]
        rdz = np.stack([ra, dec, z], axis=-1)

        if 'weight' in f:
            weight = f['weight'][...]
        else:
            weight = np.ones(len(rdz))

    return rdz, weight


def save_summary(outpath, name, **kwargs):
    with h5py.File(outpath, 'a') as f:
        group = f.create_group(name)
        for key, value in kwargs.items():
            group.create_dataset(key, data=value)


@timing_decorator
def load_randoms(wdir):
    path = join(wdir, 'obs', 'random0_DR12v5_CMASS_North_PRECOMPUTED.npy')
    if os.path.exists(path):
        return np.load(path)
    randoms = gen_randoms()
    np.save(path, randoms)
    return randoms
