"""
Helper Functions of Chaos Characterization.
"""

import numpy as np
from array_helpers import intersect2d
from g2_helpers import twopi

def co_intersection(cpid, opid, cm, om):
    """
    Current/Old Intersection.

    Compute Intersection of Indices and Masses.

    @params
    cpid = array of current particle IDs
    opid =          old
    cm   = array of current particle masses
    co   =          old

    @returns
    subset of particle IDs where IDs and masses match
    """

    current = np.transpose(np.array([cpid, cm])).copy()
    old = np.transpose(np.array([opid, om])).copy()
    return intersect2d(current, old)[:,0].astype(int)

def compute_lyapunov(ds, istep0, nsteps, tout):
    """
    Compute Lyapunov Exponents.

    @params
    ds = [0:nsteps,0:npartmax]      Particle separation
    nstep0 = [0,npartmax]           First step where separation occured
    t = [0,nsteps]                  Time array

    @returns
    lce = [0:nsteps,0:npartmax]
    """

    lce = np.ones_like(ds) * np.nan
    for istep, nstep in enumerate(nsteps):
        for ipart in range(ds.shape[1]):
            # Only compute if we are above t where first separation occured
            if istep > istep0[ipart]:
                # Compute lambda wrt to last step
                lce[istep,ipart] = np.log(ds[istep,ipart] / ds[istep0[ipart],ipart]) / \
                                   ( tout[istep] - tout[istep0[ipart]] )
                # Normalization foo?
                # ...

    return lce

def bin_lyapunov(a0, lce):
    """
    Create bins in semi-major axis.
    Then sort Lyapunov Characteristics Exponent (LCE) into bins.
    Per bin, compute median, mean, and standard deviation of LCE.
    Save.

    @parameters - lce - LCE(a0,t)
                - a0  - Initial semi-major axes of particles
    @return     - lce_med       - Mean value per bin ~ lce(a0i,t)
                - lce_median    - Median per bin
                - lce_std       - Standard deviation per bin
                - a0_bin_mids   - Bin midpoints
    """

    # Set up bins
    nbins = 32; a0_lo = 0.0; a0_hi = 5.5

    # Create bins
    a0_bin_edges = np.linspace(a0_lo, a0_hi, nbins+1)           # len=nbins+1
    a0_bin_mids = (a0_bin_edges[1:] + a0_bin_edges[:-1]) / 2.0  # len=nbins

    # Digitize
    # Cf. http://docs.scipy.org/doc/numpy/reference/generated/numpy.digitize.html#numpy.digitize
    digitalism = np.digitize(a0, a0_bin_edges)
    digitalism -= 1

    #   - Some Error Checking
    if len(a0_bin_edges)-1 in digitalism:
        raise Exception("Particle Outside Binning Range (a0 > a_hi). Adjust Limits!")

    # Allocate target arrays
    lce_mean = np.zeros([lce.shape[0], a0_bin_mids.shape[0]])
    lce_median = np.zeros_like(lce_mean)
    lce_std = np.zeros_like(lce_mean)

    # Sweep arrays first by time, then by particle
    # Compute mean/median/std and write back
    for itslice in range(0,lce.shape[0]):
        lce_binned = [[] for _ in range(0,nbins)]
        for ipslice in range(0,lce.shape[1]):
            lce_binned[digitalism[ipslice]].append(lce[itslice,ipslice])
        for ibin in range(0,nbins):
            if len(lce_binned[ibin]) == 0:
                lce_mean[itslice,ibin] = np.nan
                lce_median[itslice,ibin]= np.nan
                lce_std[itslice,ibin] = np.nan
            else:
                lce_mean[itslice,ibin] = np.mean(lce_binned[ibin])
                lce_median[itslice,ibin] = np.median(lce_binned[ibin])
                lce_std[itslice,ibin] = np.std(lce_binned[ibin])

    # Return values
    return lce_mean, lce_median, lce_std, a0_bin_mids
