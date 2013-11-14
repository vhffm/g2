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
