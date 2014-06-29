"""
Read Metric Distance. Compute Different Fits. Note Easter Egg.

1) First Exponential
2) Second Exponential
3) Final Power Law
4) Time To Hill Radius

Computing Fits For:
1) Full Disk (Semi-Major Axis)
2)     a<1.0
3) 1.0<a<2.0
4) 2.0<a<3.0
5)     a>3.0
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit
import scipy.stats.mstats as mstats
import constants as C
import kepler_helpers as kh

# Define Fitting Functions
def fit_func_exp(t, a, b, c):
    return (a*np.exp(b*t))+c

def fit_func_pow(t, a, b, c):
    return (a*t**b)+c

# Other Functions
def orbital_distance(dx):
    GM = 1.0
    x1 = 1.0; y1 = 0.0; z1 = 0.0
    x2 = 1.0 + dx; y2 = 0.0; z2 = 0.0
    v1 = np.sqrt(GM/x1)
    v2 = np.sqrt(GM/x2)
    vx1 = 0.0; vy1 = v1; vz1 = 0.0
    vx2 = 0.0; vy2 = v2; vz2 = 0.0
    return np.sqrt(kh.cart2metricX(x1, y1, z1, vx1, vy1, vz1, \
                                   x2, y2, z2, vx2, vy2, vz2))

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='Distances.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='Fits.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Read
print "// Reading %s" % args.infile
npz = np.load(args.infile)
tout = npz["tout"]
a1 = npz["a1"]
rho2 = npz["rho2"]
# ttx = np.repeat([tout], rho2.shape[1], axis=0).T
ttx = tout
rho = np.sqrt(rho2/rho2[0,:])
rho = mstats.gmean(rho, axis=1)

# #############################################################################
# Semi-Major Axis Splitting
# #############################################################################
print "// Splitting In Semi-Major Axis"

# Loop Bins
# Check If First A1 In Bin
# If So, Add To List

rho_abins = []; ratio_gmean_abins = []
for ii in [ 0, 1, 2, 3 ]:

    # Memcpy
    rho_loc = np.sqrt(npz["rho2"]).copy()
    ratio_loc = rho_loc/rho_loc[0,:]
    
    # Set Filters
    if ii == 0:
        mybool = a1<1.0
    elif ii == 1:
        mybool = np.logical_and(a1>=1.0, a1<2.0) 
    elif ii == 2:
        mybool = np.logical_and(a1>=2.0, a1<3.0) 
    elif ii == 3:
        mybool = a1>=3.0
        
    # Loop Particles
    first = True
    for iparticle in range(rho_loc.shape[1]):
        # Is First Step In This Bin?
        if mybool[0,iparticle]:
            if first:
                rho_loc_abin = rho_loc[:,iparticle]
                ratio_loc_abin = ratio_loc[:,iparticle]
                first = False
            else:
                rho_loc_abin = np.vstack((rho_loc_abin, \
                                          rho_loc[:,iparticle]))
                ratio_loc_abin = np.vstack((ratio_loc_abin, \
                                            ratio_loc[:,iparticle]))

    # Compute Geometric Mean
    if not first:
        ratio_gmean_abin = mstats.gmean(np.atleast_2d(ratio_loc_abin), axis=0)
    else:
        ratio_gmean_abin = np.zeros_like(tout) * np.nan
        rho_loc_abin = np.atleast_2d(np.zeros_like(tout)) * np.nan 
    
    # Append
    rho_abins.append(np.atleast_2d(rho_loc_abin))
    ratio_gmean_abins.append(ratio_gmean_abin)

# #############################################################################
# Fit 01
# #############################################################################
print "// Fitting First Exponential (Global)"
ttx_loc = ttx[rho<1.0e6]
rho_loc = rho[rho<1.0e6]
try:
    fitParams01, _ = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                               p0 = [ 40000.0, 0.01, -10000.0 ])
    te01 = 1.0/fitParams01[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    te01 = np.nan
print "   E-Folding Time = %.2e yr" % te01

print "// Fitting First Exponential (Semi-Major Axis Bin)"
te01_abins = []
for ii in [ 0, 1, 2, 3 ]:
    ttx_loc = tout[ratio_gmean_abins[ii]<1.0e6]
    ratio_loc = ratio_gmean_abins[ii][ratio_gmean_abins[ii]<1.0e6]
    try:
        fitParams01, _ = curve_fit(fit_func_exp, ttx_loc, ratio_loc, \
                                   p0 = [ 40000.0, 0.01, -10000.0 ])
        te01_abins.append(1.0/fitParams01[1])
    except (RuntimeError, TypeError):
        "!! Fit Failed"
        te01_abins.append(np.nan)
    print "   E-Folding Time = %.2e yr" % te01_abins[-1]

# #############################################################################
# Fit 02
# #############################################################################
print "// Fitting Second Exponential (Global)"
bool_01 = rho>1.0e7
bool_02 = rho<=1.0e12
bool_cb = np.logical_and(bool_01, bool_02)
ttx_loc = ttx[bool_cb]
rho_loc = rho[bool_cb]
try:
    fitParams02, _ = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                               p0 = [ 40000.0, 0.01, -10000.0 ])
    te02 = 1.0/fitParams02[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    te02 = np.nan
print "   E-Folding Time = %.2e yr" % te02

print "// Fitting Second Exponential (Semi-Major Axis Bin)"
te02_abins = []
for ii in [ 0, 1, 2, 3 ]:
    bool_01 = ratio_gmean_abins[ii]>1.0e7
    bool_02 = ratio_gmean_abins[ii]<=1.0e12
    bool_cb = np.logical_and(bool_01, bool_02)

    ttx_loc = tout[bool_cb]
    ratio_loc = ratio_gmean_abins[ii][bool_cb]
    try:
        fitParams02, _ = curve_fit(fit_func_exp, ttx_loc, ratio_loc, \
                                   p0 = [ 40000.0, 0.01, -10000.0 ])
        te02_abins.append(1.0/fitParams02[1])
    except (RuntimeError, TypeError):
        "!! Fit Failed"
        te02_abins.append(np.nan)
    print "   E-Folding Time = %.2e yr" % te02_abins[-1]

# #############################################################################
# Fit 03
# #############################################################################
print "// Fitting Power Law (Global)"
ttx_loc = ttx[rho>1.0e12]
rho_loc = rho[rho>1.0e12]
try:
    fitParams03, _ = curve_fit(fit_func_pow, ttx_loc, rho_loc, \
                               p0 = [ 5.0e12, 0.1, -1.0e13 ])
    slope03 = fitParams03[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    slope03 = np.nan
print "   Slope = %.2e" % slope03

print "// Fitting Power Law (Semi-Major Axis Bin)"
slope03_abins = []
for ii in [ 0, 1, 2, 3 ]:
    ttx_loc = tout[ratio_gmean_abins[ii]>2.0e12]
    ratio_loc = ratio_gmean_abins[ii][ratio_gmean_abins[ii]>2.0e12]
    try:
        fitParams03, _ = curve_fit(fit_func_pow, ttx_loc, ratio_loc, \
                               p0 = [ 5.0e12, 0.1, -1.0e13 ])
        slope03_abins.append(fitParams03[1])
    except (RuntimeError, TypeError):
        "!! Fit Failed"
        slope03_abins.append(np.nan)
    print "   Slope = %.2e" % slope03_abins[-1]

# #############################################################################
# Time To Hill Radius
# #############################################################################
m = 5.0 * C.mearth / rho2.shape[1]
Rh = 1.0 * (m / 3.0 / C.msun)**(1.0/3.0)
rho_hill = orbital_distance(Rh)
print "   (N, m, R_hill) = (%i, %.2e kg, %.2e AU)" % (rho2.shape[1], m, Rh)

print "// Computing Time To Hill Radius (Global)"
thill = []
for iparticle in range(rho2.shape[1]):
    rho_x00 = np.sqrt(rho2[:,iparticle]) - rho_hill
    sign = np.sign(rho_x00).astype(int)
    diff = np.diff(sign)
    idxdiff = np.where(diff > 0)[0]
    if len(idxdiff) > 0:
        thill.append(tout[idxdiff[0]])
thill = np.array(thill)
print "   Median Time To Hill Radius = %.2e yr" % np.median(thill)

print "// Computing Time To Hill Radius (Semi-Major Axis Bin)"
thill_abins = []
for ii in [ 0, 1, 2, 3 ]:
    thill_loc = []
    for iparticle in range(rho_abins[ii].shape[0]):
        rho_x00 = rho_abins[ii][iparticle,:] - rho_hill
        sign = np.sign(rho_x00).astype(int)
        diff = np.diff(sign)
        idxdiff = np.where(diff)[0]
        if len(idxdiff) > 0:
            thill_loc.append(tout[idxdiff[0]])
    thill_loc = np.asarray(thill_loc)
    thill_abins.append(thill_loc)
    print "   Median Time To Hill Radius = %.2e yr" % \
        np.median(thill_abins[-1])

# Write
print "// Writing %s" % args.outfile
np.savez(args.outfile, \
    te01 = te01, \
    te02 = te02, \
    slope03 = slope03, \
    te01_abins = te01_abins, \
    te02_abins = te02_abins, \
    slope03_abins = slope03_abins, \
    rho = rho, \
    thill = thill, \
    thill_abins = thill_abins )

# Done
print "// Done"
