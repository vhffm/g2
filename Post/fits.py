"""
Read Metric Distance. Compute Different Fits. Note Easter Egg.

1) First Exponential
2) Second Exponential
3) Final Power Law
4) Time To Hill Radius
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
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='Fits.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Read
print "// Reading %s" % args.infile
npz = np.load(args.infile)
tout = npz["tout"]
rho2 = npz["rho2"]
# ttx = np.repeat([tout], rho2.shape[1], axis=0).T
ttx = tout
rho = np.sqrt(rho2/rho2[0,:])
rho = mstats.gmean(rho, axis=1)

# #############################################################################
# Fit 01
# #############################################################################
print "// Fitting First Exponential"
ttx_loc = ttx[rho<1.0e7]
rho_loc = rho[rho<1.0e7]
try:
    fitParams01, _ = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                               p0 = [ 40000.0, 0.01, -10000.0 ])
    te01 = 1.0/fitParams01[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    te01 = np.nan
print "   E-Folding Time = %.2e yr" % te01

# #############################################################################
# Fit 02
# #############################################################################
print "// Fitting Second Exponential"
bool_01 = rho>1.0e7
bool_02 = rho<=1.0e11
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

# #############################################################################
# Fit 03
# #############################################################################
print "// Fitting Power Law"
ttx_loc = ttx[rho>1.0e11]
rho_loc = rho[rho>1.0e11]
try:
    fitParams03, _ = curve_fit(fit_func_pow, ttx_loc, rho_loc, \
                               p0 = [ 1.0e12, 0.3, -8.0e12 ])
    slope03 = fitParams03[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    slope03 = np.nan
print "   Slope = %.2e" % slope03

# #############################################################################
# Time To Hill Radius
# #############################################################################
print "// Computing Time To Hill Radius"
m = 5.0 * C.mearth / rho2.shape[1]
Rh = 1.0 * (m / 3.0 / C.msun)**(1.0/3.0)
rho_hill = orbital_distance(Rh)
print "   (N, m, R_hill) = (%i, %.2e kg, %.2e AU)" % (rho2.shape[1], m, Rh)

thill = []
for iparticle in range(rho2.shape[1]):
    rho_x00 = np.sqrt(rho2[:,iparticle]) - rho_hill
    sign = np.sign(rho_x00).astype(int)
    diff = np.diff(sign)
    idxdiff = np.where(diff)[0]
    if len(idxdiff) > 0:
        thill.extend(tout[idxdiff[0]])
thill = np.array(thill)

print "   Median Time To Hill Radius = %.2e yr" % np.median(thill)

# Write
print "// Writing %s" % args.outfile
np.savez(args.outfile, \
    te01 = te01, \
    te02 = te02, \
    slope03 = slope03, \
    rho = rho, \
    thill = thill )

# Done
print "// Done"
