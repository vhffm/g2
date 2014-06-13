"""
Read Metric Distance. Compute Different Fits.

1) First Exponential
2) Second Exponential
3) Final Power law
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit
import scipy.stats.mstats as mstats

# Define Fitting Functions
def fit_func_exp(t, a, b, c):
    return (a*np.exp(b*t))+c

def fit_func_pow(t, a, b, c):
    return (a*t**b)+c

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='XChaos_Distances.npz', \
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

# Fit 01
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
print "   E-Folding Time = %.2e" % te01

# Fit 02
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
print "   E-Folding Time = %.2e" % te02

# Fit 03
print "// Fitting Power Law"
ttx_loc = ttx[rho>1.0e11]
rho_loc = rho[rho>1.0e11]
try:
    fitParams03, _ = curve_fit(fit_func_pow, ttx_loc, rho_loc, \
                               p0 = [ 40000.0, 0.01, -10000.0 ])
    slope03 = 1.0/fitParams03[1]
except (RuntimeError, TypeError):
    "!! Fit Failed"
    slope03 = np.nan
print "   Slope = %.2e" % slope03

# Write
print "// Writing %s" % args.outfile
np.savez(args.outfile, \
    te01 = te01, \
    te02 = te02, \
    slope03 = slope03 )

# Done
print "// Done"
