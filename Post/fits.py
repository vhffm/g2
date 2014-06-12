"""
Read Metric Distance. Compute Different Fits.

1) First Exponential
2) Second Exponential
3) Final Power law
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit

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
ttx = np.repeat([tout], rho2.shape[1], axis=0).T
rho = np.sqrt(rho2/rho2[0,:])

# Fit 01
print "// Fitting First Exponential"
ttx_loc = ttx[rho<1.0e7]
rho_loc = rho[rho<1.0e7]
fitParams01, fitCovariances01 = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                                          p0 = [ 40000.0, 0.01, -10000.0 ])
te01 = 1.0/fitParams01[1]
print "   E-Folding Time = %.2e" % te01

# Fit 02
print "// Fitting Second Exponential"
bool_01 = rho>1.0e7
bool_02 = rho<=1.0e11
bool_cb = np.logical_and(bool_01, bool_02)
ttx_loc = ttx[bool_cb]
rho_loc = rho[bool_cb]
fitParams02, fitCovariances02 = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                                          p0 = [ 4.0e9, 0.01, -1.0e10 ])
te02 = 1.0/fitParams02[1]
print "   E-Folding Time = %.2e" % te02

# Fit 03
print "// Fitting Power Law"
ttx_loc = ttx[rho>1.0e11]
rho_loc = rho[rho>1.0e11]
fitParams03, fitCovariances03 = curve_fit(fit_func_pow, ttx_loc, rho_loc, \
                                          p0 = [ 1.0e12, 0.3, -8.0e12 ])
slope03 = fitParams03[1]
print "   Slope = %.2e" % slope03

# Write
print "// Writing %s" % args.outfile
np.savez(args.outfile, \
    fitParams01 = fitParams01, fitCovariances01, \
    fitParams02 = fitParams02, fitCovariances02, \
    fitParams03 = fitParams03, fitCovariances03, \
    te01 = te01, \
    te02 = te02, \
    slope03 = slope03 )

# Done
print "// Done"
