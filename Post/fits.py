"""
Read Metric Distance. Compute Different Fits.

1) First Exponential
2) Second Exponential
3) Final Power law
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit
import scipy.stats as sps

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

# Init Arrays
te01 = []
te02 = []
slope03 = []

# Debug Info
print "// Fitting %i Particles" % rho.shape[1]

# Fit 01
print "// Fitting First Exponential"
for ipart in range(rho.shape[1]):
    ttx_loc = ttx[:,ipart]
    rho_loc = rho[:,ipart]
    ttx_loc = ttx_loc[rho_loc<1.0e7]
    rho_loc = rho_loc[rho_loc<1.0e7]
    try:
        fitParams01, _ = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                                   p0 = [ 40000.0, 0.01, -10000.0 ])
        te01.append(1.0/fitParams01[1])
    except (RuntimeError, TypeError):
        print "!! Fit Failed"
        te01.append(np.nan)

te01 = np.array(te01)
print "   Median E-Folding Time = %.2e" % sps.nanmedian(te01)

# Fit 02
print "// Fitting Second Exponential"
for ipart in range(rho.shape[1]):
    ttx_loc = ttx[:,ipart]
    rho_loc = rho[:,ipart]
    bool_01 = rho_loc>1.0e7
    bool_02 = rho_loc<=1.0e11
    bool_cb = np.logical_and(bool_01, bool_02)
    ttx_loc = ttx_loc[bool_cb]
    rho_loc = rho_loc[bool_cb]
    try:
        fitParams02, _ = curve_fit(fit_func_exp, ttx_loc, rho_loc, \
                                   p0 = [ 4.0e9, 0.01, -1.0e10 ])
        te02.append(1.0/fitParams02[1])
    except (RuntimeError, TypeError):
        print "!! Fit Failed"
        te02.append(np.nan)

te02 = np.array(te02)
print "   Median E-Folding Time = %.2e" % sps.nanmedian(te02)

# Fit 03
print "// Fitting Power Law"
for ipart in range(rho.shape[1]):
    ttx_loc = ttx[:,ipart]
    rho_loc = rho[:,ipart]
    ttx_loc = ttx_loc[rho_loc>1.0e11]
    rho_loc = rho_loc[rho_loc>1.0e11]
    try:
        fitParams03, _ = curve_fit(fit_func_pow, ttx_loc, rho_loc, \
                                   p0 = [ 1.0e12, 0.3, -8.0e12 ])
        slope03.append(fitParams03[1])
    except (RuntimeError, TypeError):
        print "!! Fit Failed"
        slope03.append(np.nan)

slope03 = np.array(slope03)
print "   Median Slope = %.2e" % sps.nanmedian(slope03)

# Write
print "// Writing %s" % args.outfile
np.savez(args.outfile, \
    te01 = te01, \
    te02 = te02, \
    slope03 = slope03 )

# Done
print "// Done"
