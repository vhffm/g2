"""
Compute Geometric Mean Over All Distances.
"""

import numpy as np
import argparse
import sys
import scipy.stats.mstats as mstats

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--outfile', default='Distances_Means.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# List of Directories
if sys.stdin.isatty():
    print "!! No File List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    files = []
    for line in lines:
        files.append(line)
    print "// Reading %i Files" % len(files)

# Read All Files
print "// Loading Files"
first = True
for file in files:
    npz = np.load(file)
    if first:
        rho2 = npz["rho2"]
        dv2 = npz["dv2"]
        ds2 = npz["ds2"]
        tout = npz["tout"]
        first = False
    else:
        rho2 = np.concatenate((rho2, npz["rho2"]), axis=1)
        dv2 = np.concatenate((dv2, npz["dv2"]), axis=1)
        ds2 = np.concatenate((ds2, npz["ds2"]), axis=1)

# Compute Ratios
print "// Computing Ratios"
ratio_rho2 = rho2/rho2[0,:]
ratio_ds2 = ds2/ds2[0,:]
ratio_dv2 = dv2[1:]/dv2[1,:]

# Compute Means
print "// Computing Means"
rho2_mean = mstats.gmean(ratio_rho2, axis=1)
dv2_mean = mstats.gmean(ratio_dv2, axis=1)
ds2_mean = mstats.gmean(ratio_ds2, axis=1)

# Save
print "// Saving Files"
np.savez(args.outfile, tout = tout, \
         rho2_mean = rho2_mean, \
         ds2_mean = ds2_mean, dv2_mena = dv2_mean )

# Done
print "// Done"
