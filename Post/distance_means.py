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
parser.add_argument('--semi_major_axis', action='store_true', \
                    help="Filter/Bin/Save in Semi-Major Axis")
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
        a1 = npz["a1"]
        rho2 = npz["rho2"]
        dv2 = npz["dv2"]
        ds2 = npz["ds2"]
        tout = npz["tout"]
        first = False
    else:
        a1 = np.concatenate((a1, npz["a1"]), axis=1)
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

# Filter/Bin Semi-Major Axis
rho2_means_a = []
if args.semi_major_axis:
    print "// Filtering Distances. Be Patient."

    tags = [ "a<1", "1<a<2", "2<a<3", "a>3" ]
    for ii in [ 0, 1, 2, 3 ]:
        print "   %s" % tags[ii]
        
        # Copy, Compute Ratio
        a1_loc = a1.copy()
        ratio_rho2_loc = ratio_rho2.copy()
        
        # Set Filters
        if ii == 0:
            mybool = a1_loc<1.0
        elif ii == 1:
            mybool = np.logical_and(a1_loc>=1.0, a1_loc<2.0) 
        elif ii == 2:
            mybool = np.logical_and(a1_loc>=2.0, a1_loc<3.0) 
        elif ii == 3:
            mybool = a1_loc>=3.0

        # Apply Filters
        ratio_rho2_loc[~mybool] = np.nan
        ratio_rho2_loc = np.ma.masked_invalid(ratio_rho2_loc)
        
        # Mean
        rho2_mean_loc = mstats.gmean(ratio_rho2_loc, axis=1)
        
        # Append
        rho2_means_a.append(rho2_mean_loc)

# Save
print "// Saving Files"
np.savez(args.outfile, tout = tout, \
         precision = npz["precision"][()], \
         version = int(npz["version"][()])+1, \
         rho2_mean = rho2_mean, \
         rho2_means_a = rho2_means_a, \
         ds2_mean = ds2_mean, dv2_mean = dv2_mean )

# Done
print "// Done"
