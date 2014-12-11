"""
Plot XY. 6 Panels. Test and Massive Particles.
Plot Orbits. Massive Particles.
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import constants as C
import brewer2mpl as b2m
import argparse
import sys
import os
import h5py
from glob import glob
from time import gmtime, strftime
import kepler_helpers as kh

# Load Colormaps
c3 = b2m.get_map('Dark2', 'Qualitative', 3)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--tag", \
                    help='Title Tag. Defaults Current Working Dir.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--all', action='store_true', \
                   help="Use All Snapshots.")
group1.add_argument('--custom', type=int, nargs='+', \
                   help="Use Custom Snapshot Range.")
args = parser.parse_args()

# List of Directories.
# Format:
# run_name_01,/path/to/director/01
# run_name_02,/path/to/director/02
# ...
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    run_names = []
    for line in lines:
        dirs.append(line.split(",")[1])
        run_names.append(line.split(",")[0])
    print "// Reading %i Directories" % len(dirs)

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Output set must be defined by three numbers."
        sys.exit()

# Full Set. Based On First Directory.
if args.all:
    nsteps = []
    globs = glob("%s/Snapshot_*.hdf5" % dirs[0])
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)
    print "// Reading %i Outputs Per Directory" % len(nsteps)

# Custom Set
if args.custom:
    # Build Output Number Array (From Input)
    nsteps = \
        np.mgrid[args.custom[0]:args.custom[1]+args.custom[2]:args.custom[2]]
    print "// Using Outputs %012d:%012d:%012d" % \
        ( args.custom[0], args.custom[1], args.custom[2] )

# Set Unset Tag
if not args.tag:
    args.tag = os.getcwd()

# Test Set (Deprecated)
# Runs (from:to:steps)
# nsteps = np.mgrid[0:1000000+100000:100000]
# nsteps = np.mgrid[1000000-100000:1000000+100000:100000]

#
# Loop Me, Baby
# https://www.youtube.com/watch?v=DVpq59F9I9o
# 
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(2,3)
    fig.set_size_inches(12,6)
    axflat = axarr.flatten()

    # Loop Runs
    for idir, cdir in enumerate(dirs):
        try:
            with h5py.File("%s/Snapshot_%012d.hdf5" % (cdir, nstep), "r") as f5:
                # Record Time
                tout = f5.attrs["tout"]
                
                # Plot Test Particles
                axflat[idir].scatter(f5["particles/x"][9:], \
                                     f5["particles/y"][9:], \
                                     s=1, alpha=0.5, \
                                     edgecolor=c3.mpl_colors[2], \
                                     color=c3.mpl_colors[2])
                
                # Compute Ellipses
                xell, yell, zell = \
                    kh.compute_ellipseX(f5["particles/a"][:9], \
                                        f5["particles/e"][:9], \
                                        f5["particles/i"][:9], \
                                        f5["particles/Omega"][:9], \
                                        f5["particles/omega"][:9])
                    
                # Plot Planets
                axflat[idir].plot(xell.T, yell.T, color="k", linewidth=1.0)
                axflat[idir].scatter(f5["particles/x"][:9], \
                                     f5["particles/y"][:9], \
                                     c="k", marker="+", s=12**2, \
                                     alpha=1.0, zorder=3)
                
                # Title
                axflat[idir].set_title("%s" % run_names[idir])
        except:
            pass

        # Set Scale
        # for ax in axflat:
        #     ax.set_xscale("symlog", xlinthresh=0.1)
        #     ax.set_yscale("symlog", ylinthresh=0.1)

        # Set Limits
        for ax in axflat:
            ax.set_aspect("equal")
            ax.set_xlim([-60,60])
            ax.set_ylim([-60,60])
            
        # Remove Useless Labels
        for ii in [ 1, 2, 4, 5 ]:
            plt.setp(axflat[ii].get_yticklabels(), visible=False)

        for ii in [ 0, 1, 2 ]:
            plt.setp(axflat[ii].get_xticklabels(), visible=False)
            
        # Some Labels
        axarr[1,0].set_xlabel("X (AU)")
        axarr[1,0].set_ylabel("Y (AU)")
            
        # Fiture Title
        fig.suptitle("%s *** %.2e yr *** %012d steps" % \
            (os.getcwd(), tout, nstep))
    
    fig.savefig("xy_panels_%012d.png" % nstep)
    plt.close(fig)

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
