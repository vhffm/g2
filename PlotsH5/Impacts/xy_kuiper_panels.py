"""
Plot XY/Orbits for Inner, Outer Planets + Kuiper Belt Disk.
Panels Are Different Inner Disk Edges.
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
import kepler_helpers as kh
from glob import glob
from time import gmtime, strftime

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

# Fixed Planet IDs
inner_ids = np.array([0,1,2,3])
outer_ids = np.array([4,5,6,7])

#
# Loop Me, Baby
# https://www.youtube.com/watch?v=DVpq59F9I9o
# 
last_valid_step = [ 0, 0, 0, 0, 0, 0, 0, 0 ]
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(2,4)
    fig.set_size_inches(12,6)
    axflat = axarr.flatten()

    # Loop Runs
    for idir, cdir in enumerate(dirs):
        try:
            with h5py.File("%s/Snapshot_%012d.hdf5" % (cdir, nstep), "r") as f5:
                # Record Time
                tout = f5.attrs["tout"]

                # Record Last Valid Step
                last_valid_step[idir] = nstep

                # Ellipses, Outer Planets
                outer_mask = np.in1d(f5["particles/pid"][()], outer_ids)
                outer_indices = np.where(outer_mask == True)[0]
                xell, yell, zell = \
                    kh.compute_ellipseX(f5["particles/a"][outer_mask], \
                                        f5["particles/e"][outer_mask], \
                                        f5["particles/i"][outer_mask], \
                                        f5["particles/Omega"][outer_mask], \
                                        f5["particles/omega"][outer_mask])
                axflat[idir].plot(xell.T, yell.T, \
                                  lw=1.0, alpha=0.6, \
                                  c="black")

                # Kuiper Belt
                kuiper_mask = ~np.in1d(f5["particles/pid"][()], \
                                       np.array([0,1,2,3,4,5,6,7,8]))
                x = f5["particles/x"][kuiper_mask]
                y = f5["particles/y"][kuiper_mask]
                axflat[idir].scatter(x, y, \
                                     s=1, alpha=1.0, \
                                     edgecolor=c3.mpl_colors[2], \
                                     color=c3.mpl_colors[2])

                # Title
                axflat[idir].set_title("%s" % run_names[idir])

        except:
            # Load Last Valid Snapshot
            with h5py.File("%s/Snapshot_%012d.hdf5" % \
                           (cdir, last_valid_step[idir]), "r") as f5:
                
                # Ellipses, Outer Planets
                outer_mask = np.in1d(f5["particles/pid"][()], outer_ids)
                outer_indices = np.where(outer_mask == True)[0]
                xell, yell, zell = \
                    kh.compute_ellipseX(f5["particles/a"][outer_mask], \
                                        f5["particles/e"][outer_mask], \
                                        f5["particles/i"][outer_mask], \
                                        f5["particles/Omega"][outer_mask], \
                                        f5["particles/omega"][outer_mask])
                axflat[idir].plot(xell.T, yell.T, \
                                  lw=1.0, alpha=0.6, \
                                  c="black")

                # Kuiper Belt
                kuiper_mask = ~np.in1d(f5["particles/pid"][()], \
                                       np.array([0,1,2,3,4,5,6,7,8]))
                x = f5["particles/x"][kuiper_mask]
                y = f5["particles/y"][kuiper_mask]
                axflat[idir].scatter(x, y, \
                                     s=1, alpha=1.0, \
                                     edgecolor=c3.mpl_colors[1], \
                                     color=c3.mpl_colors[1])
                
                # Title
                axflat[idir].set_title("%s / %012d" % \
                    (run_names[idir], last_valid_step[idir]))

        # Set Scale
        # for ax in axflat:
        #     ax.set_xscale("symlog", xlinthresh=0.1)
        #     ax.set_yscale("symlog", ylinthresh=0.1)

        # Set Limits
        for ax in axflat:
            ax.set_aspect("equal")
            ax.set_xlim([-50,50])
            ax.set_ylim([-50,50])

        # Remove Useless Labels
        for ii in [ 1, 2, 3, 5, 6, 7 ]:
            plt.setp(axflat[ii].get_yticklabels(), visible=False)

        for ii in [ 0, 1, 2, 3 ]:
            plt.setp(axflat[ii].get_xticklabels(), visible=False)

        # Some Labels
        axarr[1,0].set_xlabel("X (AU)")
        axarr[1,0].set_ylabel("Y (AU)")
        
        # Grids
        for ax in axflat:
            ax.grid(True, linestyle="-", alpha=0.1)
            
        # Fiture Title
        fig.suptitle("%s *** %.2e yr *** %012d steps" % \
            (os.getcwd(), tout, nstep))
    
    fig.savefig("xy_kuiper_panels_%012d.png" % nstep)
    plt.close(fig)

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
