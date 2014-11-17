"""
Plot XY/Orbits for Inner, Outer Planets + Kuiper Belt Disk.
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import constants as C
import brewer2mpl as b2m
import argparse
import sys
import os
from glob import glob
from time import gmtime, strftime

# Load Colormaps
c3 = b2m.get_map('Dark2', 'Qualitative', 3)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--tag", \
                    help='Title Tag.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--all', action='store_true', \
                   help="Reduce Full Set of Snapshots.")
group1.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Output set must be defined by three numbers."
        sys.exit()

# Full Set
if args.all:
    nsteps = []
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)
    print "// Reading %i Outputs" % len(nsteps)

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

#
# Loop Me, Baby
# https://www.youtube.com/watch?v=DVpq59F9I9o
# 
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    npz = np.load("Snapshot_%012d.npz" % nstep)
    snap = npz["snapshot"][()]
    plist = snap.particles

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,aspect='equal')

    # Loop Particles
    for p in plist:
        # Inner Planets
        if int(p.id) in [ 0, 1, 2, 3 ]:
            ax.plot(p.xell, p.yell, lw=1.0, alpha=0.6, c=c3.mpl_colors[0])
        # Outer Planets
        elif int(p.id) in [ 4, 5, 6, 7, 8 ]:
            ax.plot(p.xell, p.yell, lw=2.0, alpha=0.6, c=c3.mpl_colors[1])
        # Kuiper Belt
        else:
            ax.scatter(p.x, p.y, s=1, alpha=1.0, \
                       edgecolor=c3.mpl_colors[2], color=c3.mpl_colors[2])

    # Timer
    ax.text(0.03, 0.92, r"%.2e yr" % snap.tout, \
            horizontalalignment='left', color='black', size='small', \
            transform=ax.transAxes)
            
    # Title
    ax.set_title("%s" % args.tag, fontsize="small")
        
    # Limits, Etc
    ax.set_xlim([-50,50])
    ax.set_ylim([-50,50])
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')

    # Save
    fig.savefig("xy_kuiper_%012d.png" % nstep)
    plt.close(fig)

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
