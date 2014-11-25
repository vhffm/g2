"""
Plot Apo- & Perihelia. Stacked.
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

# Load Colormaps
c3 = b2m.get_map('Dark2', 'Qualitative', 3)
cx = b2m.get_map('BuPu', 'Sequential', 9, reverse=True)

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

# Surpress Invalid Filtering Errors (e.g., np.nan>0)
np.seterr(invalid='ignore')

# Test Set (Deprecated)
# Runs (from:to:steps)
# nsteps = np.mgrid[0:1000000+100000:100000]
# nsteps = np.mgrid[1000000-100000:1000000+100000:100000]

# Load Collision & Ejection Files, Abort If Absent
ctime_all = []
cid01_all = []
cid02_all = []
etime_all = []
etype_all = []
print ""
for idir, cdir in enumerate(dirs):
    try:
        print "// Loading Collisions (%s)" % cdir
        c = np.loadtxt("%s/Collisions%s.dat" % (cdir, run_names[idir]), \
                       ndmin=2)
        print "// Loading Ejections (%s)" % cdir
        e = np.loadtxt("%s/Ejections%s.dat" % (cdir, run_names[idir]), \
                       ndmin=2)
    except:
        print "// Missing Collion or Ejection File (%s). Abort." % cdir
        sys.exit()

    # Prep Collisions, Ejections To Filter
    # Set empty 2d array if there are no collisions/ejections
    if c.shape[0]==0:
        ctime = np.array([[],[]], dtype=int)
        cid01 = np.array([[],[]], dtype=int)
        cid02 = np.array([[],[]], dtype=int)
    else:
        ctime = c[:,0]
        cid01 = c[:,1].astype(int)
        cid02 = c[:,13].astype(int)
    if e.shape[0]==0:
        etime = np.array([[],[]], dtype=int)
        etype = np.array([[],[]], dtype=int)
    else:
        etime = e[:,0]
        etype = e[:,13].astype(int)

    # Append
    ctime_all.append(ctime)
    cid01_all.append(cid01)
    cid02_all.append(cid02)
    etime_all.append(etime)
    etype_all.append(etype)

#
# Loop Me, Baby
# https://www.youtube.com/watch?v=DVpq59F9I9o
# 
print ""
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Loop Directories
    for idir, cdir in enumerate(dirs):
        print "//                %s" % cdir

        # Reset Counters
        if idir == 0:
            Nejected, Ninfall, Ncoll = 0, 0, 0
            coll_01, coll_02, coll_03, coll_04 = 0, 0, 0, 0

        with h5py.File("%s/Snapshot_%012d.hdf5" % (cdir, nstep), "r") as f5:
            aloc = f5["particles"]["a"][()]
            eloc = f5["particles"]["e"][()]
            mloc = f5["particles"]["m"][()]
            pidx_loc = f5["particles"]["pid"][()]
            tout = f5.attrs["tout"]

        # Global IDs
        if idir == 0:
            pidx = pidx_loc
        else:
            pidx = np.concatenate([pidx, pidx_loc])

        # Massive/Massless
        # Massive From First Directory Only
        if idir == 0:
            bool_test = mloc==0.0
            bool_mass = mloc>0.0
        else:
            bool_test = np.concatenate([bool_test, mloc==0])

        # Compute Apo-/Perihelia
        eloc[np.logical_or(eloc>1, eloc<0)] = np.nan
        if idir == 0:
            rp = (1.0 - eloc) * aloc
            ra = (1.0 + eloc) * aloc
        else:
            rp = np.concatenate([rp, (1.0 - eloc) * aloc])
            ra = np.concatenate([ra, (1.0 + eloc) * aloc])

        # Keep Masses
        if idir == 0:
            mall = mloc
        else:
            mall = np.concatenate([mall, mloc])

        #
        # Slice & Count Collisions, Ejections
        #

        # Apply Time Slicing
        cid01_loc = cid01_all[idir][ctime_all[idir]<=tout]
        cid02_loc = cid02_all[idir][ctime_all[idir]<=tout]
        etype_loc = etype_all[idir][etime_all[idir]<=tout]

        # Inner SS Collisions [Mercury, Venus, Earth, Mars]
        coll_01 += np.sum(np.logical_or(cid01_loc==0, cid02_loc==0))
        coll_02 += np.sum(np.logical_or(cid01_loc==1, cid02_loc==1))
        coll_03 += np.sum(np.logical_or(cid01_loc==2, cid02_loc==2))
        coll_04 += np.sum(np.logical_or(cid01_loc==3, cid02_loc==3))

        # Total Collisions, Ejected, Infall
        Nejected += np.sum(etype_loc==-3)
        Ninfall += np.sum(etype_loc==-2)
        Ncoll += np.sum(ctime<=tout)

    # Centr Covers Jupiter/Saturn
    bool_inner = np.logical_and(bool_test, pidx<25000)
    bool_centr = np.logical_and(bool_test, \
                                np.logical_and(pidx>=25000, pidx<80000))
    bool_outer = np.logical_and(bool_test, pidx>=80000)

    # 
    # Plot
    #

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Inner Solar System
    ax.fill_between([1.0e-1, 1.0e6], [1.0e-1, 1.0e-1], [1.7, 1.7], \
                    facecolor=c3.mpl_colors[0], alpha=0.05, lw=0.5)
    ax.fill_between([1.0e-1, 1.0e6], [1.0e-1, 1.0e-1], [1.1, 1.1], \
                    facecolor=c3.mpl_colors[0], alpha=0.05, lw=0.5)
    ax.fill_between([1.0e-1, 1.7], [1.0e-1, 1.0e-1], [1.0e2, 1.0e2], \
                    facecolor=c3.mpl_colors[0], alpha=0.05, lw=0.5)
    ax.fill_between([1.0e-1, 1.1], [1.0e-1, 1.0e-1], [1.0e2, 1.0e2], \
                    facecolor=c3.mpl_colors[0], alpha=0.05, lw=0.5)

    # Massive
    ax.scatter(ra[bool_mass], rp[bool_mass], \
               s=(mall[bool_mass]/(C.mmercury/C.msun)+30)**(2./3.), \
               c=c3.mpl_colors[1], alpha=0.8, lw=0.5)

    # Data Points
    ax.scatter(ra[bool_test], rp[bool_test], \
               s=1.0, \
               c=c3.mpl_colors[0], alpha=0.8, edgecolor=c3.mpl_colors[0])

    # Outer Test Particles
    ax.scatter(ra[bool_outer], rp[bool_outer], \
               s=1.0, \
               c=c3.mpl_colors[0], alpha=0.8, edgecolor=cx.mpl_colors[4])

    # Middle Test Particles
    ax.scatter(ra[bool_centr], rp[bool_centr], \
               s=1.0, \
               c=c3.mpl_colors[0], alpha=0.8, edgecolor=cx.mpl_colors[2])

    # Inner Test Particles
    ax.scatter(ra[bool_inner], rp[bool_inner], \
               s=1.0, \
               c=c3.mpl_colors[0], alpha=0.8, edgecolor=cx.mpl_colors[0])

    # Reference Circular Orbit Line
    ax.plot([0.1,1,10,100], [0.1,1,10,100], c='k', alpha=0.2, lw=0.5)

    # Timer
    ax.text(0.03, 0.92, r"%.2e yr" % tout, \
            horizontalalignment='left', color='black', size='small', \
            transform=ax.transAxes)

    # Number of Particles
    ax.text(0.75, 0.10, r"$N(r_p<1.7) = %i$" % \
            rp[np.logical_and(rp<1.7, bool_test)].shape[0], \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)
    ax.text(0.75, 0.04, r"$N(r_p<1.1) = %i$" % \
            rp[np.logical_and(rp<1.1, bool_test)].shape[0], \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)

    # Impact/Ejection/Infall Counter
    ax.text(0.75, 0.94, r"$N_{\mathrm{Ejected}} = %i$" % \
            Nejected, \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)
    ax.text(0.75, 0.88, r"$N_{\mathrm{Infall}} = %i$" % \
            Ninfall, \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)
    ax.text(0.75, 0.82, r"$N_{\mathrm{Collision}} = %i$" % \
            Ncoll, \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)

    # MVEM Collisions
    ax.text(0.3, 0.92, r"$N_\mathrm{Impacts}(M,V,E,M) = (%i,%i,%i,%i)$" % \
            (coll_01, coll_02, coll_03, coll_04), \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)

    # Remaining Particle Counter
    ax.text(0.03, 0.8, r"$N_{\mathrm{Mass}} = %i$" % \
            np.sum(bool_mass), \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)
    ax.text(0.03, 0.74, r"$N_{\mathrm{Test}} = %i$" % \
            np.sum(bool_test), \
            horizontalalignment='left', color='black', \
            transform=ax.transAxes)

    # Labels, Scales, Limits
    ax.set_xlim([1.0e-1,1.0e6])
    ax.set_ylim([1.0e-1,1.0e2])
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlabel('r_a = a (1+e)')
    # ax.set_ylabel('r_p = a (1-e)')
    ax.set_xlabel('Apohelion (Farthest) (AU)')
    ax.set_ylabel('Perihelion (Closest) (AU)')
    ax.set_title("%s" % args.tag, fontsize="small")
    
    fig.savefig("rpra_%012d.png" % nstep)
    plt.close(fig)
    print ""

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
