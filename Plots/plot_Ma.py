"""
Plot Semi-Major Axis vs. Mass.
Use List to Specifiy Directories.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['lines.linewidth'] = 1.0
import matplotlib.pyplot as plt
from g2_helpers import twopi, r2d
import os
import sys
from copy import copy
from matplotlib.ticker import MaxNLocator
from time import gmtime, strftime

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--coffset', type=int, default=0, \
                    help="Marker Colour Offset")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshots.")
args = parser.parse_args()

# Colors
dgrey = (0.5, 0.5, 0.5, 1.0)
grey  = (0.7, 0.7, 0.7, 1.0)
mred  = (1.0, 0.7, 0.7, 1.0)
mblue = (0.7, 0.7, 1.0, 1.0)
lgrey = (0.9, 0.9, 0.9, 1.0)
black = (0.0, 0.0, 0.0, 1.0)

# Lines
mpl.rcParams['lines.linewidth'] = 0.5

# Legend
# mpl.rcParams['legend.handlelength']  = 2.9
# mpl.rcParams['legend.handlelength']  = 0.5
mpl.rcParams['legend.frameon']       = True
mpl.rcParams['legend.numpoints']     = 1
mpl.rcParams['legend.scatterpoints'] = 1

# Lighten labels
# http://matplotlib.org/users/customizing.html
mpl.rcParams['axes.labelcolor']  = dgrey
mpl.rcParams['xtick.color']      = dgrey
mpl.rcParams['ytick.color']      = dgrey
mpl.rcParams['axes.edgecolor']   = dgrey

# Adjust ticks
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.minor.size'] = 0
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.minor.size'] = 0

# Adjust Font Size
mpl.rcParams['xtick.labelsize']  = 'small'
mpl.rcParams['ytick.labelsize']  = 'small'
mpl.rcParams['axes.labelsize']   = 'medium'
mpl.rcParams['axes.titlesize']   = 'medium'

# Masses in kg
m_earth   = 5.972e24
m_sun     = 1.99e30
m_jupiter = 1.89e27

# Cutoff mass (in solar masses)
m_cutoff = 2.0e23 # kg
m_cutoff = m_cutoff / m_sun

print "!!"
print "!! NOTICE - PARTICLE IDS 2000 AND 2001 ARE IGNORED !!"
print "!!          THEY ARE USUALLY JUPITER AND SATURN    !!"
print "!!"

# Style Dictionary
# http://matplotlib.org/api/markers_api.html#module-matplotlib.markers
snap_c = [ 'g', 'r', 'b', 'm', 'k', \
           'g', 'r', 'b', 'm', 'k', \
           'g', 'r', 'b', 'm', 'k', \
           'g', 'r', 'b', 'm', 'k', \
           'g', 'r', 'b', 'm', 'k' ]
snap_s = [ 'o', 'o', 'o', 'o', 'o', \
           'd', 'd', 'd', 'd', 'd', \
           'v', 'v', 'v', 'v', 'v', \
           '+', '+', '+', '+', '+', \
           '^', '^', '^', '^', '^', \
           's', 's', 's', 's', 's' ]

# List of Directories
# @todo Check for existence
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    tags = []
    for line in lines:
        dirs.append(line.split(",")[-1])
        tags.append(line.split(",")[0])
    print "// Reading %i Directories" % len(dirs)

#
# BIG WARNING
#
# We assume that all directories have the same snapshot numbers sitting around.
# The timestep of the simulations is also equal.
#

# Build Snapshot Number Array (From First Dir)
print "// Building Snapshot Array"
if args.all:
    globs = glob(dirs[0] + "/" + "Snapshot_*.npz")
    globs = sorted(globs)
    nsteps = np.zeros(len(globs))
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.npz')[0].split('/')[-1].split('_')[1])
if args.test:
    nsteps = np.mgrid[3600000000:3630000000:1000000]
if args.custom:
    # Build Snapshot Number Array (From Input)
    nsteps = range(args.custom[0], \
                   args.custom[1] + args.custom[2], \
                   args.custom[2])
print "// Found %i Snapshots" % len(nsteps)

# Loop Snapshots, Draw Figures
print "// Start Processing Snapshots"
for nstep in nsteps:
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1)
    for idir, dirchar in enumerate(dirs):
        nblocks = 2
        dirs_per_block = len(dirs)/nblocks
        nsweep = idir/dirs_per_block
        try:
            # Load Snapshot Data
            npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
            snapshot = npz['snapshot'][()]
            pa = np.zeros(snapshot.nparticles)
            pm = np.zeros(snapshot.nparticles)
            for ipart, particle in enumerate(snapshot.particles):
                if particle.m > m_cutoff:
                    pa[ipart] = particle.a
                    pm[ipart] = particle.m * m_sun / m_earth
                else:
                    pa[ipart] = np.nan
                    pm[ipart] = np.nan
                if particle.id >= 2000:
                    pa[ipart] = np.nan
                    pm[ipart] = np.nan
            # Plot Snapshot
            # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.scatter
            # Actual Data
            ax1.scatter(pa, pm, s=75, \
                        c=snap_c[nsweep+args.coffset], marker=snap_s[nsweep], \
                        edgecolors='none', \
                        alpha=0.5)
            # Fake Datapoint For Legend
            if nsweep != (idir+1)/dirs_per_block:
                ax1.scatter([-10], [-10], s=75, \
                            c=snap_c[nsweep+args.coffset], \
                            marker=snap_s[nsweep], \
                            edgecolors='none', \
                            alpha=0.5, label=tags[idir])
        except IOError:
            print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                  (dirs[idir], nstep)
    # Style Figures
    ax1.grid(False)
    ax1.set_xlim([0,5])
    ax1.set_ylim([0,1.8])
    ax1.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
    ax1.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
    ax1.set_xlabel('Semi-Major Axis (AU)')
    ax1.set_ylabel('Mass (Earth Masses)', labelpad=20)
    ax1.set_title('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    # Disable Legend For Now
    ax1.legend(prop={'size':'medium'}, loc='upper right')
    # Save Figures
    fig1.savefig('Ma_%012d.png' % snapshot.nstep)
    plt.close(fig1)
print "// Done"
