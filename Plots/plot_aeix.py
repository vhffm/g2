"""
Plot Semi-Major Axis vs. Eccentricity, Inclination.
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

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# Style Dictionary
# http://matplotlib.org/api/markers_api.html#module-matplotlib.markers
snap_c = [ 'r', 'g', 'b', 'm', 'k', \
           'r', 'g', 'b', 'm', 'k', \
           'r', 'g', 'b', 'm', 'k', \
           'r', 'g', 'b', 'm', 'k', \
           'r', 'g', 'b', 'm', 'k' ]
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
    for line in lines:
        dirs.append(line)
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
    nsteps = np.array([args.custom])
print "// Found %i Snapshots" % len(nsteps)

# Scan Limits
print "// Scanning Limits"
first = True
for nstep in nsteps:
    for idir, dirchar in enumerate(dirs):
        try:
            npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
            snapshot = npz['snapshot'][()]
            for particle in snapshot.particles:
                if first:
                    amax = particle.a; amin = particle.a
                    eccmax = particle.ecc; eccmin = particle.ecc
                    incmax = particle.inc; incmin = particle.inc
                    first = False
                else:
                    if particle.a > amax: amax = particle.a
                    if particle.a < amin: amin = particle.a
                    if particle.ecc > eccmax: eccmax = particle.ecc
                    if particle.ecc < eccmin: eccmin = particle.ecc
                    if particle.inc > incmax: incmax = particle.inc
                    if particle.inc < incmin: incmin = particle.inc
        except IOError:
            print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                  (dirs[idir], nstep)

# Loop Snapshots, Draw Figures
for nstep in nsteps:
    print "// Processing Snapshot %012d/%012d" % (nstep, nsteps[-1])
    fig1 = plt.figure(); fig2 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1); ax2 = fig2.add_subplot(1,1,1)
    for idir, dirchar in enumerate(dirs):
        try:
            # Load Snapshot Data
            npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
            snapshot = npz['snapshot'][()]
            pa = np.zeros(snapshot.nparticles)
            pecc = np.zeros(snapshot.nparticles)
            pinc = np.zeros(snapshot.nparticles)
            for ipart, particle in enumerate(snapshot.particles):
                pa[ipart] = particle.a
                pecc[ipart] = particle.ecc
                pinc[ipart] = particle.inc * r2d
            # Plot Snapshot
            ax1.plot(pa, pecc, snap_c[idir] + snap_s[idir], \
                     label=dirs[idir], alpha=0.5, markersize=3)
            ax2.plot(pa, pinc, snap_c[idir] + snap_s[idir], \
                     label=dirs[idir], alpha=0.5, markersize=3)
        except IOError:
            print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                  (dirs[idir], nstep)
    # Style Figures
    ax1.grid(True); ax2.grid(True)
    ax1.set_xlim([0,amax]); ax2.set_xlim([0,amax]); 
    ax1.set_ylim([0,eccmax]); ax2.set_ylim([0,incmax*r2d])
    ax1.set_xlabel('a [AU]', size='small')
    ax2.set_xlabel('a [AU]', size='small')
    ax1.set_ylabel('e [-]', rotation='horizontal', size='small')
    ax2.set_ylabel('inc [deg]', rotation='horizontal', size='small')
    ax1.set_title('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    ax2.set_title('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    ax1.legend(prop={'size':'xx-small'}, loc='best')
    ax2.legend(prop={'size':'xx-small'}, loc='best')
    # Save Figures
    fig1.savefig('aex_%012d.png' % snapshot.nstep)
    fig2.savefig('aix_%012d.png' % snapshot.nstep)
    plt.close(fig1); plt.close(fig2)
