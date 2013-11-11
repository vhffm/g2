"""
Plot Semi-Major Axis vs. Eccentricity, Inclination.
Use List to Specify Directories.
Can Handle Arbitrary Number of Directories. 
Displayed on 4x2=8 Grid. Then Stacked.
"""

from glob import glob
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['lines.linewidth'] = 1.0
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys
import argparse
from g2_helpers import r2d, mkline

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, \
                   help="Plot Custom Snapshot.")
parser.add_argument('--scale', action='store_true', \
                    help="Scale Marker Size with Particle Mass")
args = parser.parse_args()

# Throw NoScale Warning
if not args.scale:
    print "!!"
    print "!! NOTICE - NOT SCALING MARKERS WITH MASS !!"
    print "!!"

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
                    mmax = particle.m; mmin = particle.m
                    first = False
                else:
                    if particle.a > amax: amax = particle.a
                    if particle.a < amin: amin = particle.a
                    if particle.ecc > eccmax: eccmax = particle.ecc
                    if particle.ecc < eccmin: eccmin = particle.ecc
                    if particle.inc > incmax: incmax = particle.inc
                    if particle.inc < incmin: incmin = particle.inc
                    if particle.m > mmax: mmax = particle.m
                    if particle.m < mmin: mmin = particle.m
        except IOError:
            print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                  (dirs[idir], nstep)

# Compute Line for Marker Size(Mass)
m, n = mkline(mmin, 1.0, mmax, 36.0)

for nstep in nsteps:
    print "// Processing Snapshot %012d/%012d" % (nstep, nsteps[-1])
    fig1 = plt.figure(figsize=(16.0, 12.0)); fig2 = plt.figure(figsize=(16.0, 12.0))
    # Sweep Subplots in Steps of 8
    nsweeps = (len(dirs) - 1) / 8 + 1
    for nsweep in range(0,nsweeps):
        for irow in range(0,4):
            for icol in range(1,3,1):
                ii = irow*2 + icol
                if icol == 1: jj = irow*2 + icol + 1 
                if icol == 2: jj = irow*2 + icol - 1
                # print "%i,%i,%i,%i" % (irow, icol, ii, jj)

                # Always Generate Subplots
                # Can Be Empty
                ax1 = fig1.add_subplot(4,2,jj)
                ax2 = fig2.add_subplot(4,2,jj)

                # Directory Index
                idir = ii + nsweep * 8 - 1

                # Load, Plot Snapshot Data
                if idir < len(dirs):
                    try:
                        npz = np.load('%s/Snapshot_%012d.npz' % (dirs[idir], nstep))
                        snapshot = npz['snapshot'][()]
                        pa = np.zeros(snapshot.nparticles)
                        pecc = np.zeros(snapshot.nparticles)
                        pinc = np.zeros(snapshot.nparticles)
                        pm = np.zeros(snapshot.nparticles)
                        for ipart, particle in enumerate(snapshot.particles):
                            pa[ipart] = particle.a
                            pecc[ipart] = particle.ecc
                            pinc[ipart] = particle.inc * r2d
                            pm[ipart] = particle.m
                        # Plot Snapshot
                        if args.scale: s = pm * m + n
                        else: s = 3
                        ax1 = fig1.add_subplot(4,2,jj)
                        ax2 = fig2.add_subplot(4,2,jj)
                        ax1.scatter(pa, pecc, s=s**2., \
                                    c=snap_c[nsweep+1], marker=snap_s[nsweep+1], \
                                    edgecolors='none', \
                                    alpha=0.5, label=dirs[idir])
                        ax2.scatter(pa, pinc, s=s**2., \
                                    c=snap_c[nsweep+1], marker=snap_s[nsweep+1], \
                                    edgecolors='none', \
                                    alpha=0.5, label=dirs[idir])
                        ax1.legend(prop={'size':'x-small'}, loc='best', scatterpoints=1)
                        ax2.legend(prop={'size':'x-small'}, loc='best', scatterpoints=1)
                    except IOError:
                        print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                              (dirs[idir], nstep)

                # Style Subplots
                # Can Be Empty
                ax1.grid(True); ax2.grid(True)
                ax1.set_xlim([0,5]); ax2.set_xlim([0,5])
                ax1.set_ylim([0,0.2]); ax2.set_ylim([0,incmax*r2d])
                ax1.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax2.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax1.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax2.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                if irow == 3:
                    ax1.set_xlabel('a [au]')
                    ax2.set_xlabel('a [au]')
                if icol == 2:
                    ax1.set_ylabel('ecc [-]', labelpad=20)
                    ax2.set_ylabel('inc [deg]', labelpad=20)

    fig1.subplots_adjust(hspace=0,wspace=0)
    fig2.subplots_adjust(hspace=0,wspace=0)
    fig1.suptitle('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    fig2.suptitle('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    fig1.savefig('aex_%012d.png' % snapshot.nstep)
    fig2.savefig('aix_%012d.png' % snapshot.nstep)
    plt.close(fig1); plt.close(fig2)
