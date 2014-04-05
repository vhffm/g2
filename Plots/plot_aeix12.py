"""
Plot Semi-Major Axis vs. Eccentricity, Inclination.
Use List to Specify Directories.
Can Handle Arbitrary Number of Directories. 
Displayed on 4x3=12 Grid. Then Stacked.
Jupiter & Saturn (2000, 2001) Ignored.
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
from time import gmtime, strftime

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--quickscan', action='store_true', \
                    help="Only Scan First And Last Snapshot [M_min, M_max].")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshots.")
parser.add_argument('--scale', action='store_true', \
                    help="Scale Marker Size with Particle Mass")
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
mpl.rcParams['legend.frameon']       = False
mpl.rcParams['legend.numpoints']     = 1
mpl.rcParams['legend.scatterpoints'] = 1

# Lighten labels
# http://matplotlib.org/users/customizing.html
mpl.rcParams['axes.labelcolor']  = grey
mpl.rcParams['xtick.color']      = grey
mpl.rcParams['ytick.color']      = grey
mpl.rcParams['axes.edgecolor']   = grey

# Adjust ticks
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.minor.size'] = 0
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.minor.size'] = 0

# Adjust Font Size
mpl.rcParams['xtick.labelsize']  = 'x-small'
mpl.rcParams['ytick.labelsize']  = 'x-small'
mpl.rcParams['axes.labelsize']   = 'small'
mpl.rcParams['axes.titlesize']   = 'medium'

# Throw NoScale Warning
if not args.scale:
    print "!!"
    print "!! NOTICE - NOT SCALING MARKERS WITH MASS !!"
    print "!!"

if not args.quickscan:
    print "!!"
    print "!! NOTICE - SCANNING ENTIRE SNAPSHOT RANGE !!"
    print "!!          POSSIBLE WASTE OF TIME         !!"
    print "!!"

print "!!"
print "!! NOTICE - PARTICLE IDS 2000 AND 2001 ARE IGNORED !!"
print "!!          THEY ARE USUALLY JUPIER AND SATURN     !!"
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
    tags = []
    for line in lines:
        dirs.append(line.split(",")[1])
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

# Scan Limits
if args.quickscan:
    print "// Quick Scanning Limits"
    print "// (%s UTC) Start Scanning Limits" % strftime("%H:%M:%S", gmtime())
    first = True
    for nstep in [ nsteps[0], nsteps[-1] ]:
        print "// (%s UTC) Scanning Snapshot %012d/%012d" % \
            (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
        for idir, dirchar in enumerate(dirs):
            try:
                npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
                snapshot = npz['snapshot'][()]
                for particle in snapshot.particles:
                    if particle.id < 2000:
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
    print "// (%s UTC) Done Scanning Limits" % strftime("%H:%M:%S", gmtime())

else:
    print "// Scanning Limits"
    print "// (%s UTC) Start Scanning Limits" % strftime("%H:%M:%S", gmtime())
    first = True
    for nstep in nsteps:
        print "// (%s UTC) Scanning Snapshot %012d/%012d" % \
            (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
        for idir, dirchar in enumerate(dirs):
            try:
                npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
                snapshot = npz['snapshot'][()]
                for particle in snapshot.particles:
                    if particle.id < 2000:
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
    print "// (%s UTC) Done Scanning Limits" % strftime("%H:%M:%S", gmtime())

# Compute Line for Marker Size(Mass)
m, n = mkline(mmin, 1.0, mmax, 36.0)

print "// (%s UTC) Start Processing Snapshots" % strftime("%H:%M:%S", gmtime())
for nstep in nsteps:
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    fig1 = plt.figure(figsize=(16.0, 12.0)); fig2 = plt.figure(figsize=(16.0, 12.0))
    # Sweep Subplots in Steps of 12
    nsweeps = (len(dirs) - 1) / 12 + 1
    for nsweep in range(0,nsweeps):
        for irow in range(0,4):
            for icol in range(1,4,1):
                ii = irow*3 + icol
                # print "%i,%i,%i,%i" % (irow, icol, ii, jj)

                # Always Generate Subplots
                # Can Be Empty
                ax1 = fig1.add_subplot(4,3,ii)
                ax2 = fig2.add_subplot(4,3,ii)

                # Directory Index
                idir = ii + nsweep * 12 - 1

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
                            if particle.id >= 2000:
                                pa[ipart] = np.nan
                                pecc[ipart] = np.nan
                                pinc[ipart] = np.nan
                                pm[ipart] = np.nan
                        # Plot Snapshot
                        if args.scale: s = pm * m + n
                        else: s = 3
                        ax1 = fig1.add_subplot(4,3,ii)
                        ax2 = fig2.add_subplot(4,3,ii)
                        ax1.scatter(pa, pecc, s=s**2., \
                                    c=snap_c[nsweep+1], marker=snap_s[nsweep+1], \
                                    edgecolors='none', \
                                    alpha=0.5, label=tags[idir])
                        ax2.scatter(pa, pinc, s=s**2., \
                                    c=snap_c[nsweep+1], marker=snap_s[nsweep+1], \
                                    edgecolors='none', \
                                    alpha=0.5, label=tags[idir])
                        if nsweep > 0:
                            ax1.legend(prop={'size':'x-small'}, loc='best')
                            ax2.legend(prop={'size':'x-small'}, loc='best')
                    except IOError:
                        print "!! Could Not Open %s/Snapshot_%012d.npz" % \
                              (dirs[idir], nstep)

                # Style Subplots
                # Can Be Empty
                ax1.grid(False); ax2.grid(False)
                ax1.set_xlim([0,5]); ax2.set_xlim([0,5])
                ax1.set_ylim([0,0.2]); ax2.set_ylim([0,incmax*r2d])
                ax1.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax2.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax1.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                ax2.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='5'))
                if irow == 3:
                    if icol == 1:
                        ax1.set_xlabel('Semi-Major Axis (AU)')
                        ax1.set_xlabel('Semi-Major Axis (AU)')
                else:
                    plt.setp(ax1.get_xticklabels(), visible=False)
                    plt.setp(ax2.get_xticklabels(), visible=False)
                if icol == 1:
                    if irow == 3:
                        ax1.set_ylabel('Eccentricity', labelpad=20)
                        ax2.set_ylabel('Inclination (Degree)', labelpad=20)
                else:
                    plt.setp(ax1.get_yticklabels(), visible=False)
                    plt.setp(ax2.get_yticklabels(), visible=False)

    fig1.subplots_adjust(hspace=0,wspace=0)
    fig2.subplots_adjust(hspace=0,wspace=0)
    fig1.suptitle('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    fig2.suptitle('t=%.2e yr / nstep=%012d' % (snapshot.tout, snapshot.nstep))
    fig1.savefig('aex_%012d.png' % snapshot.nstep)
    fig2.savefig('aix_%012d.png' % snapshot.nstep)
    plt.close(fig1); plt.close(fig2)
print "// (%s UTC) Done Processing Snapshots" % strftime("%H:%M:%S", gmtime())
print "// Done"
