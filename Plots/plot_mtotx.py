"""
Plot Total Disk Mass. Use List to Specifiy Directories.
Also Compute Mean, and Standard Deviation.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['lines.linewidth'] = 1.0
import matplotlib.pyplot as plt
from g2_helpers import twopi, MEarth
import os
from scipy import stats as sps
import sys
from copy import copy

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
parser.add_argument('--chopchop', action='store_true', \
                    help='Only Show Last 2 Dirs of Path in Legend')
args = parser.parse_args()

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

# Tweak Dirs For Legend
dirs_leg = copy(dirs)
if args.chopchop:
    for idir, dirchar in enumerate(dirs_leg):
        dirs_leg[idir] = "/".join(dirchar.split("/")[-2:])

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
print "// Found %i Snapshots" % len(nsteps)

# Loop Dirs, Snaps, Particles
tout = np.zeros_like(nsteps)
mass = np.zeros([len(dirs),nsteps.shape[0]])
for idir, dirchar in enumerate(dirs):
    print "// Processing Directory %s [%i/%i]" % (dirchar, idir+1, len(dirs))
    for istep, nstep in enumerate(nsteps):
        try:
            npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
            snapshot = npz['snapshot'][()]
            if idir == 0:
                tout[istep] = snapshot.tout
            for particle in snapshot.particles:
                mass[idir,istep] += particle.m
        except IOError:
            mass[idir,istep] = np.nan

# Compute Lower, Upper Bounds
print "// Computing Stats"
mass_mean = sps.nanmean(mass, axis=0)
mass_std = sps.stats.nanstd(mass, axis=0)
mass_lo = mass_mean - mass_std
mass_hi = mass_mean + mass_std

# Normalize to Earth Masses
mass      /= MEarth
mass_mean /= MEarth
mass_std  /= MEarth
mass_lo   /= MEarth
mass_hi   /= MEarth

print "// Rendering Figures"

# Plot Mean + Standard Deviation
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ax.fill_between(tout/1.0e6, mass_lo, mass_hi, alpha=0.5, lw=0.5)
ax.plot(tout/1.0e6, mass_mean, 'k', lw=0.5)
ax.grid()
ax.set_xlabel('t [Myr]')
ax.set_ylabel('M [M_Earth]')
ax.set_ylim([0, np.nanmax([np.nanmax(mass), 5])])
ax.set_title('Total Mass Disk Mass (Mean+Std)')
fig.savefig('mtot_stats.pdf')
plt.close()
plt.clf()

# Plot Mean Curves
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ttout = np.repeat(np.array([tout]), len(dirs), axis=0)
ax.plot(np.rot90(ttout/1.0e6), np.rot90(mass), lw=1.0)
ax.grid()
ax.set_xlabel('t [Myr]')
ax.set_ylabel('M [M_Earth]')
ax.set_ylim([0, np.nanmax([np.nanmax(mass), 5])])
ax.set_title('Total Disk Mass')
ax.legend(dirs_leg, prop={'size':'x-small'}, loc='best')
fig.savefig('mtot_all.pdf')
plt.close()
plt.clf()

print "// Done"
