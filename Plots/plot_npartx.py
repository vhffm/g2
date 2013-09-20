"""
Plot Number of Particles. Use List to Specifiy Directories.
Also Compute Mean, and Standard Deviation at Snapshots.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl; mpl.rcParams['lines.linewidth'] = 1.0
import matplotlib.pyplot as plt
from g2_helpers import twopi
import os
from scipy import stats as sps
import sys

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
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

# Build Absolute Paths
for ii, dd in enumerate(dirs):
    dirs[ii] = os.getcwd() + "/" + dd

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
npart = np.zeros([len(dirs),nsteps.shape[0]])
for idir, dirchar in enumerate(dirs):
    print "// Processing Directory %s [%i/%i]" % (dirchar, idir+1, len(dirs))
    for istep, nstep in enumerate(nsteps):
        try:
            npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nstep))
            snapshot = npz['snapshot'][()]
            if idir == 0:
                tout[istep] = snapshot.tout
            npart[idir, istep] = snapshot.nparticles
        except IOError:
            npart[idir, istep] = np.nan

# Compute Lower, Upper Bounds
print "// Computing Stats"
npart_mean = sps.nanmean(npart, axis=0)
npart_std = sps.stats.nanstd(npart, axis=0)
npart_lo = npart_mean - npart_std
npart_hi = npart_mean + npart_std

print "// Rendering Figures"

# Plot Mean + Standard Deviation
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ax.fill_between(tout/1.0e6, npart_lo, npart_hi, alpha=0.5, lw=0.5)
ax.plot(tout/1.0e6, npart_mean, 'k', lw=0.5)
ax.grid()
ax.set_xlabel('t [Myr]')
ax.set_ylabel('N', rotation='horizontal')
ax.set_ylim([0, np.nanmax(npart_hi)])
ax.set_title('Remaining Number of Particles (Mean+Std)')
fig.savefig('npartx_stats.pdf')
plt.close()
plt.clf()

# Plot Mean Curves
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ttout = np.repeat(np.array([tout]), len(dirs), axis=0)
ax.plot(np.rot90(ttout/1.0e6), np.rot90(npart), lw=0.5)
ax.grid()
ax.set_xlabel('t [Myr]')
ax.set_ylabel('N', rotation='horizontal')
ax.set_ylim([0, np.nanmax(npart)])
ax.set_title('Remaining Number of Particles')
fig.savefig('npartx_all.pdf')
plt.close()
plt.clf()

print "// Done"
