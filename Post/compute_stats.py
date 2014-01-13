"""
Compute and Aggregate Statistics.
"""

import numpy as np
import argparse
import sys
from glob import glob

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Aggregative Full Set of Snapshots.")
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

# Loop Dirs, Snaps, Particles
tout = np.zeros_like(nsteps)
mass = np.zeros([len(dirs),nsteps.shape[0]])
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
            for particle in snapshot.particles:
                # Do not count Jupiter and Saturn
                if particle.id != 2000 and particle.id != 2001:
                    mass[idir,istep] += particle.m
        except IOError:
            mass[idir,istep] = np.nan

# Process Mass
print "// Computing Mass Stats"
mass_avg = np.mean(mass, axis=0)
mass_med = np.median(mass, axis=0)
mass_std = np.std(mass, axis=0)

# Process Particle Number
print "// Computing Particle Number Stats"
npart_avg = np.mean(npart, axis=0)
npart_med = np.median(npart, axis=0)
npart_std = np.std(npart, axis=0)

# Store
print "// Saving"
np.savez("Stats.npz", \
    mass = mass, npart = npart, tout = tout, \
    mass_avg = mass_avg, mass_med = mass_med, mass_std = mass_std, \
    npart_avg = npart_avg, npart_med = npart_med, npart_std = npart_std \
    )

# Done
print "// Done"