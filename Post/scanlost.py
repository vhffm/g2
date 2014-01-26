"""
Scan which particles are lost between two given outputs.
"""

import numpy as np
import argparse
import sys
from glob import glob
from copy import copy

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# Generate Snapshot Array
if args.all:
    # Build Snapshot Number Array (From First Dir)
    print "// Building Snapshot Array"
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    nsteps = np.zeros(len(globs))
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.npz')[0].split('/')[-1].split('_')[1])
    print "// Found %i Snapshots" % len(nsteps)
if args.custom:
    # Build Snapshot Number Array (Form Input)
    nsteps = \
        np.mgrid[args.custom[0]:args.custom[1]+args.custom[2]:args.custom[2]]
    print "// Using Snapshots %012d:%012d:%012d" % \
        ( args.custom[0], args.custom[1], args.custom[2] )

# Load particle IDs for first snapshot
npz_first = np.load('Snapshot_%012d.npz' % nsteps[0])
snap_first = npz_first['snapshot'][()]
pids_first = np.zeros(snap_first.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_first.particles):
    pids_first[iparticle] = particle.id

# Load particle IDs for last snapshot
npz_last = np.load('Snapshot_%012d.npz' % nsteps[-1])
snap_last = npz_last['snapshot'][()]
pids_last = np.zeros(snap_last.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_last.particles):
    pids_last[iparticle] = particle.id

# Compare 
last_in_first = np.in1d(pids_last, pids_first)
# Print if we lost particles
if (~last_in_first).any():
    print "// Lost particles %s" % pid_last[~last_in_first]
else:
    print "// No particles lost"

print "// Done"
