"""
Check which particles are lost between two given outputs.
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
                   help="Full Snapshot Range.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Custom Snapshot Range.")
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
    # Build Snapshot Number Array (From Input)
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
first_in_last = np.in1d(pids_first, pids_last)
# Print if we lost particles
if (~first_in_last).any():
    print "// Lost particles %s" % pids_first[~first_in_last]
else:
    print "// No particles lost"

print "// Done"
