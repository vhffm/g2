"""
Scan which particles are lost from successive outputs.
"""

import numpy as np
import argparse
import sys
from glob import glob

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

# Loop snapshots
first = True
for istep, nstep in enumerate(nsteps):
    # Load Snapshot
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snap = npz['snapshot'][()]
    npart = snap.nparticles
    
    # Write particles IDs, then sort
    pid = np.zeros(npart, dtype=int)
    for iparticle, particle in enumerate(snap.particles):
        pid[iparticle] = particle.id

    # Compare
    if not first:
        last_in_now = np.in1d(pid_last, pid)
        # Print if we lost particles
        if (~last_in_now).any():
            print "// Step %i / Lost IDs %s" % (nstep, pid_last[~last_in_now])

    # Prepare next iteration
    pid_last = pid.copy()

    # Unset first
    if first:
        first = False

print "// Done"
