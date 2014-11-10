"""
Load and store timeseries for some particle.

Usage:
python ./mkseries --all --pid 123
"""

import numpy as np
import argparse
import sys
from time import gmtime, strftime
from glob import glob

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
parser.add_argument('--pid', type=int, required=True, \
                    help="Particle to store in time series.")
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

# Allocate arrays
x = np.ones(len(nsteps)) * np.nan
y = np.ones_like(x) * np.nan
z = np.ones_like(x) * np.nan
t = np.zeros_like(x)

# Load first snapshot, make sure particle is present
npz = np.load('Snapshot_%012d.npz' % nsteps[0])
snap = npz['snapshot'][()]
found = False
for iparticle, particle in enumerate(snap.particles):
    if int(particle.id) == args.pid:
        found = True

# Not Found?
if not found:
    print "// Particle Not Found. Quitting."
    sys.exit()

# Loop snapshots
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    # Status
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    # Load Snapshots
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snap = npz['snapshot'][()]
    t[istep] = snap.tout
    for iparticle, particle in enumerate(snap.particles):
        if int(particle.id) == args.pid:
            x[istep] = particle.x
            y[istep] = particle.y
            z[istep] = particle.z
            break
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())

# Save Data
print "// Saving Data"
np.savez("TimeSeries_%04d.npz" % args.pid, \
    x = x, y = y, z = z, t = t )

# Done
print "// Done"
