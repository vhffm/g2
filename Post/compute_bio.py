"""
Tabulate a(t) and e(t).
"""

from glob import glob
import numpy as np
from Structs import Particle
import argparse
from time import gmtime, strftime

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

# Build Snapshot Number Array (From First Dir)
print "// Building Snapshot Array"
if args.all:
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    nsteps = np.zeros(len(globs))
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.npz')[0].split('/')[-1].split('_')[1])
if args.test:
    nsteps = np.mgrid[3600000000:3630000000:1000000]
if args.custom:
    nsteps = np.array([args.custom])
print "// Found %i Snapshots" % len(nsteps)

# Determine Maximum Number of Particles
npz = np.load('Snapshot_%012d.npz' % nsteps[0])
snapshot = npz['snapshot'][()]
npartmax = snapshot.nparticles
print "// Found %i Particles" % npartmax 

# Loop Snapshots
a = np.ones([nsteps.shape[0], npartmax]) * np.nan
e = np.ones([nsteps.shape[0], npartmax]) * np.nan
m = np.ones([nsteps.shape[0], npartmax]) * np.nan
t = np.zeros(nsteps.shape[0])
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    t[istep] = snapshot.tout
    for p in snapshot.particles:
        a[istep,p.id] = p.a
        e[istep,p.id] = p.ecc
        m[istep,p.id] = p.m
    # Compute Apoapsis & Periapsis
    rp = (1.0 - e) * a
    ra = (1.0 + e) * a
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())

print "// Saving Data"
np.savez('XBio.npz', \
    t = t, \
    a = a, \
    e = e, \
    rp = rp, \
    ra = ra, \
    m = m )

print "// Done"
