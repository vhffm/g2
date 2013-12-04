"""
Compute Particle Spacing, Lyapuynov Characteristic Exponents.
"""

from glob import glob
import numpy as np
from Structs import Particle
import sys
from copy import copy
from chaos_helpers import co_intersection, compute_lyapunov
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
parser.add_argument('--scale', action='store_true', \
                    help="Scale Marker Size with Particle Mass")
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

# Only proceed if we have two directories
if len(dirs) != 2:
    print "// Must Pass 2 Directories."
    sys.exit()

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

# Determine maximum number of particles we can have
npartmax = 0
for idir, dirchar in enumerate(dirs):
    try:
        npz = np.load('%s/Snapshot_%012d.npz' % (dirchar, nsteps[0]))
        snapshot = npz['snapshot'][()]
        if idir == 0:
            npartmax = snapshot.nparticles
        else:
            if snapshot.nparticles > npartmax:
                npartmax = snapshot.nparticles
    except IOError:
        print "!! Could Not Open %s/Snapshot_%012d.npz" % \
              (dirs[idir], nstep)

print "// Found %i Particles." % npartmax 

# Loop steps
first = True
ds = np.zeros([nsteps.shape[0], npartmax])
istep0 = np.zeros(npartmax, dtype=int)
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    # Status
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    # Load Snapshots
    npz1 = np.load('%s/Snapshot_%012d.npz' % (dirs[0], nstep))
    npz2 = np.load('%s/Snapshot_%012d.npz' % (dirs[1], nstep))
    snap1 = npz1['snapshot'][()]
    snap2 = npz2['snapshot'][()]

    # Allocate Arrays
    c1pid = np.zeros(snap1.nparticles, dtype='int'); c1m = np.zeros(snap1.nparticles)
    c2pid = np.zeros(snap2.nparticles, dtype='int'); c2m = np.zeros(snap2.nparticles)
    c1x = np.zeros(snap1.nparticles)
    c1y = np.zeros(snap1.nparticles)
    c1z = np.zeros(snap1.nparticles)
    c2x = np.zeros(snap2.nparticles)
    c2y = np.zeros(snap2.nparticles)
    c2z = np.zeros(snap2.nparticles)

    # Allocate arrays for initial particle IDs
    if first:
        c1pid0 = np.zeros(snap1.nparticles, dtype='int')
        c2pid0 = np.zeros(snap1.nparticles, dtype='int')

    # Allocate arrays for initial semi-major axes
    if first:
        c1a0 = np.zeros(snap1.nparticles)
        c2a0 = np.zeros(snap1.nparticles)

    # Fill particle IDs from snapshots
    for ip, p in enumerate(snap1.particles):
        c1pid[ip] = p.id; c1x[ip] = p.x; c1y[ip] = p.y; c1z[ip] = p.z; c1m[ip] = p.m
        if first:
            c1a0[ip] = p.a
    for ip, p in enumerate(snap2.particles):
        c2pid[ip] = p.id; c2x[ip] = p.x; c2y[ip] = p.y; c2z[ip] = p.z; c2m[ip] = p.m
        if first:
            c2a0[ip] = p.a
    if first:
        o1pid = copy(c1pid); o1m = copy(c1m)
        o2pid = copy(c2pid); o2m = copy(c2m)
        c1pid0 = c1pid.copy()
        c2pid0 = c2pid.copy()
        ds[0,:] = np.zeros(npartmax)
        first = False
    else:
        # All particle IDs where the mass has not changed between old & current
        pid_co1 = co_intersection(c1pid, o1pid, c1m, o1m)
        pid_co2 = co_intersection(c2pid, o2pid, c2m, o2m)
        pid_co = np.intersect1d(pid_co1, pid_co2)
        # Return all particle IDs present in both current snapshots
        pid_both = np.intersect1d(c1pid, c2pid)
        # Return all particle IDs
        # (i) present in both snapshots
        # (ii) with unchanged mass
        pid_relevant = np.intersect1d(pid_both, pid_co)
        # Order by increasing particle ID
        pid_relevant = np.sort(pid_relevant)
        # Find indices corresponding to particle IDs
        c1idx = np.zeros_like(pid_relevant)
        c2idx = np.zeros_like(pid_relevant)
        for ii in range(pid_relevant.shape[0]):
            c1idx[ii] = int(np.where(c1pid == pid_relevant[ii])[0])
            c2idx[ii] = int(np.where(c2pid == pid_relevant[ii])[0])
        # Compute ds
        for ii in range(pid_relevant.shape[0]):
            ds[istep,c1idx[ii]] = \
                np.sqrt( \
                         (c1x[c1idx[ii]] - c2x[c2idx[ii]])**2. + \
                         (c1y[c1idx[ii]] - c2y[c2idx[ii]])**2. + \
                         (c1z[c1idx[ii]] - c2z[c2idx[ii]])**2. \
                       )
            if ds[istep-1,c1idx[ii]] == 0 and ds[istep,c1idx[ii]] > 0:
                istep0[c1idx[ii]] = istep
        # Copy current particle IDs to old particle IDs; dito for masses
        o1pid = copy(c1pid); o1m = copy(c1m)
        o2pid = copy(c2pid); o2m = copy(c2m)
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())

# Compute Lyapuynov Exponents
print "// Computing LCEs"
tout = 6.0 * nsteps / 365.25
lce = compute_lyapunov(ds, istep0, nsteps, tout)

# Save Relevant Arrays
print "// Saving Data"
np.savez('XChaos.npz', \
    lce = lce, ds = ds, istep0 = istep0, nsteps = nsteps, tout = tout, \
    c1pid0 = c1pid0, c2pid0 = c2pid0, c1a0 = c1a0, c2a0 = c2a0)

# Done
print "// Done"
