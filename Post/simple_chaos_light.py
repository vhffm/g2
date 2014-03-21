"""
Compute
- Separation (ds)
- Lyapunov Exponent (lce)
- Lyapynov Time (ltime)

Simple Way.

Light Way.
- Load, Save Only Separation (ds) & Semi-Major Axis (a1)
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
parser.add_argument('--outfile', default='XChaos_Light.npz', \
                    help="Output File Name.")
parser.add_argument('--ignore', type=int, nargs='+', \
                    help="Particle IDs to Ignore.")
parser.add_argument('-v', '--verbose', action="store_true", \
                    help="Show Verbose Output (Particle IDs).")
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('--genga', action='store_true', \
                    help="Genga Outputs.")
group2.add_argument('--pkd', action='store_true', \
                    help="Pkdgrav_Planet Outputs.")
args = parser.parse_args()

# Light Info
print "**"
print "** THIS IS THE LIGHT VERSION."
print "** THE OUTPUT FILE CONTAINTS ONLY THE FOLLOWING"
print "**"
print "** TIME,"
print "** SEPARATION,"
print "** SEMI-MAJOR AXIS."
print "**"

# Info On Output
if args.genga:
    print "// Assuming Genga Output"
elif args.pkd:
    print "// Assuming Pkdgrav_Planet Output"

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Snapshot set must be defined by three numbers."
        sys.exit()

# List of Directories
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
    print "!! Must Pass 2 Directories."
    sys.exit()

# Generate Snapshot Array
if args.all:
    # Build Snapshot Number Array (From First Dir)
    print "// Building Snapshot Array"
    globs = glob(dirs[0] + "/" + "Snapshot_*.npz")
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

#
# BEGIN IGNORE AUTO-DETECT
# 
# Detect particles we should ignore
# Adapted from check_lost.py
#

# Load first snapshot in first dir
npz_d0_first  = np.load('%s/Snapshot_%012d.npz' % (dirs[0], nsteps[0]))
snap_d0_first = npz_d0_first['snapshot'][()]
pids_d0_first = np.zeros(snap_d0_first.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_d0_first.particles):
    pids_d0_first[iparticle] = particle.id

# Load last snapshot in first dir
npz_d0_last  = np.load('%s/Snapshot_%012d.npz' % (dirs[0], nsteps[-1]))
snap_d0_last = npz_d0_last['snapshot'][()]
pids_d0_last = np.zeros(snap_d0_last.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_d0_last.particles):
    pids_d0_last[iparticle] = particle.id

# Load first snapshot in second dir
npz_d1_first  = np.load('%s/Snapshot_%012d.npz' % (dirs[1], nsteps[0]))
snap_d1_first = npz_d1_first['snapshot'][()]
pids_d1_first = np.zeros(snap_d1_first.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_d1_first.particles):
    pids_d1_first[iparticle] = particle.id

# Load last snapshot in second dir
npz_d1_last  = np.load('%s/Snapshot_%012d.npz' % (dirs[1], nsteps[-1]))
snap_d1_last = npz_d1_last['snapshot'][()]
pids_d1_last = np.zeros(snap_d1_last.nparticles, dtype=int)
for iparticle, particle in enumerate(snap_d1_last.particles):
    pids_d1_last[iparticle] = particle.id

# Compare first and last snapshot in first dir
first_in_last_mask_d0 = np.in1d(pids_d0_first, pids_d0_last)

# Compare first and last snapshot in second dir
first_in_last_mask_d1 = np.in1d(pids_d1_first, pids_d1_last)

# Make list of detected particles to ignore
ignore_auto = []
if (~first_in_last_mask_d0).any():
    print "// Lost %i Particle(s) in Dir 00" % \
        len(pids_d0_first[~first_in_last_mask_d0])
    ignore_auto.extend(pids_d0_first[~first_in_last_mask_d0])
    if args.verbose:
        print ""
        print pids_d0_first[~first_in_last_mask_d0]
        print ""
if (~first_in_last_mask_d1).any():
    print "// Lost %i Particle(s) in Dir 01" % \
        len(pids_d1_first[~first_in_last_mask_d1])
    ignore_auto.extend(pids_d1_first[~first_in_last_mask_d1])
    if args.verbose:
        print ""
        print pids_d1_first[~first_in_last_mask_d1]
        print ""

#
# END IGNORE AUTO-DETECT
# 

# Scan number of particles
npz0 = np.load('%s/Snapshot_%012d.npz' % (dirs[0], nsteps[0]))
snap0 = npz0['snapshot'][()]
found_saturn = False
found_jupiter = False
if args.genga:
    for particle in snap0.particles:
        if particle.id == 2000:
            found_jupiter = True
        if particle.id == 2001:
            found_saturn = True
nparts = npz0['snapshot'][()].nparticles
print "// %i Particles Found" % nparts
ignore_manual = []
if found_saturn and found_jupiter:
    print "// Found Jupiter and Saturn"
    ignore_manual.extend([2000])
    ignore_manual.extend([2001])
if args.ignore:
    ignore_manual.extend(args.ignore)

# Build Master Ignore Set
ignore = ignore_auto + ignore_manual
ignore = set(ignore)
nparts -= len(ignore)

# Some Info
print "// Auto Ignoring %i Particle(s)" % len(ignore_auto)
if args.verbose:
    print ""
    print ignore_auto
    print ""
print "// Manually Ignoring %i Particle(s)" % len(ignore_manual)
if args.verbose:
    print ""
    print ignore_manual
    print ""
print "// Ignoring %i Particle(s) in Total" % len(ignore)
if args.verbose:
    print ""
    print ignore
    print ""

# Compute time array
dt = 6.0
tout = nsteps * dt / 365.25
t0 = tout[0]

# Allocate some arrays
ds = np.zeros([len(nsteps),nparts])
a1 = np.zeros_like(ds)

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
    x1_loc = np.zeros(nparts)
    x2_loc = np.zeros_like(x1_loc)
    y1_loc = np.zeros_like(x1_loc)
    y2_loc = np.zeros_like(x1_loc)
    z1_loc = np.zeros_like(x1_loc)
    z2_loc = np.zeros_like(x1_loc)
    i1_loc = np.ones_like(x1_loc, dtype=int) * -1
    i2_loc = np.ones_like(x1_loc, dtype=int) * -1
    a1_loc = np.zeros_like(x1_loc)
    ireduce = 0
    for iparticle, particle in enumerate(snap1.particles):
        iparticle -= ireduce
        # Skip Particles in Ignore List
        if not int(particle.id) in ignore:
            i1_loc[iparticle] = int(particle.id)
            x1_loc[iparticle] = particle.x
            y1_loc[iparticle] = particle.y
            z1_loc[iparticle] = particle.z
            a1_loc[iparticle] = particle.a
        else:
            ireduce += 1
    ireduce = 0
    for iparticle, particle in enumerate(snap2.particles):
        iparticle -= ireduce
        # Skip Particles in Ignore List
        if not int(particle.id) in ignore:
            i2_loc[iparticle] = int(particle.id)
            x2_loc[iparticle] = particle.x
            y2_loc[iparticle] = particle.y
            z2_loc[iparticle] = particle.z
        else:
            ireduce += 1
    # Fix Order 01
    x1_loc = x1_loc[i1_loc.argsort()]
    y1_loc = y1_loc[i1_loc.argsort()]
    z1_loc = z1_loc[i1_loc.argsort()]
    a1_loc = a1_loc[i1_loc.argsort()]
    i1_loc = i1_loc[i1_loc.argsort()]
    # Fix Order 02
    x2_loc = x2_loc[i2_loc.argsort()]
    y2_loc = y2_loc[i2_loc.argsort()]
    z2_loc = z2_loc[i2_loc.argsort()]
    i2_loc = i2_loc[i2_loc.argsort()]
    # Safeguard 01
    if (i1_loc==-1).any() or (i2_loc==-1).any():
        print "!! Lost Particles"
        i1_in_i2 = np.in1d(i1_loc[i1_loc != -1], i2_loc[i2_loc != -1])
        i2_in_i1 = np.in1d(i2_loc[i2_loc != -1], i1_loc[i1_loc != -1])
        print "!! IDs in Snap01, Not in Snap 02 -- %s" % i1_loc[~i1_in_i2]
        print "!! IDs in Snap02, Not in Snap 01 -- %s"  % i2_loc[~i2_in_i1]
        sys.exit()
    # Safeguard 02
    for iidx, idx in enumerate(i1_loc):
        if i1_loc[iidx] != i2_loc[iidx]:
            print "!! Particle IDs Unaligned (%i != %i)" % \
                (i1_loc[iidx], i2_loc[iidx])
            sys.exit()
    # Compute Separation
    dx2_loc=(x1_loc-x2_loc)**2.
    dy2_loc=(y1_loc-y2_loc)**2.
    dz2_loc=(z1_loc-z2_loc)**2.
    ds2_loc=dx2_loc+dy2_loc+dz2_loc
    ds_loc=np.sqrt(ds2_loc)
    # Append
    ds[istep,:] = ds_loc
    a1[istep,:] = a1_loc
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())

# Shift, and tile time array
tt = tout-tout[0]
tt = np.tile(tt,[nparts,1]).T

# Save Data
print "// Saving Data"
np.savez(args.outfile, \
    tt = tt, \
    t0 = t0, \
    ds = ds, \
    a1 = a1 )

# Done
print "// Done"
