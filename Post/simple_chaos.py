"""
Compute
- Separation (ds)
- Lyapunov Exponent (lce)
- Lyapynov Time (ltime)

Simple Way.
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
parser.add_argument('--outfile', default='XChaos.npz', \
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

# Allocate some arrays
tout = np.zeros(len(nsteps))

x1 = np.zeros([len(nsteps),nparts])
x2 = np.zeros_like(x1)

y1 = np.zeros_like(x1)
y2 = np.zeros_like(x1)

z1 = np.zeros_like(x1)
z2 = np.zeros_like(x1)

m1 = np.zeros_like(x1)
m2 = np.zeros_like(x1)

a1 = np.zeros_like(x1)
a2 = np.zeros_like(x1)

phase1 = np.zeros_like(x1)
phase2 = np.zeros_like(x2)

eccen1 = np.zeros_like(x1)
eccen2 = np.zeros_like(x1)

Omega1 = np.zeros_like(x1)
Omega2 = np.zeros_like(x1)

omega1 = np.zeros_like(x1)
omega2 = np.zeros_like(x1)

incli1 = np.zeros_like(x1)
incli2 = np.zeros_like(x1)

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
    # Make Sure Output Times Align
    if snap1.tout != snap2.tout:
        print "!! Output Times Misaligned"
        sys.exit()
    tout_loc = snap1.tout
    x1_loc = np.zeros(nparts)
    x2_loc = np.zeros_like(x1_loc)
    y1_loc = np.zeros_like(x1_loc)
    y2_loc = np.zeros_like(x1_loc)
    z1_loc = np.zeros_like(x1_loc)
    z2_loc = np.zeros_like(x1_loc)
    i1_loc = np.ones_like(x1_loc, dtype=int) * -1
    i2_loc = np.ones_like(x1_loc, dtype=int) * -1
    m1_loc = np.zeros_like(x1_loc)
    m2_loc = np.zeros_like(x1_loc)
    a1_loc = np.zeros_like(x1_loc)
    a2_loc = np.zeros_like(x1_loc)
    phase1_loc = np.zeros_like(x1_loc)
    phase2_loc = np.zeros_like(x1_loc)
    eccen1_loc = np.zeros_like(x1_loc)
    eccen2_loc = np.zeros_like(x1_loc)
    Omega1_loc = np.zeros_like(x1_loc)
    Omega2_loc = np.zeros_like(x1_loc)
    omega1_loc = np.zeros_like(x1_loc)
    omega2_loc = np.zeros_like(x1_loc)
    incli1_loc = np.zeros_like(x1_loc)
    incli2_loc = np.zeros_like(x1_loc)
    ireduce = 0
    for iparticle, particle in enumerate(snap1.particles):
        iparticle -= ireduce
        # Skip Particles in Ignore List
        if not int(particle.id) in ignore:
            i1_loc[iparticle] = int(particle.id)
            x1_loc[iparticle] = particle.x
            y1_loc[iparticle] = particle.y
            z1_loc[iparticle] = particle.z
            m1_loc[iparticle] = particle.m
            a1_loc[iparticle] = particle.a
            phase1_loc[iparticle] = particle.M0
            eccen1_loc[iparticle] = particle.ecc
            Omega1_loc[iparticle] = particle.Omega
            omega1_loc[iparticle] = particle.omega
            incli1_loc[iparticle] = particle.inc
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
            m2_loc[iparticle] = particle.m
            a2_loc[iparticle] = particle.a
            phase2_loc[iparticle] = particle.M0
            eccen2_loc[iparticle] = particle.ecc
            Omega2_loc[iparticle] = particle.Omega
            omega2_loc[iparticle] = particle.omega
            incli2_loc[iparticle] = particle.inc
        else:
            ireduce += 1
    # Fix Order 01
    x1_loc = x1_loc[i1_loc.argsort()]
    y1_loc = y1_loc[i1_loc.argsort()]
    z1_loc = z1_loc[i1_loc.argsort()]
    m1_loc = m1_loc[i1_loc.argsort()]
    a1_loc = a1_loc[i1_loc.argsort()]
    phase1_loc = phase1_loc[i1_loc.argsort()]
    eccen1_loc = eccen1_loc[i1_loc.argsort()]
    Omega1_loc = Omega1_loc[i1_loc.argsort()]
    omega1_loc = omega1_loc[i1_loc.argsort()]
    incli1_loc = incli1_loc[i1_loc.argsort()]
    i1_loc = i1_loc[i1_loc.argsort()]
    # Fix Order 02
    x2_loc = x2_loc[i2_loc.argsort()]
    y2_loc = y2_loc[i2_loc.argsort()]
    z2_loc = z2_loc[i2_loc.argsort()]
    m2_loc = m2_loc[i2_loc.argsort()]
    a2_loc = a2_loc[i2_loc.argsort()]
    phase2_loc = phase2_loc[i2_loc.argsort()]
    eccen2_loc = eccen2_loc[i2_loc.argsort()]
    Omega2_loc = Omega2_loc[i2_loc.argsort()]
    omega2_loc = omega2_loc[i2_loc.argsort()]
    incli2_loc = incli2_loc[i2_loc.argsort()]
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
    # Append
    tout[istep] = tout_loc
    x1[istep,:] = x1_loc; y1[istep,:] = y1_loc; z1[istep,:] = z1_loc; m1[istep,:] = m1_loc
    phase1[istep,:] = phase1_loc
    eccen1[istep,:] = eccen1_loc
    Omega1[istep,:] = Omega1_loc
    omega1[istep,:] = omega1_loc
    incli1[istep,:] = incli1_loc
    x2[istep,:] = x2_loc; y2[istep,:] = y2_loc; z2[istep,:] = z2_loc; m2[istep,:] = m2_loc
    phase2[istep,:] = phase2_loc
    eccen2[istep,:] = eccen2_loc
    Omega2[istep,:] = Omega2_loc
    omega2[istep,:] = omega2_loc
    incli2[istep,:] = incli2_loc
    a1[istep,:] = a1_loc; a2[istep,:] = a2_loc
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())

# Shift, and tile time array
tt = tout-tout[0]
tt = np.tile(tt,[nparts,1]).T

# Compute separations
print "// Computing Separations"
dx2=(x1-x2)**2.
dy2=(y1-y2)**2.
dz2=(z1-z2)**2.
ds2=dx2+dy2+dz2
ds=np.sqrt(ds2)

# Compute Lyapunov time
print "// Computing Lyapunov Times"
alpha=ds[1:,:]/ds[1,:]
# lce=1.0/tt[1:,:]*np.log(alpha)
lce=1.0/tt[2:,:]*np.log(alpha[1:])
ltime=1.0/lce

# Create semi-major axis bins
print "// Generating Bins in Semi-Major Axis"
bin_edges = np.linspace(0,6.0,13)
bin_cents = (bin_edges[1:]+bin_edges[:-1])/2.0
digitalism = np.zeros([len(nsteps), nparts], dtype=int)
for istep, nstep in enumerate(nsteps):
    digitalism[istep,:] = np.digitize(a1[istep,:], bin_edges)

# --------------
# Bin Separation
# --------------
# For details, see XChaos-PaintingByNumbers Notebook
print "// Binning Separations"

ds_mean = np.zeros([len(nsteps),len(bin_cents)])
ds_std = np.zeros_like(ds_mean)
ds_median = np.zeros_like(ds_mean)
ds_q1 = np.zeros_like(ds_mean)
ds_q3 = np.zeros_like(ds_mean)

# 0 - outer loop over all timesteps
for istep, nstep in enumerate(nsteps):
    # in each bin, allocate empty lists to store all values in this list
    ds_binned = [[] for _ in range(0,len(bin_cents))]
    # 1 - inner loop over all particles per timestep
    for ipart in range(0,nparts):
        # for relevant particle @ timestep
        # identify bin, append value to array above
        ds_binned[digitalism[istep,ipart]].append(ds[istep,ipart])
    # compute mean/med/std for each bin
    for ibin in range(0,len(bin_cents)):
        if len(ds_binned[ibin]) > 0:
            ds_mean[istep,ibin] = np.mean(ds_binned[ibin])
            ds_std[istep,ibin] = np.std(ds_binned[ibin])
            ds_median[istep,ibin] = np.median(ds_binned[ibin])
            ds_q1[istep,ibin] = np.percentile(ds_binned[ibin], 25)
            ds_q3[istep,ibin] = np.percentile(ds_binned[ibin], 75)
        else:
            ds_mean[istep,ibin] = np.nan
            ds_std[istep,ibin] = np.nan
            ds_median[istep,ibin] = np.nan
            ds_q1[istep,ibin] = np.nan
            ds_q3[istep,ibin] = np.nan

# -----------------
# Bin Lyapunov Time
# -----------------
# For details, see XChaos-PaintingByNumbers Notebook
print "// Binning Lyapunov Time"

# Skip first two steps
ltime_mean = np.zeros([len(nsteps)-2,len(bin_cents)])
ltime_std = np.zeros_like(ltime_mean)
ltime_median = np.zeros_like(ltime_mean)
ltime_q1 = np.zeros_like(ltime_mean)
ltime_q3 = np.zeros_like(ltime_mean)

# 0 - outer loop over all timesteps
# for istep in range(1,2):
for istep, nstep in enumerate(nsteps):
    # Skip first two steps!
    if istep > 1:
        # in each bin, allocate empty lists to store all values in this list
        ltime_binned = [[] for _ in range(0,len(bin_cents))]
        # 1 - inner loop over all particles per timestep
        for ipart in range(0,nparts):
            # for relevant particle @ timestep
            # identify bin, append value to array above
            ltime_binned[digitalism[istep,ipart]].append(ltime[istep-2,ipart])
        # compute mean/med/std for each bin
        for ibin in range(0,len(bin_cents)):
            if len(ltime_binned[ibin]) > 0:
                ltime_median[istep-2,ibin] = np.median(ltime_binned[ibin])
                ltime_q1[istep-2,ibin] = np.percentile(ltime_binned[ibin], 25)
                ltime_q3[istep-2,ibin] = np.percentile(ltime_binned[ibin], 75)
            else:
                ltime_median[istep-2,ibin] = np.nan
                ltime_q1[istep-2,ibin] = np.nan
                ltime_q3[istep-2,ibin] = np.nan

# Save Data
print "// Saving Data"
np.savez(args.outfile, \
    tt = tt, \
    t0 = tout[0], \
    ds = ds, \
    lce = lce, \
    ltime = ltime, \
    x1 = x1, y1 = y1, z1 = z1, \
    x2 = x2, y2 = y2, z2 = z2, \
    a1 = a1, m1 = m1, i1_loc = i1_loc, \
    a2 = a2, m2 = m2, i2_loc = i2_loc, \
    phase1 = phase1, ecc1 = eccen1, inc1 = incli1, \
    phase2 = phase2, ecc2 = eccen2, inc2 = incli2, \
    Omega1 = Omega1, omega1 = omega1, \
    Omega2 = Omega2, omega2 = omega2, \
    a0_bin_edges = bin_edges, \
    a0_bin_cents = bin_cents, \
    ds_mean = ds_mean, ds_std = ds_std, \
    ds_median = ds_median, ds_q1 = ds_q1, ds_q3 = ds_q3, \
    ltime_mean = ltime_mean, ltime_std = ltime_std, \
    ltime_median = ltime_median, ltime_q1 = ltime_q1, ltime_q3 = ltime_q3 )

# Done
print "// Done"
print ""
print ">> BEGIN NOTE <<"
print ""
print "Arrays lce, ltime, ltime_* skip first two time steps"
print "(a) ds=0 at zeroth step"
print "(b) lce=0 -> 1/lce = inf at first step"
print ""
print ">> END NOTE <<"
print ""
