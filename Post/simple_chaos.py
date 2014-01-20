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
args = parser.parse_args()

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

# Scan number of particles
npz0 = np.load('%s/Snapshot_%012d.npz' % (dirs[0], nsteps[0]))
nparts = npz0['snapshot'][()].nparticles
print "// %i Particles Found" % nparts

# Compute time array
dt = 6.0
tout = nsteps * dt / 365.25
t0 = tout[0]

# Allocate some arrays
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
    i1_loc = np.zeros_like(x1_loc, dtype=int)
    i2_loc = np.zeros_like(x1_loc, dtype=int)
    m1_loc = np.zeros_like(x1_loc)
    m2_loc = np.zeros_like(x1_loc)
    a1_loc = np.zeros_like(x1_loc)
    a2_loc = np.zeros_like(x1_loc)
    for iparticle, particle in enumerate(snap1.particles):
        i1_loc[iparticle] = int(particle.id)
        x1_loc[iparticle] = particle.x
        y1_loc[iparticle] = particle.y
        z1_loc[iparticle] = particle.z
        m1_loc[iparticle] = particle.m
        a1_loc[iparticle] = particle.a
    for iparticle, particle in enumerate(snap2.particles):
        i2_loc[iparticle] = int(particle.id)
        x2_loc[iparticle] = particle.x
        y2_loc[iparticle] = particle.y
        z2_loc[iparticle] = particle.z
        m2_loc[iparticle] = particle.m
        a2_loc[iparticle] = particle.a
    # Safeguards
    for iidx, idx in enumerate(i1_loc):
        if i1_loc[iidx] != i2_loc[iidx]:
            print "!! Particle ID Mismatch"
            sys.exit()
    # Append
    x1[istep,:] = x1_loc; y1[istep,:] = y1_loc; z1[istep,:] = z1_loc; m1[istep,:] = m1_loc
    x2[istep,:] = x2_loc; y2[istep,:] = y2_loc; z2[istep,:] = z2_loc; m2[istep,:] = m2_loc
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
alpha=np.sqrt(ds2[1:,:]/ds2[1,:])
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
        ds_binned[digitalism[istep,ipart]].append(np.sqrt(ds2[istep,ipart]))
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
            ltime_binned[digitalism[istep,ipart]].append(1.0/lce[istep-2,ipart])
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
    t0 = t0, \
    ds = ds, \
    lce = lce, \
    ltime = ltime, \
    a1 = a1, m1 = m1, i1_loc = i1_loc, \
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
