"""
Compute and Aggregate Statistics.
"""

import numpy as np
import argparse
import sys
from glob import glob
from time import gmtime, strftime
import io_helpers as ioh

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Aggregative Full Set of Snapshots.")
parser.add_argument('--maxout', type=int, \
                    help="Maximum Output Step.")
args = parser.parse_args()

# Masses in kg
m_earth   = 5.972e24
m_sun     = 1.99e30
m_jupiter = 1.89e27

# Cutoff Mass
m_cutoff = 2.0e23 # kg
m_cutoff = m_cutoff / m_earth

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

# Build Snapshot Number Array (From First Dir)
print "// Building Snapshot Array"
if args.all:
    globs = glob(dirs[0] + "/" + "Out_*.dat")
    globs = sorted(globs)
    nsteps = np.zeros(len(globs), dtype=np.int64)
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.dat')[0].split('/')[-1].split('_')[3])
        globs[ii] = gg.split("/")[-1]

# Remove Outputs > Max Output
if args.maxout:
    nsteps = nsteps[nsteps<=args.maxout]
    globs = globs[:len(nsteps)]

print "// Verifying Snapshots"
if args.all:
    xglobs = globs
    for ii, gg in enumerate(xglobs):
        xglobs[ii] = gg.split("/")[-1].split("_")[-1]
    for dir_loc in dirs:
        globs_loc = glob(dir_loc + "/" + "Out_%s_*.dat" % \
            dir_loc.split("/")[-1])
        globs_loc = sorted(globs_loc)
        # Remove Outputs > Max Output
        if args.maxout:
            globs_loc = globs_loc[:len(nsteps)]
        for ii, gg in enumerate(globs_loc):
            globs_loc[ii] = gg.split("/")[-1].split("_")[-1]
        if not (globs == globs_loc):
            print "!! Snapshots Differ. Terminating."
            print "!! %s" % dir_loc
            sys.exit(0)

# Loop Dirs, Snaps, Particles
tout = np.zeros_like(nsteps)
mass = np.zeros([len(dirs),nsteps.shape[0]])
mass_above_cutoff = np.zeros_like(mass)
mass_below_cutoff = np.zeros_like(mass)
npart = np.zeros([len(dirs),nsteps.shape[0]])
npart_above_cutoff = np.zeros_like(npart)
npart_below_cutoff = np.zeros_like(npart)
for idir, dirchar in enumerate(dirs):
    print "// (%s UTC) Processing Directory %s [%i/%i]" % \
        (strftime("%H:%M:%S", gmtime()), dirchar, idir+1, len(dirs))
    for istep, nstep in enumerate(nsteps):
        try:
            fname = "%s/Out_%s_%012d.dat" % (dirchar, \
                                             dirchar.split("/")[-1], \
                                             nstep)
            df = ioh.read_output(fname, frame="heliocentric")
            df = df[df.pid<2000]
            if idir == 0:
                tout[istep] = df.time.iloc[0]
            # Number of Particles
            npart[idir, istep] = len(df)
            npart_above_cutoff[idir,istep] = np.sum(df.mass >= m_cutoff)
            npart_below_cutoff[idir,istep] = np.sum(df.mass < m_cutoff)
            # Count Mass
            mass[idir,istep] = np.sum(df.mass)
            mass_above_cutoff[idir,istep] = \
                np.sum(df[df.mass >= m_cutoff].mass)
            mass_below_cutoff[idir,istep] = \
                np.sum(df[df.mass < m_cutoff].mass)
        except IOError:
            mass[idir,istep] = np.nan

# Process Mass
print "// Computing Mass Stats"
mass_avg = np.mean(mass, axis=0)
mass_med = np.median(mass, axis=0)
mass_std = np.std(mass, axis=0)
mass_q25 = np.percentile(mass, 25, axis=0)
mass_q75 = np.percentile(mass, 75, axis=0)
mass_q10 = np.percentile(mass, 10, axis=0)
mass_q90 = np.percentile(mass, 90, axis=0)

mass_above_cutoff_avg = np.mean(mass_above_cutoff, axis=0)
mass_above_cutoff_med = np.median(mass_above_cutoff, axis=0)
mass_above_cutoff_std = np.std(mass_above_cutoff, axis=0)
mass_above_cutoff_q25 = np.percentile(mass_above_cutoff, 25, axis=0)
mass_above_cutoff_q75 = np.percentile(mass_above_cutoff, 75, axis=0)
mass_above_cutoff_q10 = np.percentile(mass_above_cutoff, 10, axis=0)
mass_above_cutoff_q90 = np.percentile(mass_above_cutoff, 90, axis=0)

mass_below_cutoff_avg = np.mean(mass_below_cutoff, axis=0)
mass_below_cutoff_med = np.median(mass_below_cutoff, axis=0)
mass_below_cutoff_std = np.std(mass_below_cutoff, axis=0)
mass_below_cutoff_q25 = np.percentile(mass_below_cutoff, 25, axis=0)
mass_below_cutoff_q75 = np.percentile(mass_below_cutoff, 75, axis=0)
mass_below_cutoff_q10 = np.percentile(mass_below_cutoff, 10, axis=0)
mass_below_cutoff_q90 = np.percentile(mass_below_cutoff, 90, axis=0)

# Process Particle Number
print "// Computing Particle Number Stats"
npart_avg = np.mean(npart, axis=0)
npart_med = np.median(npart, axis=0)
npart_std = np.std(npart, axis=0)
npart_q25 = np.percentile(npart, 25, axis=0)
npart_q75 = np.percentile(npart, 75, axis=0)
npart_q10 = np.percentile(npart, 10, axis=0)
npart_q90 = np.percentile(npart, 90, axis=0)

npart_above_cutoff_avg = np.mean(npart_above_cutoff, axis=0)
npart_above_cutoff_med = np.median(npart_above_cutoff, axis=0)
npart_above_cutoff_std = np.std(npart_above_cutoff, axis=0)
npart_above_cutoff_q25 = np.percentile(npart_above_cutoff, 25, axis=0)
npart_above_cutoff_q75 = np.percentile(npart_above_cutoff, 75, axis=0)
npart_above_cutoff_q10 = np.percentile(npart_above_cutoff, 10, axis=0)
npart_above_cutoff_q90 = np.percentile(npart_above_cutoff, 90, axis=0)

npart_below_cutoff_avg = np.mean(npart_below_cutoff, axis=0)
npart_below_cutoff_med = np.median(npart_below_cutoff, axis=0)
npart_below_cutoff_std = np.std(npart_below_cutoff, axis=0)
npart_below_cutoff_q25 = np.percentile(npart_below_cutoff, 25, axis=0)
npart_below_cutoff_q75 = np.percentile(npart_below_cutoff, 75, axis=0)
npart_below_cutoff_q10 = np.percentile(npart_below_cutoff, 10, axis=0)
npart_below_cutoff_q90 = np.percentile(npart_below_cutoff, 90, axis=0)

# Store
print "// Saving"
np.savez("Stats.npz", \
    mass = mass, npart = npart, tout = tout, \
    mass_avg = mass_avg, mass_med = mass_med, mass_std = mass_std, \
    mass_above_cutoff = mass_above_cutoff, \
    mass_below_cutoff = mass_below_cutoff, \
    mass_above_cutoff_avg = mass_above_cutoff_avg, \
    mass_above_cutoff_med = mass_above_cutoff_med, \
    mass_above_cutoff_std = mass_above_cutoff_std, \
    mass_above_cutoff_q25 = mass_above_cutoff_q25, \
    mass_above_cutoff_q75 = mass_above_cutoff_q75, \
    mass_above_cutoff_q10 = mass_above_cutoff_q10, \
    mass_above_cutoff_q90 = mass_above_cutoff_q90, \
    mass_below_cutoff_avg = mass_below_cutoff_avg, \
    mass_below_cutoff_med = mass_below_cutoff_med, \
    mass_below_cutoff_std = mass_below_cutoff_std, \
    mass_below_cutoff_q25 = mass_below_cutoff_q25, \
    mass_below_cutoff_q75 = mass_below_cutoff_q75, \
    mass_below_cutoff_q10 = mass_below_cutoff_q10, \
    mass_below_cutoff_q90 = mass_below_cutoff_q90, \
    npart_avg = npart_avg, npart_med = npart_med, npart_std = npart_std, \
    npart_above_cutoff_avg = npart_above_cutoff_avg, \
    npart_above_cutoff_med = npart_above_cutoff_med, \
    npart_above_cutoff_std = npart_above_cutoff_std, \
    npart_below_cutoff_avg = npart_below_cutoff_avg, \
    npart_below_cutoff_med = npart_below_cutoff_med, \
    npart_below_cutoff_std = npart_below_cutoff_std, \
    mass_q25 = mass_q25, mass_q75 = mass_q75, \
    mass_q10 = mass_q10, mass_q90 = mass_q90, \
    npart_q25 = npart_q25, npart_q75 = npart_q75, \
    npart_q10 = npart_q10, npart_q90 = npart_q90 \
    )

# Done
print "// Done"
