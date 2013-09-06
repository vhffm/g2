"""
Plot Number of Particles.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl; mpl.rcParams['lines.linewidth'] = 1.0
import matplotlib.pyplot as plt
from g2_helpers import twopi
import os

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--d1', help="First Directory.")
parser.add_argument('--d2', help="Second Directory.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
args = parser.parse_args()

# Sanity Check
if (args.d1 and not args.d2) or (not args.d1 and args.d2):
    raise Exception("Must Pass Two Directories or None.")

# Single or Dual Run?
if args.d1 and args.d2:
    run_type = 'double'
else:
    run_type = 'single'

# Build Absolute Paths
# @todo Check Existence
if run_type == 'double':
    dir1 = os.getcwd() + "/" + args.d1
    dir2 = os.getcwd() + "/" + args.d2

# Full Set
if args.all:
    if run_type == 'single':
        nsteps = []
        globs = glob("Snapshot_*.npz")
        globs = sorted(globs)
        for g in globs:
            nstep = int(g.split("_")[1].split(".")[0])
            nsteps.append(nstep)
    elif run_type == 'double':
        nsteps1 = []
        nsteps2 = []
        globs1 = glob(dir1 + "/" + "Snapshot_*.npz")
        globs2 = glob(dir2 + "/" + "Snapshot_*.npz")
        globs1 = sorted(globs1)
        globs2 = sorted(globs2)
        for g in globs1:
            nstep = int(g.split('.npz')[0].split('/')[-1].split('_')[1])
            nsteps1.append(nstep)
        for g in globs2:
            nstep = int(g.split('.npz')[0].split('/')[-1].split('_')[1])
            nsteps2.append(nstep)

# Test Set
if args.test:
    if run_type == 'single':
        nsteps = np.mgrid[70000000:73000000:100000]
        # nsteps = np.mgrid[3600000000:3630000000:1000000]
    elif run_type == 'double':
        nsteps1 = np.mgrid[70000000:73000000:100000]
        nsteps2 = nsteps1

if run_type == 'single':
    tout = np.zeros(len(nsteps)); npart = np.zeros(len(nsteps))
elif run_type == 'double':
    tout1 = np.zeros(len(nsteps1)); npart1 = np.zeros(len(nsteps1))
    tout2 = np.zeros(len(nsteps2)); npart2 = np.zeros(len(nsteps2))

if run_type == 'single':
    for istep, nstep in enumerate(nsteps):
        npz = np.load('Snapshot_%012d.npz' % nstep)
        snapshot = npz['snapshot'][()]
        tout[istep] = snapshot.tout
        npart[istep] = snapshot.nparticles
elif run_type == 'double':
    for istep, nstep in enumerate(nsteps1):
        npz = np.load(dir1 + '/' + 'Snapshot_%012d.npz' % nstep)
        snapshot = npz['snapshot'][()]
        tout1[istep] = snapshot.tout
        npart1[istep] = snapshot.nparticles
    for istep, nstep in enumerate(nsteps2):
        npz = np.load(dir2 + '/' + 'Snapshot_%012d.npz' % nstep)
        snapshot = npz['snapshot'][()]
        tout2[istep] = snapshot.tout
        npart2[istep] = snapshot.nparticles

# Plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
if run_type == 'single':
    ax.plot(tout/1.0e6, npart)
elif run_type == 'double':
    ax.plot(tout1/1.0e6, npart1, 'b', label=args.d1)
    ax.plot(tout2/1.0e6, npart2, 'r', label=args.d2)
    ax.legend(loc='best', prop={'size':'small'})
ax.grid(True)
ax.set_xlabel('t [Myr]', size='small')
ax.set_ylabel('N [-]', rotation='horizontal', size='small')
fig.savefig('npart.pdf')
plt.close()
plt.clf()
