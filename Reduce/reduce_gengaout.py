import Loaders
from glob import glob
import numpy as np
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--run_name", default='gasrun', \
                    help='Name of Simulation Run.')
parser.add_argument("--ellipses", action='store_true', \
                    help='Compute & Store Orbit Ellipses for Particles.')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Reduce Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Reduce Test Set of Snapshots.")
args = parser.parse_args()

# Full Set
if args.all:
    nsteps = []
    globs = glob("Out%s_*.dat" % args.run_name)
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[1])
        nsteps.append(nstep)

# Test Set
if args.test:
    nsteps = np.mgrid[3600000000:3630000000:1000000]

# Load, Reduce, Save
for istep, nstep in enumerate(nsteps):
    print "Reducing %i/%i." % (istep+1, len(nsteps))
    loader = Loaders.GengaOut(nstep, args.ellipses)
    loader.load()
    np.savez('Snapshot_%012d.npz' % nstep, snapshot=loader.snapshot)
