import Loaders
from glob import glob
import numpy as np
import argparse
from time import gmtime, strftime

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--run_name", default='Out', \
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
    globs = glob("%s.*.dat" % args.run_name)
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[1])
        nsteps.append(nstep)

# Test Set
if args.test:
    nsteps = np.mgrid[70000000:73000000:100000]

# Load, Reduce, Save
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()) nstep, nsteps[-1])
    loader = Loaders.SSAscii(nstep, args.ellipses)
    loader.load()
    np.savez('Snapshot_%012d.npz' % nstep, snapshot=loader.snapshot)
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())
