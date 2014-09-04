import Loaders
from glob import glob
import numpy as np
import argparse
from time import gmtime, strftime

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--run_name", default='gasrun', \
                    help='Name of Simulation Run.')
parser.add_argument("--ellipses", action='store_true', \
                    help='Compute & Store Orbit Ellipses for Particles.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--all', action='store_true', \
                   help="Reduce Full Set of Snapshots.")
group1.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('--quad', action='store_true', \
                    help="Store Kepler Elements in Quad Precision.")
group2.add_argument('--double', action='store_true', \
                    help="Store Kepler Elements in Double Precision.")
args = parser.parse_args()

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Output set must be defined by three numbers."
        sys.exit()

# Precision?
if args.double:
    precision = 'double'
if args.quad:
    precision = 'quad'

# Full Set
if args.all:
    nsteps = []
    globs = glob("Out%s_*.dat" % args.run_name)
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)

# Custom Set
if args.custom:
    # Build Output Number Array (From Input)
    nsteps = \
        np.mgrid[args.custom[0]:args.custom[1]+args.custom[2]:args.custom[2]]
    print "// Using Outputs %012d:%012d:%012d" % \
        ( args.custom[0], args.custom[1], args.custom[2] )

# Load, Reduce, Save
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])
    loader = Loaders.GengaOut(nstep, args.ellipses, args.run_name, precision)
    loader.load()
    np.savez('Snapshot_%012d.npz' % nstep, snapshot=loader.snapshot)
print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())
