"""
Read Full XChaos.npz.
Write Light XChaos.npz.
"""

import numpy as np
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='XChaos_Light.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Read
print "// Reading %s" % (args.infile,)
npz_in = np.load(args.infile)
ds = npz_in["ds"]
tt = npz_in["tt"]

# Write
print "// Writing %s" % (args.outfile,)
np.savez(args.outfile, ds = ds, tt = tt[:,0])

# Done
print "// Done"
