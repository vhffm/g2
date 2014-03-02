"""
Read Full XChaos.npz.
Write Light XChaos.npz.
With XYZ!
"""

import numpy as np
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='XChaos_Light_XYZ.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Read
print "// Reading %s" % (args.infile,)
npz_in = np.load(args.infile)
ds = npz_in["ds"]
a1 = npz_in["a1"]
tt = npz_in["tt"]
x1 = npz_in["x1"]
x2 = npz_in["x2"]
y1 = npz_in["y1"]
y2 = npz_in["y2"]
z1 = npz_in["z1"]
z2 = npz_in["z2"]

# Write
print "// Writing %s" % (args.outfile,)
np.savez(args.outfile, ds = ds, tt = tt[:,0], a1 = a1[0,:], \
    x1 = x1, y1 = y1, z1 = z1, \
    x2 = x2, y2 = y2, z2 = z2)

# Done
print "// Done"
