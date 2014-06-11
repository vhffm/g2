"""
Read XChaos.npz. Compute Distances. Write Distances.
"""

import numpy as np
import argparse
import kepler_helpers as kh

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Output File Name.")
parser.add_argument('--outfile', default='XChaos_Distances.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Read
print "// Reading %s" % (args.infile,)
npz_in = np.load(args.infile)
ds = npz_in["ds"]
tt = npz_in["tt"]
x1 = npz_in["x1"]; x2 = npz_in["x2"]
y1 = npz_in["y1"]; y2 = npz_in["y2"]
z1 = npz_in["z1"]; z2 = npz_in["z2"]
vx1 = npz_in["vx1"]; vx2 = npz_in["vx2"]
vy1 = npz_in["vy1"]; vy2 = npz_in["vy2"]
vz1 = npz_in["vz1"]; vz2 = npz_in["vz2"]

# Compute Distances
print "// Computing New Distances"
ds2 = ds**2.0
dv2 = (vx1-vx2)**2.0 + (vy1-vy2)**2.0 + (vz1-vz2)**2.0
rho2 = kh.cart2metricX(x1, y1, z1, vx1, vy1, vz1, \
                       x2, y2, z2, vx2, vy2, vz2)

# Write
print "// Writing %s" % (args.outfile,)
np.savez(args.outfile, tt = tt, ds2 = ds2, rho2 = rho2, dv2 = dv2)

# Done
print "// Done"
