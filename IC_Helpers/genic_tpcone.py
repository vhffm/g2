"""
Generate Test Particle Cone.

Replaces ID (--pid) with cone of Test Particles.
Their velocity vectors centered around original particle.
Massive Particle cloned as Test Particle w/ ID 10000.
Cone Test Particle have ID >= 10001.
"""

import numpy as np
import argparse
import sys

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--pid", type=int, required=True, \
                    help="Particle ID to Replace.")
parser.add_argument("--nr", type=int, default=32, \
                    help="Cone: Number of Radial Samples")
parser.add_argument("--nalpha", type=int, default=128, \
                    help="Cone: Number of Angular Samples")
parser.add_argument("--aspect", type=float, default=0.1, \
                    help="Cone: Aspect Ratio (Circle vs. Ellipse)")
args = parser.parse_args()

# Coding Convenience
nr = args.nr
nalpha = args.nalpha
aspect = args.aspect

# Read Lines from Stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Loop Lines
lines_out = []; replaced = False
for line in lines_in:
    pid = int(line.strip().split()[1])
    if pid == args.pid:
        # CheckMe
        replaced = True
        # Read Massive Particle
        line_loc = line.strip().split()
        r_in = float(line_loc[3])
        x_in = np.array([ float(line_loc[4]), \
                          float(line_loc[5]), \
                          float(line_loc[6]) ])
        v_in = np.array([ float(line_loc[7]), \
                          float(line_loc[8]), \
                          float(line_loc[9]) ])

        # Massive Particle => Test Particle (ID 10000)
        # line_new = "0.0 %i 0.0 %.16e " % (10000, r_in)
        # line_new += "%.16e %.16e %.16e " % (x_in[0], x_in[1], x_in[2])
        # line_new += "%.16e %.16e %.16e " % (v_in[0], v_in[1], v_in[2])
        # line_new += "0.0 0.0 0.0"
        # lines_out.append(line_new)
        lines_out.append(line)

        # Generate Cone Ranges
        xr, dr = np.linspace(np.pi/64.0, 0.0, nr, endpoint=False, retstep=True)
        xalpha = np.linspace(-np.pi, np.pi, nalpha, endpoint=False)
        xr = xr[::-1]; dr = -dr

        # Allocate Arrays
        vold = np.tile(v_in[np.newaxis,:], (nr*nalpha,1))
        vnew = np.zeros([nr*nalpha,3])

        # Generate Velocity Vectors
        # @todo - Clearly, this code can be vectorized, etc.
        #       - It would probably take an hour of coding for neglible gain.
        #       - So, meh.
        for ir, r in enumerate(xr):
            for ialpha, alpha in enumerate(xalpha):
                # Index Juggling
                ii = ir * nalpha + ialpha
                # XY in Cone (=Theta/Phi)
                theta = r * np.cos(alpha)
                phi = r * np.sin(alpha) * aspect
                # Velocity Normalization
                vnorm = np.sin(1.0 - (r / (np.pi/64.0 + dr)))
                # Rotations
                Rz = np.array([[ np.cos(theta), -np.sin(theta), 0.0 ], \
                               [ np.sin(theta),  np.cos(theta), 0.0 ], \
                               [ 0.0, 0.0, 1.0 ]])
                Ry = np.array([[ np.cos(phi), 0.0, np.sin(phi) ], \
                               [ 0.0, 1.0, 0.0 ], \
                               [ - np.sin(phi), 0.0, np.cos(phi) ] ])
                # Rotate, Scale, Replace
                vnew[ii,:] = np.dot(Ry, np.dot(Rz, vold[ii,:])) * vnorm

        # Generate IC Lines
        for ii in range(vnew.shape[0]):
            line_new = "0.0 %i 0.0 %.16e " % (ii + 10001, r_in)
            line_new += "%.16e %.16e %.16e " % (x_in[0], x_in[1], x_in[2])
            line_new += "%.16e %.16e %.16e " % (vnew[ii,0], \
                                                vnew[ii,1], \
                                                vnew[ii,2])
            line_new += "0.0 0.0 0.0"
            lines_out.append(line_new)

    else:
        lines_out.append(line)

# Removed Stats
if not replaced:
    sys.stderr.write("Particle %i Not Found. Nothing Replaced.\n" % args.pid)

# Stdout Dump
for line in lines_out:
    print line
