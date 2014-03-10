"""
Reads output file (<stdin).
Modifies one particle.
Writes new output file (>stdout).

Usage:
python modpart.py --x --v --m --pid 123 --factor 1.1 < $infile > $outfile

Note1:
Perturbation factor is range (1.0 - args.factor, 1.0 + args.factor)

Note2:
If --m, perturbs mass.
If --x, perturbs positions.
If --v, perturbs velocities.
If --pid < 0, all particles are modified.

Format:
line[0]        - t
line[1]        - i
line[2]        - m
line[3]        - r
line[4,5,6]    - x,y,z
line[7,8,9]    - vx,vy,vz
line[10,11,12] - Sx,Sy,Xz
"""

import sys
import argparse
import numpy as np

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pid', type=int, \
                    help="Particle ID.", required=True)
parser.add_argument('--factor', type=float, \
                    help="Mass Adjustment Factor.", required=True)
parser.add_argument('--ic', action='store_true', \
                    help="Initial Conditions?")
parser.add_argument("--m", action="store_true", \
                    help="Perturb Mass")
parser.add_argument("--x", action="store_true", \
                    help="Perturb Position")
parser.add_argument("--v", action="store_true", \
                    help="Perturb Velocity")
args = parser.parse_args()

# Sanity Check
if not (args.m or args.x or args.v):
    raise Exception("Perturb what variable?")

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Perturb Selected Particle (ID>0)
# Perturb All Particles (ID<0)
lines_out = []
for line in lines_in:
    line = line.split()
    if args.pid >= 0:
        if int(line[1]) == args.pid:
            factor = np.random.uniform(1.0 - args.factor, 1.0 + args.factor, 7)
            if args.m:
                line[2] = str(float(line[2]) * factor[0])
            if args.x:
                line[4] = str(float(line[4]) * factor[1])
                line[5] = str(float(line[5]) * factor[2])
                line[6] = str(float(line[6]) * factor[3])
            if args.v:
                line[7] = str(float(line[7]) * factor[4])
                line[8] = str(float(line[8]) * factor[5])
                line[9] = str(float(line[9]) * factor[6])
    else:
        factor = np.random.uniform(1.0 - args.factor, 1.0 + args.factor, 7)
        if args.m:
            line[2] = str(float(line[2]) * factor[0])
        if args.x:
            line[4] = str(float(line[4]) * factor[1])
            line[5] = str(float(line[5]) * factor[2])
            line[6] = str(float(line[6]) * factor[3])
        if args.v:
            line[7] = str(float(line[7]) * factor[4])
            line[8] = str(float(line[8]) * factor[5])
            line[9] = str(float(line[9]) * factor[6])
    line = " ".join(line)
    # Append space for outputs.
    if not args.ic:
        line = line + " "
    lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
