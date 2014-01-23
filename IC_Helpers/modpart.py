"""
Reads output file (<stdin).
Modifies one particle.
Writes new output file (>stdout).

Usage:
python modpart.py --pid 123 --factor 1.1 < $infile > $outfile

Note:
If --pid < 0, all particles are modified.
If you set --milkshake, the entire box, and m,x,y,z,vx,vy,vz are modified.
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
parser.add_argument('--milkshake', action='store_true', \
                    help="Shake It, Baby! Requires --pid < 0.")
args = parser.parse_args()

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Change mass of selected particle
lines_out = []
for line in lines_in:
    line = line.split()
    if args.pid >= 0:
        if int(line[1]) == args.pid:
            line[2] = str(float(line[2]) * args.factor)
    else:
        if args.milkshake:
            # Milkshake mode
            # Modify m,x,y,z,vx,vy,vz by a random value
            # Values is between factor and 1-(1-factor).
            # Factor must be > 1.
            # E.g. - Factor = 1.1, 2 - Factor = 2 - 1.1 = 0.9
            factor = np.random.uniform(args.factor, 2-args.factor, 7)
            line[2] = str(float(line[2]) * factor[0])
            line[4] = str(float(line[4]) * factor[1])
            line[5] = str(float(line[5]) * factor[2])
            line[6] = str(float(line[6]) * factor[3])
            line[7] = str(float(line[7]) * factor[4])
            line[8] = str(float(line[8]) * factor[5])
            line[9] = str(float(line[9]) * factor[6])
        else:
            line[2] = str(float(line[2]) * args.factor)
    line = " ".join(line) + " "
    lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
