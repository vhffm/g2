"""
Chop off energy file after a given step.

There's some rounding trickery involved,
so too steps too close together in time
should not be checked on. Sorry!

Usage:
python chopenergy.py --maxstep 1000000 --dt 6 < $infile > $outfile
"""

import sys
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--maxstep', type=int, \
                    help="Particle ID.", required=True)
parser.add_argument('--dt', type=int, \
                    help="Timestep [Days]", required=True)
args = parser.parse_args()

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Compute maximum timestep
maxstep = args.maxstep * args.dt / 365.25

# Change mass of selected particle
lines_out = []
last = False
for iline, line in enumerate(lines_in):
    if last:
        break
    nstep = int(round(float(line.split()[0])))
    if nstep == int(round(maxstep)):
        last = True
    lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
