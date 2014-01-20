"""
Reads output file (<stdin).
Modifies one particle.
Writes new output file (>stdout).

Usage:
python modpart.py --pid 123 < $infile > $outfile
"""

import sys
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--pid', type=int, \
                    help="Particle ID.", required=True)
args = parser.parse_args()

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Change mass of selected particle
lines_out = []
for line in lines_in:
    line = line.split()
    if int(line[1]) == args.pid:
        line[2] = str(float(line[2])*1.01)
    line = " ".join(line) + " "
    lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
