"""
Modify Genga ICs. Scale Up Mass/Radius.
"""

import sys
import argparse

# Read Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--mass", action="store_true", help="Scale By Mass")
group.add_argument("--radius", action="store_true", help="Scale By Radius")
parser.add_argument("--factor", type=float, required=True, \
                    help="Scale Factor.")
args = parser.parse_args()

# Read Lines from Stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Loop Lines
lines_out = []
for line_in in lines_in:
    line = line_in.strip().split(" ")
    mass = float(line[2])
    radius = float(line[3])
    if args.radius:
        radius *= args.factor
        mass *= args.factor**(3.0)
    elif args.mass:
        radius *= args.factor**(1.0/3.0)
        mass *= args.factor
    line[2] = "%.16e" % mass
    line[3] = "%.16e" % radius
    lines_out.append(" ".join(line))

# Output Lines
for line_out in lines_out:
    print line_out
