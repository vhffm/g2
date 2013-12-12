"""
Remove moons from stdin (requires default Genga column ordering).
python killmoons2.py --ce 95000 < $output_to_change > $changed_output
"""

import numpy as np
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ce', type=int, required=True, \
                    help="Number of close encounters required for removal")
args = parser.parse_args()

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Remove particles with too many CEs
lines_out = []
for line in lines_in:
    if not int(line.strip().split()[19]) >= args.ce:
        lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
