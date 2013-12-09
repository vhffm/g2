"""
Remove Moons from ICs.
python nomoons.py --ic_old XXX --ic_new YYY --output_file ZZZ.
"""

import numpy as np
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ic_old', type=str, required=True, \
                    help="Old IC file with moons.")
parser.add_argument('--ic_new', type=str, required=True, \
                    help="New IC file without moons.")
parser.add_argument('--output_file', type=str, required=True, \
                    help="Output file to scan for CEs at step 1000. ")
args = parser.parse_args()

# Read output
# Decide which particle IDs to remove
print "// Reading output                 -- %s" % args.output_file
with open(args.output_file, "r") as f:
    lines = f.readlines()
    pids_to_remove = []
    for line in lines:
        line = line.strip().split()
        # Column 19 counts number of close encounters (CEs)
        # If at step 1000 we had 1000 CEs, we probably have a moon
        if int(line[19]) == 1000:
            pids_to_remove.append(int(line[1]))

# Read initial conditions file
# Filter out particles we should remove
print "// Reading old initial conditions -- %s" % args.ic_old
with open(args.ic_old, "r") as f:
    lines_in = f.readlines()
    lines_out = []
    for line in lines_in:
        # int(line.strip().split()[1]) contains the particle ID
        if not int(line.strip().split()[1]) in pids_to_remove:
            lines_out.append(line)

# Write cleaned up initial conditions file
print "// Writing new initial conditions -- %s" % args.ic_new
with open(args.ic_new, "w") as f:
    for line in lines_out:
        f.write(line)
