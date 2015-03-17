"""
Modify Genga Outputs. Convert Gas/Ice Giants to Solar2 Model.
"""

import ic_helpers as ih
import sys

# Read Lines from Stdin
lines = sys.stdin.read().rstrip("\n").split("\n")

# Load Solar System
# Relevant Indices 4, 5, 6, 7 - Jupiter, Saturn, Uranus, Neptune
solar2, _ = ih.Solar2()

# Modify Relevant Lines
for iline, line in enumerate(lines):
    line = line.strip().split()
    # Careful, in newer ICs particle IDs are increased by 1, i.e.
    # 5, 6, 7, 8 - Jupiter, Saturn, Uranus, Neptune
    pid = int(line[1])
    if pid in [ 4, 5, 6, 7, 8]:
        line_solar2 = solar2[pid].strip().split()
        # Update Position, Velocities
        for ii in [ 4, 5, 6, 7, 8, 9 ]:
            line[ii] = line_solar2[ii]
        lines[iline] = " ".join(line)

# Output Lines
for line in lines:
    print line
