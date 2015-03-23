"""
Modify Genga Outputs. Convert to Initial Conditions Format.
"""

import sys
import numpy as np

# Read Lines from Stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Convert to IC Format
lines_out = []
for iline, line in enumerate(lines_in):
    lines_out[iline] = " ".join(line.strip().split(" ")[:13])

# Output Lines
for line_out in lines_out:
    print line_out
