"""
Fix Mass, Radius, Position, Velocity Format.
Fix Velocity Units On Request.

Only fix velocity units when we import from NASA HORIZONS or similar.

!!! Some Notes:

Genga & Pkdgrav use (almost) the same unit set.
G = M = 1 // 1yr = 2 Pi // 1 Day = 2 Pi / 365.25 ~ 0.017202423838958484.

The day in Genga is different from twpi/365.25 in the 6th decimal place. Why?
No adjustment is needed after running readinput.c for conversion.

Assume File Format:
<< t i m r x y z vx vy vz Sx Sy Sz >>
   0 1 2 3 4 5 6  7  8  9 10 11 12
"""

import sys
import argparse
import numpy as np

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--fixvel", action="store_true", \
                    help="Fix Velocity Units (AU/Day => AU/Day * Genga DayUnit")
args = parser.parse_args()

# Genga DayUnit
# Cf. define.h in genga/source
# Why is this not EXACTLTY = 2.0 * np.pi / 365.25 = 0.017202423838958484
dayUnit = 0.01720209895

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Loop lines. 
lines_out = []
for line in lines_in:
    line = line.split()
    line[2] = "%.16e" % float(line[2]) # Mass
    line[3] = "%.16e" % float(line[3]) # Radius
    for ii in [ 4, 5, 6 ]:
        line[ii] = "%+.16e" % float(line[ii]) # X,Y,Z
    for ii in [ 7, 8, 9 ]:
        if args.fixvel:
            line[ii] = "%+.16e" % (float(line[ii]) / dayUnit) # VX, VY, VZ
        else:
            line[ii] = "%+.16e" % float(line[ii]) # VX, VY, VZ
    line = " ".join(line)
    lines_out.append(line)

# Dump to stdout
for line in lines_out:
    print line
