"""
Modify Genga Outputs. Retain Subset of Particles.
"""

import sys
import numpy as np
import argparse

# Read Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--toic', action='store_true', \
                    help="Convert to Genga IC Format.")
parser.add_argument('--planets', action='store_true', \
                    help="Keep (Only) Planets.")
args = parser.parse_args()

# Read Lines from Stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Define Subset
if args.planets:
    keep_ids = \
        np.array([0,1,2,3,4,5,6,7,8,9])
else:
    keep_ids = \
        np.array( \
            [10279, 10132, 12164, 10326, 10393, 10254, 10455, 10471, 10394,
             10098, 10486, 10239, 10315, 10399, 10331, 10307, 10408, 10388,
             10281, 10321, 10324, 10477, 10090, 10360, 10385, 10294, 10401,
             10458, 10342, 10359, 10172, 10459, 10353, 10448, 10339, 10263,
             10428, 10320, 10521, 10478, 10290, 10499, 10362, 10427, 10201,
             10425, 10233, 10197, 10228, 10225, 10271, 10505, 10122, 10431,
             10343, 10210, 10184, 10364, 10224, 10479, 10165, 10211, 10348,
             10332, 10430, 10417, 10308, 10198, 10138, 10273, 10092, 10180,
             10469, 10247, 10484, 10089, 18643, 10310, 10169, 10283, 10170,
             10373, 12192, 10419, 10325, 10179, 10414, 10500, 10382, 10487,
             10289, 10267, 10391, 10363, 10463, 10260, 10293, 10411, 10441,
             10418, 10277, 10338, 10396, 10099, 10449, 10202, 10492, 10167,
             10227, 12186, 10405, 10432, 10403, 10127, 10282, 10317, 10174,
             10422, 10301, 10104, 10465, 10371, 10291, 10189, 10319, 10337,
             10269, 10166, 10334] \
             )

# Keep Particle Subset
lines_out = []
for line_in in lines_in:
    if int(line_in.strip().split(" ")[1]) in keep_ids:
        lines_out.append(line_in)

# Convert to IC Format
for iline, line in enumerate(lines_out):
    lines_out[iline] = " ".join(line.strip().split(" ")[:13])

# Output Lines
for line_out in lines_out:
    print line_out
