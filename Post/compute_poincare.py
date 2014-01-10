"""
Read 
"""

from glob import glob
import numpy as np
from Structs import Particle
import argparse
from time import gmtime, strftime
import sys

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--run_name", default='gasrun', \
                    help='Name of Simulation Run.')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# Build Snapshot Number Array (From First Dir)
print "// Building Snapshot Array"
if args.all:
    globs = glob("Poincare%s_*.dat" % args.run_name)
    globs = sorted(globs)
    nsteps = np.zeros(len(globs), dtype='int')
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.dat')[0].split('_')[-1])
if args.test:
    nsteps = np.mgrid[3600000000:3630000000:1000000]
if args.custom:
    nsteps = np.array([args.custom])
print "// Found %i Snapshots" % len(nsteps)

# Prescan to determine how many lines we need to allocate for data
nlines_tot = 0
for nstep in nsteps:
    fname = "Poincare%s_%012d.dat" % (args.run_name, nstep)
    with open(fname, 'r') as f:
        nlines_tot += len(f.readlines())

# Allocate data array
data = np.zeros([nlines_tot,4])

# @todo Need some sort of lookahead so we know how many Poincare crossings there will be?
# This sets the maximum size of our array
# We also need the maximum number of particles we are going to have
# So, let's loop!
iline_tot = 0
for nstep in nsteps:
    fname = "Poincare%s_%012d.dat" % (args.run_name, nstep)
    with open(fname, 'r') as f:
        lines = f.readlines()
        nlines_loc = len(lines)
        data_loc = np.zeros([nlines_loc,4])
        for iline, line in enumerate(lines):
            line = line.strip().split(' ')
            line = np.asarray(line, dtype='float64')
            data_loc[iline,:] = line
        # print data_loc
        data[iline_tot:iline_tot + nlines_loc,:] = data_loc
        iline_tot += nlines_loc
    
# Find largest ID, maximum number of crossings
data_pid = data[:,1].astype('int')
npart = np.max(data_pid)
ncrossings = np.max(np.bincount(data_pid))

print npart
print ncrossings

# print np.sort(data_pid.astype('int'))
# print np.sort(data_pid.astype('int'))
# print np.bincount(np.sort(data_pid.astype('int')))
# print np.max(np.bincount(np.sort(data_pid.astype('int'))))

# Load Pointcare Snaps
