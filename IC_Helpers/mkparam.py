"""
Make parameter file. Dump to stdout.

Usage:
python mkparam --XXX > param.dat
"""

import argparse
import sys
from sim_helpers import mkparam

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dt', type=float, \
                    help="Timestep [Days].", required=True)
parser.add_argument('--run_name', \
                    help="Run Name.", required=True)
parser.add_argument('--coordinate_interval', type=int, \
                    help="Coordinate Interval [Steps].", required=True)
parser.add_argument('--energy_interval', type=int, \
                    help="Energy Interval [Steps].", required=True)
parser.add_argument('--nsteps', type=int, \
                    help="Integration Steps From Here.", required=True)
parser.add_argument('--nrestart', type=int, \
                    help="Restart Step.", required=True)
args = parser.parse_args()

# Build output
lines = mkparam(dt = args.dt, \
                output_name = args.run_name, \
                energy_interval = args.energy_interval, \
                coordinate_interval = args.coordinate_interval, \
                restart_step = args.nrestart, \
                integration_steps = args.nsteps + args.nrestart)

# Dump to stdout
print "".join(lines)
