"""
Generate HitAndRun Initial Conditions.
Use Solar2 for Planets.
Generates Main Fragment and/or Cones of Test Particles.

If --cone and --fragment, generates Massive Fragment and Cone.
If only --cone, Main Fragment is also a Test Particle.

* Default Fragment Density:
Moon Density: 3.34 g/cc
Earth       : 5.51 g/cc
Fragment    : 4.50 g/cc (Somewhere In Between)

Debug Mode shows only Earth, Main Fragment, and 10 Test Particles.
"""

import numpy as np
import argparse
import ic_helpers as ih
import kepler_helpers as kh
import physics_helpers as ph
import constants as C
import sys


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true', \
                    help="Debug Mode.")
parser.add_argument("--rho", type=float, default=4.5, \
                    help="Fragment Density (Default: 4.5 g/cc).")
parser.add_argument("--r12sim", type=str, required=True, \
                    choices=["cC03p","cC03m","fA01p","fA01m" ], \
                    help="Select Reufer+ 2012 Simulation.")
group1 = parser.add_argument_group()
group1.add_argument('--fragment', action='store_true', \
                    help="Primary Fragment.")
group1.add_argument('--cone', action='store_true', \
                    help="Test Particle Cone .")
# group2 = parser.add_mutually_exclusive_group(required=True)
# group2.add_argument('--solar2', action='store_true', \
#                     help="Base Solar System: NASA Horizon Query 2014-01-01.")
args = parser.parse_args()

# Sanity Check
if not (args.fragment or args.cone):
    print "At least one of --fragment or --cone is required."
    sys.exit()

# Load Planets (Solar2)
ss, _ = ih.Solar2()
earth = ss[2]

# Load Fragment Orbital Elements
a, e, i, Omega, omega, M0, mass = ih.MainFragmentReufer12(args.r12sim, earth)

# Convert
x, v = kh.kep2cart(a, e, i, Omega, omega, M0, mass * C.mearth / C.msun)
r_large = ph.mass2radius(mass * C.mearth, rho=args.rho) / C.au2km
r_small = ph.mass2radius(mass/4096./100.0 * C.mearth, rho=args.rho) / C.au2km
line_fragment = "0.0 "
if args.fragment:
    line_fragment += "%05d " % 10
    line_fragment += "%.16e " % (mass * C.mearth / C.msun)
    line_fragment += "%.16e " % r_large
else:
    if args.cone:
        line_fragment += "%05d " % 10000
        line_fragment += "%.16e " % 0.0
        line_fragment += "%.16e " % r_small
line_fragment += "%+.16e %+.16e %+.16e " % (x[0], x[1], x[2])
line_fragment += "%+.16e %+.16e %+.16e " % (v[0], v[1], v[2])
line_fragment += "0.0 0.0 0.0"

# Generate Cone Test Particles
if args.cone:
    #lines_cone = ih.Cone(x, v, r_small, nr=2, nalpha=4)
    lines_cone = ih.Cone(x, v, r_small)

# Debug Mode
# Show Impactor Orbit, Earth-Impactor Distance, ID 2, and 10 Fragments
if args.debug:
    earth = earth.strip().split()
    xe = np.array([float(earth[4]), float(earth[5]), float(earth[6])]) # AU
    ds = np.sqrt(np.sum((x-xe)**2.0)) * C.au2km
    debug_impactor = "I(a,e,i,Omega,omega,M0,m) = "
    debug_impactor += "(%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)" % \
        (a, e, i, Omega, omega, M0, mass)
    print debug_impactor
    print "D(Earth-Impactor) = %.2f km" % ds

    # Print Earth, Main Fragment, First 10 Test particles
    print " ".join(earth)
    print line_fragment
    if args.cone:
        for itp, tp in enumerate(lines_cone):
            if itp < 10:
                print tp

# Default Mode
else:
    # Stdout (Planets)
    for planet in ss:
        print planet

    # Std (Main Fragments)
    print line_fragment

    # Cone?
    if args.cone:
        for tp in lines_cone:
            print tp
