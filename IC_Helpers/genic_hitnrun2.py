"""
Generate HitnRun Genga Initial Conditions from SPH Collision Outputs.

Requires SPH files to be in Tipsy format (we use Pynbody).

Requires skid to generate fragment groups.
Ex: "skid -tau 0.1 -msol 4.80443417e-8 -kpc 2.06470049e-13
          -hop -o skid-08870 < run1.08870"

Requires a dummy parameter files.
Ex: "dKpcUnit = 2.06701e-13
     dMsolUnit = 4.80438e-08"

* Default Fragment Density:
Moon Density: 3.34 g/cc
Earth       : 5.51 g/cc
Fragment    : 4.50 g/cc (Somewhere In Between)
"""

import numpy as np
import pandas as pd
import constants as C
import kepler_helpers as kh
import io_helpers as ioh
import pynbody as pb
import vector_helpers as vh
import ic_helpers as ih
import physics_helpers as ph
import sys
import argparse


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true", \
                    help="Verbose Mode")
parser.add_argument("--nodump", action="store_true", \
                    help="Do not dump ICs to Stdout")
parser.add_argument("--sph_file", default="run1.08870", \
                    help="SPH Source File (Tipsy Binary)")
parser.add_argument("--group_file", default="skid-08870.grp", \
                    help="Skid Group File (Tipsy Array)")
parser.add_argument("--rho", type=float, default=4.5, \
                    help="Fragment Density (Default: 4.5 g/cc).")
parser.add_argument("--unbound_group", type=int, default=0, \
                    help="Group ID of Unbound Particles")
parser.add_argument("--frag_group", type=int, default=-1, \
                    help="Group ID of Main Fragment")
parser.add_argument("--earth_group", type=int, default=2, \
                    help="Group ID of Earth (+ Moon)")
parser.add_argument("--phi", type=float, default=0, \
                    help="Rotate Earth Velocity Vector (Degree)")
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--cgs', action='store_true', \
                    help="CGS Units (Alexandre).")
group1.add_argument('--g1', action='store_true', \
                    help="G=1, Rearth Units (Christian).")
args = parser.parse_args()

# Info Output
if args.verbose:
    print ""
    print "* Skid Group File     : %s" % args.group_file
    print "* SPH Input File      : %s" % args.sph_file
    print "* Unbound Group       : %i" % args.unbound_group
    print "* Main Fragment Group : %i" % args.frag_group
    print "* Earth Group         : %s" % args.earth_group
    print "* Rotation Angle      : %.2f Degree" % args.phi
    if args.cgs:
        print "* Units               : CGS"
    elif args.g1:
        print "* Units               : G=1, Rearth, km/s, 9.5565e25 g"
    print ""

# Convert Rotation
args.phi *= C.d2r

###############################################################################
# LOAD DATA
###############################################################################
if args.verbose:
    print "// Loading SPH Data"

# Load SPH Data
s = pb.load("%s" % args.sph_file)

# Units
if args.cgs:
    s["pos"] /= 1000.0 * 1000.0 * C.Rearth # Rearth
    s["vel"] /= 1000.0 * 100.0             # km/s
    s["mass"] /= 1000.0                    # kg
    s["mass"] /= C.mearth                  # Earth masses
elif args.g1:
    s["mass"] *= 9.5565e25 # g
    s["mass"] *= 1000.0    # kg
    s["mass"] /= C.mearth  # Earth masses

# Load Groups from Skid
# >>> Tau should be a few times the softening <<<
# >>> Check softening with s["eps"] <<<
# >>> Apparently, doing this with 100 times(!) the softening produces 
#     the best groups <<<
# skid -tau 0.1 -msol 4.80443417e-8 -kpc 2.06470049e-13
#      -hop -o skid-08870 < run1.08870
group = np.genfromtxt("%s" % args.group_file, skiprows=1, dtype=np.int64)
unique_groups = np.unique(group)

###############################################################################
# MAKE DATAFRAME
###############################################################################

data = { "x": s["x"].view(type=np.ndarray), \
         "y": s["y"].view(type=np.ndarray), \
         "z": s["z"].view(type=np.ndarray), \
         "vx": s["vx"].view(type=np.ndarray), \
         "vy": s["vy"].view(type=np.ndarray), \
         "vz": s["vz"].view(type=np.ndarray), \
         "mass": s["mass"].view(type=np.ndarray), \
         "particle_id": \
            np.asarray(range(len(s["x"].view(type=np.ndarray))), \
                       dtype=np.float64) + 100000, \
         "group_id": group }
cols = [ "particle_id", "group_id", "mass", "x", "y", "z", "vx", "vy", "vz" ]

df = pd.DataFrame(data = data, columns = cols)

# Debug
# df.head()

###############################################################################
# COM/Velocity
###############################################################################
if args.verbose:
    print "// Computing Center of Mass"

# Compute Centre of Mass & Centre of Mass w/ Velocity, Stick into Dataframe
group_id = []; group_mass = []
com_x, com_y, com_z = [], [], []
com_vx, com_vy, com_vz = [], [], []
for ig, ng in enumerate(unique_groups):
    dfx = df[df.group_id == ng]

    # Center of Mass
    cx = np.sum(dfx["mass"]*dfx["x"])/np.sum(dfx["mass"])
    cy = np.sum(dfx["mass"]*dfx["y"])/np.sum(dfx["mass"])
    cz = np.sum(dfx["mass"]*dfx["z"])/np.sum(dfx["mass"])

    # Center of Mass Velocity
    cvx = np.sum(dfx["mass"]*dfx["vx"])/np.sum(dfx["mass"])
    cvy = np.sum(dfx["mass"]*dfx["vy"])/np.sum(dfx["mass"])
    cvz = np.sum(dfx["mass"]*dfx["vz"])/np.sum(dfx["mass"])
    
    group_id.append(ng)
    group_mass.append(np.sum(dfx["mass"]))
    com_x.append(cx)
    com_y.append(cy)
    com_z.append(cz)
    com_vx.append(cvx)
    com_vy.append(cvy)
    com_vz.append(cvz)

# Prepare Dataframe Construction
group_id = np.asarray(group_id, dtype=np.int64)
group_mass = np.asarray(group_mass)
com_x = np.asarray(com_x)
com_y = np.asarray(com_y)
com_z = np.asarray(com_z)
com_vx = np.asarray(com_vx)
com_vy = np.asarray(com_vy)
com_vz = np.asarray(com_vz)

# Make Dataframe
data = { "group_id": group_id, "group_mass": group_mass, \
         "com_x": com_x, "com_y": com_y, "com_z": com_z, \
         "com_vx": com_vx, "com_vy": com_vy, "com_vz": com_vz }
cols = [ "group_id", "group_mass", \
         "com_x", "com_y", "com_z", \
         "com_vx", "com_vy", "com_vz" ]
df_groups = pd.DataFrame(data = data, columns = cols)

# Debug
# print df_groups

###############################################################################
# Shift to Earth Centric Frame
###############################################################################
if args.verbose:
    print "// Shifting Frame"

# Shift Individual Particles
df["x"] -= df_groups[df_groups.group_id == args.earth_group]["com_x"].iloc[0]
df["y"] -= df_groups[df_groups.group_id == args.earth_group]["com_y"].iloc[0]
df["z"] -= df_groups[df_groups.group_id == args.earth_group]["com_z"].iloc[0]
df["vx"] -= df_groups[df_groups.group_id == args.earth_group]["com_vx"].iloc[0]
df["vy"] -= df_groups[df_groups.group_id == args.earth_group]["com_vy"].iloc[0]
df["vz"] -= df_groups[df_groups.group_id == args.earth_group]["com_vz"].iloc[0]

# Shift Groups
df_groups["com_x"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_x"].iloc[0]
df_groups["com_y"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_y"].iloc[0]
df_groups["com_z"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_z"].iloc[0]
df_groups["com_vx"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_vx"].iloc[0]
df_groups["com_vy"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_vy"].iloc[0]
df_groups["com_vz"] -= \
    df_groups[df_groups.group_id == args.earth_group]["com_vz"].iloc[0]

###############################################################################
# Compute Escape Velocity wrt Earth
###############################################################################
if args.verbose:
    print "// Computing Escape Velocity"

# Individual Particles
df["dist"] = np.sqrt(df["x"]**2.0 + df["y"]**2.0 + df["z"]**2.0)
df["vrel"] = np.sqrt(df["vx"]**2.0 + df["vy"]**2.0 + df["vz"]**2.0)
df["vesc"] = np.sqrt( (2.0 * C.G * C.mearth) / (df["dist"] * C.Rearth) )
df["vrat"] = df["vrel"] / df["vesc"]

# Clumps
df_groups["dist"] = np.sqrt(df_groups["com_x"]**2.0 + \
                            df_groups["com_y"]**2.0 + \
                            df_groups["com_z"]**2.0)
df_groups["vrel"] = np.sqrt(df_groups["com_vx"]**2.0 + \
                            df_groups["com_vy"]**2.0 + \
                            df_groups["com_vz"]**2.0)
df_groups["vesc"] = np.sqrt( (2.0 * C.G * C.mearth) / \
                             (df_groups["dist"] * C.Rearth) )
df_groups["vrat"] = df_groups["vrel"] / df_groups["vesc"]

###############################################################################
# Load Solar ICs
###############################################################################
if args.verbose:
    print "// Loading Solar System ICs"

s2p, _ = ih.Solar2()
earth = s2p[2].split(" ")
xe = float(earth[4]) # AU
ye = float(earth[5])
ze = float(earth[6])
vxe = float(earth[7]) * C.genga_to_kms # km/s
vye = float(earth[8]) * C.genga_to_kms
vze = float(earth[9]) * C.genga_to_kms

if args.verbose:
    print "   Using Earth @"
    print "   (%.2f,%.2f,%.2f) AU" % (xe,ye,ze)
    print "   (%.2f,%.2f,%.2f) km/s" % (vxe,vye,vze)

###############################################################################
# IC File - Massive Fragment
###############################################################################
if args.verbose:
    print "// Generating Genga ICs for Main Fragment"

# Select Massive Fragment
dfx = df_groups.copy()
dfx = dfx[dfx.group_id == args.frag_group]

# Remove Below Escape Velocity
dfx = dfx[dfx.vrat>=1.0]

# Debug
# print dfx

# Rotate
dfx["com_x"], dfx["com_y"] = vh.rotate_xy(dfx["com_x"], \
                                          dfx["com_y"], \
                                          args.phi)
dfx["com_vx"], dfx["com_vy"] = vh.rotate_xy(dfx["com_vx"], \
                                            dfx["com_vy"], \
                                            args.phi)

# Shift
dfx["com_x"] *= C.Rearth/C.au2km
dfx["com_y"] *= C.Rearth/C.au2km
dfx["com_z"] *= C.Rearth/C.au2km
dfx["com_x"] += xe
dfx["com_y"] += ye
dfx["com_z"] += ze

# Boost
dfx["com_vx"] += vxe
dfx["com_vy"] += vye
dfx["com_vz"] += vze

# Debug
# dfx

# Generate Genga IC Line
# t i m r x y z vx vy vz Sx Sy Sz
lines_frag = []
ii = 10
for _, df_row in dfx.iterrows():
    line = ""
    line += "0.0 %06d " % ii
    line += "%.16e " % (df_row["group_mass"]*C.mearth/C.msun)
    line += "%.16e " % (ph.mass2radius(df_row["group_mass"]*C.mearth, \
                                       rho=4.5)/C.au2km)
    line += "%+.16e %+.16e %+.16e " % (df_row["com_x"], \
                                       df_row["com_y"], \
                                       df_row["com_z"])
    line += "%+.16e %+.16e %+.16e " % (df_row["com_vx"] * C.kms_to_genga, \
                                       df_row["com_vy"] * C.kms_to_genga, \
                                       df_row["com_vz"] * C.kms_to_genga)
    line += "0.0 0.0 0.0"
    lines_frag.append(line)
    ii += 1

###############################################################################
# IC File - Unbound Particles
###############################################################################
if args.verbose:
    print "// Generating Genga ICs for Unbound Particles"

# Select Unbound Group
dfx = df.copy()
dfx = dfx[dfx.group_id == args.unbound_group]

# Remove Below Escape Velocity
dfx = dfx[dfx.vrat>=1.0]

# Rotate
dfx["x"], dfx["y"] = vh.rotate_xy(dfx["x"], dfx["y"], args.phi)
dfx["vx"], dfx["vy"] = vh.rotate_xy(dfx["vx"], dfx["vy"], args.phi)

# Shift
dfx["x"] *= C.Rearth/C.au2km
dfx["y"] *= C.Rearth/C.au2km
dfx["z"] *= C.Rearth/C.au2km
dfx["x"] += xe
dfx["y"] += ye
dfx["z"] += ze

# Boost
dfx["vx"] += vxe
dfx["vy"] += vye
dfx["vz"] += vze

# Debug
# dfx
# dfx[dfx.group_id == frag_group]

# Generate Genga IC Lines
# t i m r x y z vx vy vz Sx Sy Sz
lines_unbound = []
for _, df_row in dfx.iterrows():
    line = ""
    line += "0.0 %06d " % df_row["particle_id"]
    line += "%.16e " % 0.0
    line += "%.16e " % (ph.mass2radius(df_row["mass"]*C.mearth, \
                                       rho=4.5)/C.au2km)
    line += "%+.16e %+.16e %+.16e " % (df_row["x"], df_row["y"], df_row["z"])
    line += "%+.16e %+.16e %+.16e " % (df_row["vx"] * C.kms_to_genga, \
                                       df_row["vy"] * C.kms_to_genga, \
                                       df_row["vz"] * C.kms_to_genga)
    line += "0.0 0.0 0.0"
    lines_unbound.append(line)

###############################################################################
# Write ICs
###############################################################################

if not args.nodump:

    # Output
    if args.verbose:
        print ""
        print "************ STARTING OUTPUT NOW ************"
        print ""

    # Write Solar System
    for planet in s2p:
        print planet

    # Write Fragments
    for line in lines_frag:
        print line

    # Write Unbound Particles
    for line in lines_unbound:
        print line

else:
    print "// Skipped Output (--nodump)"
    print "// Done"
    print ""
