"""
Load Genga Ascii File. Save HDF5.
"""

import numpy as np
import kepler_helpers as kh
import h5py
import argparse
from glob import glob
from time import gmtime, strftime

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--run_name", default='_impacts_nice1_01', \
                    help='Name of Simulation Run.')
parser.add_argument("--ellipses", action='store_true', \
                    help='Compute & Store Orbit Ellipses for Particles.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--all', action='store_true', \
                   help="Reduce Full Set of Snapshots.")
group1.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# Output Simulation Name
print "// Run Name: %s" % args.run_name

# Ellipse Warning
if args.ellipses:
    print ""
    print "!!"
    print "!! Ellipse Computations Un-Optimized and Slow. For-Loops Suck."
    print "!!"
    print ""

# Supress Invalid Errora (e.g., Arccos[pi])
np.seterr(invalid="ignore")

# Full Set
if args.all:
    nsteps = []
    globs = glob("Out%s_*.dat" % args.run_name)
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)
    print "// Processing %i Outputs" % len(nsteps)

# Custom Set
if args.custom:
    # Build Output Number Array (From Input)
    nsteps = \
        np.mgrid[args.custom[0]:args.custom[1]+args.custom[2]:args.custom[2]]
    print "// Using Outputs %012d:%012d:%012d" % \
        ( args.custom[0], args.custom[1], args.custom[2] )

# Load, Reduce, Save
print "// Starting -- %s UTC" % strftime("%H:%M:%S", gmtime())
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    fname = "Out%s_%012d.dat" % (args.run_name, nstep)
    genga = np.loadtxt("%s" % fname, ndmin=2)

    # t i m r x y z vx vy vz Sx Sy Sz
    # 0 1 2 3 4 5 6  7  8  9 10 11 12
    tout = genga[:,0].astype("float64")
    pid = genga[:,1].astype("int")
    m = genga[:,2].astype("float64")
    x = genga[:,4].astype("float64")
    y = genga[:,5].astype("float64")
    z = genga[:,6].astype("float64")
    vx = genga[:,7].astype("float64")
    vy = genga[:,8].astype("float64")
    vz = genga[:,9].astype("float64")

    # Compute Keplerian Elements
    a, e, i, Omega, omega, M0 = \
        kh.cart2kepX(x, y, z, vx, vy, vz, m, central_mass=1.0)

    # Write File
    with h5py.File("Snapshot_%012d.hdf5" % nstep, "w") as f5:
        f5.create_dataset("particles/pid", data=pid)
        f5.create_dataset("particles/m", data=m)
        f5.create_dataset("particles/x", data=x)
        f5.create_dataset("particles/y", data=y)
        f5.create_dataset("particles/z", data=z)
        f5.create_dataset("particles/vx", data=vx)
        f5.create_dataset("particles/vy", data=vy)
        f5.create_dataset("particles/vz", data=vz)
        f5.create_dataset("particles/a", data=a)
        f5.create_dataset("particles/e", data=e)
        f5.create_dataset("particles/i", data=i)
        f5.create_dataset("particles/Omega", data=Omega)
        f5.create_dataset("particles/omega", data=omega)
        f5.create_dataset("particles/M0", data=M0)
        f5["particles/x"].attrs["units"] = "AU"
        f5["particles/y"].attrs["units"] = "AU"
        f5["particles/z"].attrs["units"] = "AU"
        f5["particles/vx"].attrs["units"] = "AU/yr [1yr=2pi]"
        f5["particles/vy"].attrs["units"] = "AU/yr [1yr=2pi]"
        f5["particles/vz"].attrs["units"] = "AU/yr [1yr=2pi]"
        f5["particles/m"].attrs["units"] = "Msun"
        f5.attrs["tout"] = tout[0]
        f5.attrs["tout_units"] = "yr"
        f5.attrs["run_name"] = args.run_name

        # Compute Ellipses?
        # @todo - Write vectorized version?
        if args.ellipses:
            for ipid, pidloc in enumerate(pid):
                xell, yell, zell = \
                    kh.compute_ellipse(a[ipid], e[ipid], i[ipid], \
                                       Omega[ipid], omega[ipid])
                f5.create_dataset("ellipses/%i/x" % pidloc, data=xell)
                f5.create_dataset("ellipses/%i/y" % pidloc, data=yell)
                f5.create_dataset("ellipses/%i/z" % pidloc, data=zell)
                f5["ellipses/%i/x" % pidloc].attrs["units"] = "AU"
                f5["ellipses/%i/y" % pidloc].attrs["units"] = "AU"
                f5["ellipses/%i/z" % pidloc].attrs["units"] = "AU"

print "// Done -- %s UTC" % strftime("%H:%M:%S", gmtime())
