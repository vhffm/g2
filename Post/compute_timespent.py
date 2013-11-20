"""
Compute the time particles spend at a given semi-major axis, and eccentricity.
Assumes all time between the current and previous outputs was spent there.
"""

from glob import glob
import numpy as np
from Structs import Particle
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--custom', type=int, \
                   help="Plot Custom Snapshot.")
parser.add_argument('--scale', action='store_true', \
                    help="Scale Marker Size with Particle Mass")
args = parser.parse_args()

# Build Snapshot Number Array (From First Dir)
print "// Building Snapshot Array"
if args.all:
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    nsteps = np.zeros(len(globs))
    for ii, gg in enumerate(globs):
        nsteps[ii] = int(gg.split('.npz')[0].split('/')[-1].split('_')[1])
if args.test:
    nsteps = np.mgrid[3600000000:3630000000:1000000]
if args.custom:
    nsteps = np.array([args.custom])
print "// Found %i Snapshots" % len(nsteps)

# Determine Maximum Number of Particles
npz = np.load('Snapshot_%012d.npz' % nsteps[0])
snapshot = npz['snapshot'][()]
npartmax = snapshot.nparticles
print "// Found %i Particles" % npartmax 

# Loop Snapshots
# abins_limits = [ 0 1 2 3 4 5 6 ]
# abins        = [  X X X X X X  ]
bin_edges_a = np.linspace(0, 6, 32)
bin_edges_e = np.linspace(0, 0.1, 32)
abins = np.zeros([bin_edges_a.shape[0] - 1, npartmax])
ebins = np.zeros([bin_edges_e.shape[0] - 1, npartmax])
for istep, nstep in enumerate(nsteps):
    if nstep > 0:
        print "// Processing Snapshot %012d/%012d" % (nstep, nsteps[-1])
        npz = np.load('Snapshot_%012d.npz' % nstep)
        snapshot = npz['snapshot'][()]
        # Bin Particles
        for p in snapshot.particles:
            # Loop Semi-Major Axis Bins
            for ibin in range(bin_edges_a.shape[0]-1):
                if bin_edges_a[ibin] < p.a and p.a <= bin_edges_a[ibin+1]:
                    # Time passed since last step is
                    # (nstep-nsteps[istep-1])*6./365.25
                    abins[ibin,p.id] += (nstep-nsteps[istep-1])*6./365.25
                    break
            # Loop Eccentricity Bins
            for ibin in range(bin_edges_e.shape[0]-1):
                if bin_edges_e[ibin] < p.ecc and p.ecc <= bin_edges_e[ibin+1]:
                    ebins[ibin,p.id] += (nstep-nsteps[istep-1])*6./365.25
                    break

print "// Saving Data"
np.savez('TimeSpent.npz', \
    bin_edges_a = bin_edges_a, \
    bin_edges_e = bin_edges_e, \
    npartmax = npartmax, \
    abins = abins, \
    ebins = ebins )

print "// Done"
