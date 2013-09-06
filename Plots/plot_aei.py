"""
Plot Semi-Major Axis vs. Eccentricity, Inclination.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl; mpl.rcParams['lines.linewidth'] = 0.5
import matplotlib.pyplot as plt
from g2_helpers import mkline, twopi, MEarth
import time

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
args = parser.parse_args()

# Full Set
if args.all:
    nsteps = []
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split("_")[1].split(".")[0])
        nsteps.append(nstep)

# Test Set
if args.test:
    nsteps = np.mgrid[70000000:73000000:100000]
    # nsteps = np.mgrid[3600000000:3630000000:1000000]

# Scan Limits
print "// Scanning Limits"
first = True
for istep, nstep in enumerate(nsteps):
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    for particle in snapshot.particles:
        if first:
            amax = particle.a; amin = particle.a
            eccmax = particle.ecc; eccmin = particle.ecc
            incmax = particle.inc; incmin = particle.inc
            mmax = particle.m; mmin = particle.m
            first = False
        else:
            if particle.a > amax: amax = particle.a
            if particle.a < amin: amin = particle.a
            if particle.ecc > eccmax: eccmax = particle.ecc
            if particle.ecc < eccmin: eccmin = particle.ecc
            if particle.inc > incmax: incmax = particle.inc
            if particle.inc < incmin: incmin = particle.inc
            if particle.m > mmax: mmax = particle.m
            if particle.m < mmin: mmin = particle.m

# rad2deg
incmax *= 360.0/twopi
incmin *= 360.0/twopi

# Compute Line for Marker Size(Mass)
# m, n = mkline(5./2000., 1.0, 1.0, 36.0)
m, n = mkline(mmin, 1.0, mmax, 36.0)

# Generate Plots
for istep, nstep in enumerate(nsteps):
    print "// Plotting %i/%i" % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    first = True
    for particle in snapshot.particles:
        if first:
            mmax_snap = particle.m
            first = False
        else:
            if particle.m > mmax_snap: mmax_snap = particle.m

    # 2x1 @ 100dpi
    fig = plt.figure(figsize=(16.0, 6.0))

    #
    # (a,e)
    #
    ax = fig.add_subplot(1,2,1)
    for particle in snapshot.particles:
        h, = ax.plot(particle.a, particle.ecc, 'b.')
        h.set_markersize(particle.m * m + n)
    ax.grid(True)
    ax.set_xlim([0,amax])
    ax.set_ylim([0,eccmax])
    ax.set_xlabel('a [AU]', size='small')
    ax.set_ylabel('e [-]', rotation='horizontal', size='small')

    #
    # (a,i)
    #
    ax = fig.add_subplot(1,2,2)
    for particle in snapshot.particles:
        h, = ax.plot(particle.a, particle.inc, 'b.')
        h.set_markersize(particle.m * m + n)
    ax.grid(True)
    ax.set_xlim([0,amax])
    ax.set_ylim([0,incmax])
    ax.set_xlabel('a [AU]', size='small')
    ax.set_ylabel('i [deg]', rotation='horizontal', size='small')

    #
    # Global
    # 
    plt.suptitle('t=%.2e yr / nstep=%012d / nparticles=%05d / Mmax=%.2e Me' %  \
                 (snapshot.tout, snapshot.nstep, \
                  snapshot.nparticles, mmax / MEarth) )
    fig.savefig('aei_%012d.png' % snapshot.nstep)
    plt.close()
    plt.clf()
