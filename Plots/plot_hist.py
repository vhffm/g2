"""
Plot Histograms of Semi-Major Axis, Eccentricity, Mean Anomaly.
@todo - For two simulation directories.
      - Superimpose lines!
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl; mpl.rcParams['lines.linewidth'] = 1.0
import matplotlib.pyplot as plt
from g2_helpers import twopi

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
print "Scanning Limits..."
first = True
for istep, nstep in enumerate(nsteps):
    # print "Scanning %i/%i" % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    for particle in snapshot.particles:
        if first:
            amax = particle.a; amin = particle.a
            eccmax = particle.ecc; eccmin = particle.ecc
            incmax = particle.inc; incmin = particle.inc
            first = False
        else:
            if particle.a > amax: amax = particle.a
            if particle.a < amin: amin = particle.a
            if particle.ecc > eccmax: eccmax = particle.ecc
            if particle.ecc < eccmin: eccmin = particle.ecc
            if particle.inc > incmax: incmax = particle.inc
            if particle.inc < incmin: incmin = particle.inc

# rad2deg
incmax *= 360.0/twopi
incmin *= 360.0/twopi

# Bins Edges
bin_edges_a = np.linspace(0, amax, 64)
bin_edges_ecc = np.linspace(0, eccmax, 64)
bin_edges_inc = np.linspace(0, incmax, 64)
bin_edges_M0 = np.linspace(0, 360, 64)

# Process Snapshots
for istep, nstep in enumerate(nsteps):
    print "Computing, Plotting %i/%i..." % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    pa = []; pecc = []; pM0 = []; pinc = []
    for particle in snapshot.particles:
        pa.append(particle.a)
        pecc.append(particle.ecc)
        pM0.append(particle.M0 * 360.0/twopi)
        pinc.append(particle.inc * 360.0/twopi)
    fig = plt.figure(figsize=(16.0, 12.0))
    #
    # Semi-Major Axis
    #
    ax = fig.add_subplot(2,2,1)
    hist_a, _ = \
        np.histogram(pa, bin_edges_a, \
                     weights=np.repeat(1./snapshot.nparticles, len(pa)))
    ax.step(bin_edges_a[:-1], hist_a, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,amax])
    ax.set_ylim([0,1.0])
    ax.set_xlabel('a [AU]')
    ax.set_ylabel('N/Ntotal [-]')
    #
    # Eccentricity
    #
    ax = fig.add_subplot(2,2,2)
    hist_ecc, _ = \
        np.histogram(pecc, bin_edges_ecc, \
                     weights=np.repeat(1./snapshot.nparticles, len(pecc)))
    ax.step(bin_edges_ecc[:-1], hist_ecc, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,eccmax])
    ax.set_ylim([0,1.0])
    ax.set_xlabel('e [AU]')
    ax.set_ylabel('N/Ntotal [-]')
    #
    # Mean Anomaly
    #
    ax = fig.add_subplot(2,2,3)
    hist_M0, _ = \
        np.histogram(pM0, bin_edges_M0, \
                     weights=np.repeat(1./snapshot.nparticles, len(pM0)))
    ax.step(bin_edges_M0[:-1], hist_M0, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,360])
    ax.set_ylim([0,1])
    ax.set_xlabel('Mean Anomaly [deg]')
    ax.set_ylabel('N/Ntotal [-]')
    #
    # Inclination
    #
    ax = fig.add_subplot(2,2,4)
    hist_inc, _ = \
        np.histogram(pinc, bin_edges_inc, \
                     weights=np.repeat(1./snapshot.nparticles, len(pinc)))
    ax.step(bin_edges_inc[:-1], hist_inc, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,incmax])
    ax.set_ylim([0,1])
    ax.set_xlabel('Inclination [deg]')
    ax.set_ylabel('N/Ntotal [-]')
    #
    # Main Title
    #
    fig.suptitle('t=%.2e [yr] / nstep=%010d / nparticles=%i' % \
                 (snapshot.tout, snapshot.nstep, snapshot.nparticles))
    #
    # Save Figure
    #
    fig.savefig('hist_%012d.png' % snapshot.nstep)
    plt.close()
    plt.clf()
