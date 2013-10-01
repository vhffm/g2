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
parser.add_argument('--weights', default='mass', choices=[ 'mass', 'npart' ], \
                    help="Histogram Weight Variable.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--t0', action='store_true', \
                   help="Plot First Snapshot.")
args = parser.parse_args()

# Some Info
if args.weights == 'mass':
    print "// Generating Mass Weighted Histograms"
elif args.weights == 'npart':
    print "// Generating Particle Number Histograms"

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

# Initial Conditions
if args.t0:
    nsteps = np.array([0])

# Scan Limits
print "// Scanning Range Limits"
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

# Histogram Limits
print "// Scanning Histogram Limits"
first = True
for istep, nstep in enumerate(nsteps):
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    pa = []; pecc = []; pinc = []
    pmass = np.zeros(snapshot.nparticles)

    for ii, particle in enumerate(snapshot.particles):
        pa.append(particle.a)
        pecc.append(particle.ecc)
        pinc.append(particle.inc * 360.0/twopi)
        pmass[ii] = particle.m

    if args.weights == 'mass':
        weights = pmass / float(np.sum(pmass))
    elif args.weights == 'npart':
        weights = np.repeat(1.0 / snapshot.nparticles, snapshot.nparticles)

    hist_a, _ = np.histogram(pa, bin_edges_a, weights = weights)
    hist_ecc, _ = np.histogram(pecc, bin_edges_ecc, weights = weights)
    hist_inc, _ = np.histogram(pinc, bin_edges_inc, weights = weights)

    if first:
        hist_a_max = np.nanmax(hist_a); hist_a_min = np.nanmin(hist_a)
        hist_ecc_max = np.nanmax(hist_ecc); hist_ecc_min = np.nanmax(hist_ecc)
        hist_inc_max = np.nanmax(hist_inc); hist_inc_min = np.nanmin(hist_inc)
        first = False
    else:
        if np.nanmax(hist_a) > hist_a_max: hist_a_max = np.nanmax(hist_a)
        if np.nanmin(hist_a) < hist_a_min: hist_a_min = np.nanmin(hist_a)
        if np.nanmax(hist_ecc) > hist_ecc_max: hist_ecc_max = np.nanmax(hist_ecc)
        if np.nanmin(hist_ecc) < hist_ecc_min: hist_ecc_min = np.nanmin(hist_ecc)
        if np.nanmax(hist_inc) > hist_inc_max: hist_inc_max = np.nanmax(hist_inc)
        if np.nanmin(hist_inc) < hist_inc_min: hist_inc_min = np.nanmin(hist_inc)

# Process Snapshots
for istep, nstep in enumerate(nsteps):
    print "// Computing, Plotting %i/%i" % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    pa = []; pecc = []; pinc = []
    pmass = np.zeros(snapshot.nparticles)

    for ii, particle in enumerate(snapshot.particles):
        pa.append(particle.a)
        pecc.append(particle.ecc)
        pinc.append(particle.inc * 360.0/twopi)
        pmass[ii] = particle.m

    if args.weights == 'mass':
        weights = pmass / float(np.sum(pmass))
        title = '(Mass Weighted)'
        ylabel = 'f_M [-]'
        fname = "hist_mass_%012d.png" % snapshot.nstep
    elif args.weights == 'npart':
        weights = np.repeat(1.0 / snapshot.nparticles, snapshot.nparticles)
        title = '(Particle Number Weighted)'
        ylabel = 'f_N [-]'
        fname = "hist_npart_%012d.png" % snapshot.nstep
    # fig = plt.figure(figsize=(16.0, 12.0))        # 2x2
    fig = plt.figure(figsize=(24.0, 6.0))           # 3x1

    #
    # Semi-Major Axis
    #
    ax = fig.add_subplot(1,3,1)
    hist_a, _ = np.histogram(pa, bin_edges_a, weights = weights)
    ax.step(bin_edges_a[:-1], hist_a, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,amax])
    ax.set_ylim([0,hist_a_max])
    ax.set_xlabel('a [AU]', size='small')
    ax.set_ylabel(ylabel, rotation='horizontal', size='small')
    ax.set_title(title, size='x-small')

    #
    # Eccentricity
    #
    ax = fig.add_subplot(1,3,2)
    hist_ecc, _ = np.histogram(pecc, bin_edges_ecc, weights = weights)
    ax.step(bin_edges_ecc[:-1], hist_ecc, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,eccmax])
    ax.set_ylim([0,hist_ecc_max])
    ax.set_xlabel('e [AU]', size='small')
    ax.set_ylabel(ylabel, rotation='horizontal', size='small')
    ax.set_title(title, size='x-small')

    #
    # Inclination
    #
    ax = fig.add_subplot(1,3,3)
    hist_inc, _ = np.histogram(pinc, bin_edges_inc, weights = weights)
    ax.step(bin_edges_inc[:-1], hist_inc, label='X')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0,incmax])
    ax.set_ylim([0,hist_inc_max])
    ax.set_xlabel('inc [deg]', size='small')
    ax.set_ylabel(ylabel, rotation='horizontal', size='small')
    ax.set_title(title, size='x-small')

    #
    # Main Title
    #
    fig.suptitle('t=%.2e [yr] / nstep=%010d / nparticles=%i' % \
                 (snapshot.tout, snapshot.nstep, snapshot.nparticles))

    #
    # Save Figure
    #
    fig.savefig(fname)
    plt.close()
    plt.clf()
