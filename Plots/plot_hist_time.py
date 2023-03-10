"""
Plot Histograms w/ Time Axis.
"""

from glob import glob
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['image.interpolation'] = 'nearest'
import matplotlib.pyplot as plt
import os
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

# Scan Limits
print "// Scanning Limits"
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

# Bin Edges
bin_edges_a = np.linspace(0, amax, 64)
bin_edges_ecc = np.linspace(0, eccmax, 64)
bin_edges_inc = np.linspace(0, incmax, 64)

# Process Snapshots
myhist_a = np.zeros([bin_edges_a.shape[0]-1, len(nsteps)])
myhist_ecc = np.zeros([bin_edges_ecc.shape[0]-1, len(nsteps)])
myhist_inc = np.zeros([bin_edges_inc.shape[0]-1, len(nsteps)])
tout = np.zeros(len(nsteps))
print "// Computing Histograms"
for istep, nstep in enumerate(nsteps):
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    tout[istep] = snapshot.tout
    pa = []; pecc = []; pinc = []; pmass = np.zeros(snapshot.nparticles)
    for ii, particle in enumerate(snapshot.particles):
        pa.append(particle.a)
        pecc.append(particle.ecc)
        pinc.append(particle.inc * 360.0/twopi)
        pmass[ii] = particle.m
    if args.weights == 'mass':
        weights = pmass / float(np.sum(pmass))
    elif args.weights == 'npart':
        weights = np.repeat(1.0 / snapshot.nparticles, snapshot.nparticles)
    hist_a, bin_edges_a = np.histogram(pa, bin_edges_a, weights = weights)
    hist_ecc, bin_edges_ecc = np.histogram(pecc, bin_edges_ecc, weights = weights)
    hist_inc, bin_edges_inc = np.histogram(pinc, bin_edges_inc, weights = weights)
    myhist_a[:,istep] = hist_a
    myhist_ecc[:,istep] = hist_ecc
    myhist_inc[:,istep] = hist_inc

# Make Figure
print "// Generating Figures"

#
# Semi-Major Axis
#
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
im = ax.imshow(np.rot90(myhist_a), cmap='hot', \
               extent=[bin_edges_a[0], bin_edges_a[-1], tout[0]/1.0e6, tout[-1]/1.0e6])
ax.set_aspect('auto')
ax.grid(True)
ax.set_ylabel('t [Myr]', rotation='horizontal', size='small')
ax.set_xlabel('a [AU]', size='small')
simpath = os.getcwd().split("/")[-4:]
title = "/".join(simpath)
if args.weights == 'mass':
    title += ' -- (Mass Weighted)'
elif args.weights == 'npart':
    title += ' -- (Particle Number Weighted)'
ax.set_title(title, size='small')
plt.colorbar(im)
if args.weights == 'mass':
    fig.savefig('hist_a_mass.pdf')
elif args.weights == 'npart':
    fig.savefig('hist_a_npart.pdf')
plt.close()
plt.clf()

#
# Eccentricity
#
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
im = ax.imshow(np.rot90(myhist_ecc), cmap='hot', \
               extent=[bin_edges_ecc[0], bin_edges_ecc[-1], \
                       tout[0]/1.0e6, tout[-1]/1.0e6])
ax.set_aspect('auto')
ax.grid(True)
ax.set_ylabel('t [Myr]', rotation='horizontal', size='small')
ax.set_xlabel('ecc [-]', size='small')
simpath = os.getcwd().split("/")[-4:]
title = "/".join(simpath)
if args.weights == 'mass':
    title += ' -- (Mass Weighted)'
elif args.weights == 'npart':
    title += ' -- (Particle Number Weighted)'
ax.set_title(title, size='small')
plt.colorbar(im)
if args.weights == 'mass':
    fig.savefig('hist_ecc_mass.pdf')
elif args.weights == 'npart':
    fig.savefig('hist_ecc_npart.pdf')
plt.close()
plt.clf()

#
# Inclination
#
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
im = ax.imshow(np.rot90(myhist_inc), cmap='hot', \
               extent=[bin_edges_inc[0], bin_edges_inc[-1], \
                       tout[0]/1.0e6, tout[-1]/1.0e6])
ax.set_aspect('auto')
ax.grid(True)
ax.set_ylabel('t [Myr]', rotation='horizontal', size='small')
ax.set_xlabel('i [deg]', size='small')
simpath = os.getcwd().split("/")[-4:]
title = "/".join(simpath)
if args.weights == 'mass':
    title += ' -- (Mass Weighted)'
elif args.weights == 'npart':
    title += ' -- (Particle Number Weighted)'
ax.set_title(title, size='small')
plt.colorbar(im)
if args.weights == 'mass':
    fig.savefig('hist_inc_mass.pdf')
elif args.weights == 'npart':
    fig.savefig('hist_inc_npart.pdf')
plt.close()
plt.clf()
