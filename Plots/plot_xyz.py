"""
3D Orbital Plots.
"""

from glob import glob
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import argparse
from g2_helpers import MEarth

# Paint The Rims Black
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['lines.color'] = 'white'
mpl.rcParams['patch.edgecolor'] = 'white'
mpl.rcParams['text.color'] = 'white'
mpl.rcParams['axes.facecolor'] = 'black'
mpl.rcParams['axes.edgecolor'] = 'white'
mpl.rcParams['axes.labelcolor'] = 'white'
# mpl.rcParams['xtick.color'] = 'white'
# mpl.rcParams['ytick.color'] = 'white'
# mpl.rcParams['grid.color'] = 'white'
mpl.rcParams['figure.facecolor'] = 'black'
mpl.rcParams['figure.edgecolor'] = 'black'
# mpl.rcParams['savefig.facecolor'] = 'black'
# mpl.rcParams['savefig.edgecolor'] = 'black'
mpl.rcParams['axes3d.grid'] = True
# mpl.rcParams['grid.alpha'] = 0.5
mpl.rcParams['grid.color'] = 'black'

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--log', action='store_true', \
                    help="Scale Colors as Log10 .")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
args = parser.parse_args()

# Full Set
if args.all:
    nsteps = []
    globs = glob("Snapshot_*.npz")
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split("_")[1].split(".")[0])
        nsteps.append(nstep)

# Compute Ellipses
for istep, nstep in enumerate(nsteps):
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    if not snapshot.ellipses:
        print "// Computing Ellipse %i/%i" % (istep+1, len(nsteps))
        for particle in snapshot.particles:
            particle.compute_ellipse()
        snapshot.ellipses = True
        np.savez('Snapshot_%012d.npz' % nstep, snapshot=snapshot)

# MinMax
print "// Scanning Limits"
first = True
for istep, nstep in enumerate(nsteps):
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    for particle in snapshot.particles:
        if first:
            mmin = particle.m
            mmax = particle.m
            first = False
        else:
            if particle.m > mmax: mmax = particle.m
            if particle.m < mmin: mmin = particle.m

# Set Colormap
# cmap = cm.get_cmap('Blues_r')
# cmap = cm.get_cmap('OrRd_r')
# cmap = cm.get_cmap('PuRd_r')
# cmap = cm.get_cmap('Purples_r')
cmap = cm.get_cmap('RdPu_r')

# Draw Snapshots
for istep, nstep in enumerate(nsteps):
    print "// Processing Snapshot %i/%i" % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]

    # Set Up Plot
    fig = plt.figure(frameon=False)
    fig.set_size_inches(8,4.5)
    ax = Axes3D(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.w_xaxis.set_pane_color([0.0, 0.0, 0.0, 0.0])
    ax.w_yaxis.set_pane_color([0.0, 0.0, 0.0, 0.0])
    ax.w_zaxis.set_pane_color([0.0, 0.0, 0.0, 0.0])
    ax.plot([-7,7],[0,0],[0,0], c=(0.5,0.5,0.5,0.4), lw=0.2)
    ax.plot([0,0],[-7,7],[0,0], c=(0.5,0.5,0.5,0.4), lw=0.2)
    ax.plot([0,0],[0,0],[-6,6], c=(0.5,0.5,0.5,0.4), lw=0.2)
    ax.set_xlim3d([-4,4])
    ax.set_ylim3d([-4,4])
    ax.set_zlim3d([-4,4])
    
    # Loop Particles
    mdisk = 0.0
    mmax_loc = 0.0
    for particle in snapshot.particles:
        mdisk += particle.m
        if particle.m > mmax_loc: mmax_loc = particle.m
        if args.log:
            pm = (np.log10(particle.m) - np.log10(mmin)) / \
                 (np.log10(mmax) - np.log10(mmin))
        else:
            pm = (particle.m - mmin) / (mmax - mmin)
        # Color Floor (If Cmap is Black at 0)
        # cm = cmap(np.max([0.025,pm]))
        cm = cmap(pm)
        ax.plot(particle.xell, particle.yell, particle.zell, \
                c=(cm[0], cm[1], cm[2], 0.7), \
                lw=0.5)
    # ax.elev += 10
    ax.azim -= 20
    txtime = "t=%.2e yr" % snapshot.tout
    txnstp = "nstep=%012d" % snapshot.nstep
    txnprt = "npart=%04d" % snapshot.nparticles
    txmdsk = "mdisk=%.2e [M_Earth]" % (mdisk / MEarth)
    txmmax = "mmax=%.2e [M_Earth]" % (mmax_loc / MEarth)
    lftext = txtime + " / " + txnstp + " / " + txnprt + " / " + txmdsk + " / " + txmmax
    ax.text2D(0.02, 0.95, lftext, transform=ax.transAxes, \
              color=(0.85,0.5,0.85,0.5), size='x-small')
    ax.grid(False)
    fig.savefig('OrbitsXYZ_%012d.png' % nstep, dpi=160)
