from glob import glob
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['lines.linewidth'] = 0.5
import matplotlib.pyplot as plt

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--ellipses", action='store_true', \
                    help='Plot Orbit Ellipses for Particles.')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--test', action='store_true', \
                   help="Plot Test Set of Snapshots.")
group.add_argument('--t0', action='store_true', \
                   help="Plot First Snapshot.")
args = parser.parse_args()

# Debug Info
if args.ellipses:
    print "// Plotting XY Positions With Orbital Ellipses"
else:
    print "// Plotting XY Positions Without Orbital Ellipses"

# Full Set
if args.all:
    print "// Scanning Snapshots"
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

for istep, nstep in enumerate(nsteps):
    print "// Plotting %i/%i" % (istep+1, len(nsteps))
    npz = np.load('Snapshot_%012d.npz' % nstep)
    snapshot = npz['snapshot'][()]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for particle in snapshot.particles:
        ax.plot(particle.x, particle.y, 'b.')
        if args.ellipses:
            if not snapshot.ellipses:
                particle.compute_ellipse()
            ax.plot(particle.xell, particle.yell, 'b-')
        ax.grid(True)
    ax.set_title('t=%.2e / nstep=%010d / nparticles=%i' % \
                 (snapshot.tout, snapshot.nstep, snapshot.nparticles))
    ax.set_xlim([-6,6])
    # ax.set_ylim([-5.5,5.5])
    ax.set_aspect('equal', 'datalim')
    ax.set_xlabel('X [AU]')
    ax.set_ylabel('Y [AU]')
    plt.close()
    plt.clf()
