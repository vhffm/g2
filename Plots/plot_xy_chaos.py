"""
Plot Three XY Planes. Also Plots Separation For Three Panels.

File List Format:
tag01,/path/to/xchaos01.npz
tag02,/path/to/xchaos01.npz
...
"""

import matplotlib as mpl
mpl.rcParams['backend'] = 'Agg'
import numpy as np
import matplotlib.pyplot as plt
from kepler_helpers import compute_ellipse
import argparse
import sys
from matplotlib.ticker import MaxNLocator

# Colors
dgrey = (0.5, 0.5, 0.5, 1.0)
grey  = (0.7, 0.7, 0.7, 1.0)
mred  = (1.0, 0.7, 0.7, 1.0)
mblue = (0.7, 0.7, 1.0, 1.0)
lgrey = (0.9, 0.9, 0.9, 1.0)
black = (0.0, 0.0, 0.0, 1.0)

# Lines
mpl.rcParams['lines.linewidth'] = 0.5

# Lighten labels
# http://matplotlib.org/users/customizing.html
mpl.rcParams['axes.labelcolor']  = grey
mpl.rcParams['xtick.color']      = grey
mpl.rcParams['ytick.color']      = grey
mpl.rcParams['axes.edgecolor']   = grey

# Adjust ticks
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.minor.size'] = 0
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.minor.size'] = 0

# Adjust Font Size
mpl.rcParams['xtick.labelsize']  = 'x-small'
mpl.rcParams['ytick.labelsize']  = 'x-small'
mpl.rcParams['axes.labelsize']   = 'small'
mpl.rcParams['axes.titlesize']   = 'medium'

# Parse Arguments
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--all', action='store_true', \
                   help="Plot Full Set of Snapshots.")
group.add_argument('--custom', type=int, nargs='+', \
                   help="Plot Custom Snapshot.")
args = parser.parse_args()

# List of Files and Tags
if sys.stdin.isatty():
    print "!! No File List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    filenames = []
    tags = []
    for line in lines:
        fname = line.split(",")[1]
        tag = line.split(",")[0]
        filenames.append(fname)
        tags.append(tag)
    print "// Reading %i Files" % len(filenames)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Load Data
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

tt = []; ds = []
x1 = []; x2 = []
y1 = []; y2 = []
a1 = []; a2 = []
e1 = []; e2 = []
i1 = []; i2 = []
O1 = []; O2 = []
o1 = []; o2 = []
for ifilename, filename in enumerate(filenames):
    first = True
    npz = np.load("%s" % filename)
    if first:
        ds_loc = npz["ds"]
        x1_loc = npz["x1"]
        x2_loc = npz["x2"]
        y1_loc = npz["y1"]
        y2_loc = npz["y2"]
        a1_loc = npz["a1"]
        a2_loc = npz["a2"]
        e1_loc = npz["ecc1"]
        e2_loc = npz["ecc2"]
        i1_loc = npz["inc1"]
        i2_loc = npz["inc2"]
        O1_loc = npz["Omega1"]
        O2_loc = npz["Omega2"]
        o1_loc = npz["omega1"]
        o2_loc = npz["omega2"]
        first = False
    else:
        ds_loc = np.concatenate((ds_loc, npz["ds"]), axis=1)
        x1_loc = np.concatenate((x1_loc, npz["x1"]), axis=1)
        x2_loc = np.concatenate((x2_loc, npz["x2"]), axis=1)
        y1_loc = np.concatenate((y1_loc, npz["y1"]), axis=1)
        y2_loc = np.concatenate((y2_loc, npz["y2"]), axis=1)
        a1_loc = np.concatenate((a1_loc, npz["a1"]), axis=1)
        a2_loc = np.concatenate((a2_loc, npz["a2"]), axis=1)
        e1_loc = np.concatenate((e1_loc, npz["ecc1"]), axis=1)
        e2_loc = np.concatenate((e2_loc, npz["ecc2"]), axis=1)
        i1_loc = np.concatenate((i1_loc, npz["inc1"]), axis=1)
        i2_loc = np.concatenate((i2_loc, npz["inc2"]), axis=1)
        O1_loc = np.concatenate((O1_loc, npz["Omega1"]), axis=1)
        O2_loc = np.concatenate((O2_loc, npz["Omega2"]), axis=1)
        o1_loc = np.concatenate((o1_loc, npz["omega1"]), axis=1)
        o2_loc = np.concatenate((o2_loc, npz["omega2"]), axis=1)
    ds.append(ds_loc)
    x1.append(x1_loc)
    x2.append(x2_loc)
    y1.append(y1_loc)
    y2.append(y2_loc)
    a1.append(a1_loc)
    a2.append(a2_loc)
    e1.append(e1_loc)
    e2.append(e2_loc)
    i1.append(i1_loc)
    i2.append(i2_loc)
    O1.append(O1_loc)
    O2.append(O2_loc)
    o1.append(o1_loc)
    o2.append(o2_loc)
    tt.append(npz["tt"])

# Determine Steps
if args.all:
    nsteps = range(len(ds[0]))
if args.custom:
    # Build Snapshot Number Array (Form Input)
    nsteps = range(args.custom[0], \
                   args.custom[1] + args.custom[2], \
                   args.custom[2])

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Plot Video
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# iparticle = np.array([1,6])
iparticle = 2

for istep in nsteps:
    fig = plt.figure()
    print "// Rendering Frame %06d/%06d" % (istep,nsteps[-1])
    for ii in [ 0, 1, 2 ]:
        ax = plt.subplot2grid((2,3), (0, ii))
        ax.set_aspect('equal')
        
        # Compute ellipses
        xell1, yell1, zell1 = compute_ellipse(a1[ii][istep,iparticle], e1[ii][istep,iparticle], \
                                              i1[ii][istep,iparticle], O1[ii][istep,iparticle], \
                                              o1[ii][istep,iparticle])
        xell2, yell2, zell2 = compute_ellipse(a2[ii][istep,iparticle], e2[ii][istep,iparticle], \
                                              i2[ii][istep,iparticle], O2[ii][istep,iparticle], \
                                              o2[ii][istep,iparticle])
        
        # Draw points
        ax.plot(x1[ii][istep,iparticle], y1[ii][istep,iparticle], 'bo', alpha=0.5, markeredgecolor='k', linewidth=0.5)
        ax.plot(x2[ii][istep,iparticle], y2[ii][istep,iparticle], 'ro', alpha=0.5, markeredgecolor='k', linewidth=0.5)
        
        # Draw ellipses
        ax.plot(xell1, yell1, linewidth=0.5, color='b', alpha=0.5)
        ax.plot(xell2, yell2, linewidth=0.5, color='r', alpha=0.5)

        # Fix Ticks
        ax.xaxis.set_major_locator(MaxNLocator(prune='both', nbins='7'))
        ax.yaxis.set_major_locator(MaxNLocator(prune='both', nbins='7'))

        if ii == 1 or ii == 2:
            ax.yaxis.set_ticklabels([])

        if ii == 0:
            ax.set_ylabel('Y (AU)')

        ax.set_xlabel('X (AU)')

        # ax.set_xlim([-5,5])
        # ax.set_ylim([-5,5])
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        ax.set_title("%s | ds=%.2e" % (tags[ii], ds[ii][istep,iparticle]), \
                     fontsize='x-small')
        ax.grid(False)

    ax = plt.subplot2grid((2,3), (1,0), colspan=3)
    colors = [ 'b', 'r', 'g' ]
    for ii in [ 0, 1, 2 ]:
        ax.plot(tt[ii][:istep,iparticle], ds[ii][:istep,iparticle], color=colors[ii], linewidth=1.0, alpha=0.5, label=tags[ii])
        ax.plot(tt[ii][istep,iparticle], ds[ii][istep,iparticle], 'o', color=colors[ii], alpha=0.5)
    ax.set_xlim([0,2000])
    ax.set_ylim([0,4])
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Separation (AU)')
    # Fix Ticks
    ax.yaxis.set_ticks([0,1,2,3,4])
    # ax.yaxis.set_ticklabels([0,1,2,3,4])
    # Legend
    ax.legend(loc='upper right', fontsize='x-small')

    fig.suptitle("t=%.2f yr" % tt[0][istep,0], fontsize='small')
    fig.savefig("%06d.png" % istep)
    plt.close(fig)

print "// Done"
