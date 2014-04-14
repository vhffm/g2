"""
Plot Drift Histograms (Radial, Arc).
"""

import sys
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Input File Name.")
args = parser.parse_args()

# Load Data
print "// Loading Data"
npz = np.load("%s" % args.infile)
ds = npz["ds"]
tt = npz["tt"]
x1 = npz["x1"]
x2 = npz["x2"]
y1 = npz["y1"]
y2 = npz["y2"]
z1 = npz["z1"]
z2 = npz["z2"]

# Compute Relative Heliocentric Distance
r1 = np.sqrt(x1**2. + y1**2. + z1**2.)
r2 = np.sqrt(x2**2. + y2**2. + z2**2.)
dr = r1 - r2

# Compute Angular Differences, Convert To Arc Length
dot_r1_r2 = x1 * x2 + y1 * y2 + z1 * z2
cos_theta = dot_r1_r2 / ( r1 * r2 )
theta = np.arccos(cos_theta)
dS = theta * (r1 + r2) / 2.0

# Define Histograms
bins_dr = np.linspace(-0.15,0.15,64)
bins_dS = np.linspace(0,14,64)

# Step?
# istep = 1000

# Loop Steps:
# rng = range(len(tt))
rng = range(0,10000,100)
for istep in rng:

    # Debug
    print "// Plotting %05d/%05d" % (istep, rng[-1])

    # Set up figure
    fig, axarr = plt.subplots(1,2)
    fig.set_size_inches(16.0, 6.0)

    # Left Panel, Radial Difference Histogram
    ax = axarr[0]
    ax.hist(dr[istep], bins_dr, \
            weights=np.ones_like(dr[istep])/2000.0, \
            facecolor='b', alpha=0.8, edgecolor='none')
    ax.set_xlabel('Difference (AU)')
    ax.set_ylabel('Fraction of Particles')
    ax.set_title('Radius | t=%0.2f Years' % tt[istep])
    ax.set_xlim([-0.15,0.15])

    # Right Panel, Radial Difference Histogram
    ax = axarr[1]
    ax.hist(dS[istep], bins_dS, \
            weights=np.ones_like(dr[istep,:])/2000.0, \
            facecolor='g', alpha=0.8, edgecolor='none')
    ax.set_xlabel('Difference (AU)')
    # ax.set_ylabel('Fraction of Particles')
    ax.set_title('Arc Length | t=%0.2f Years' % tt[istep])
    ax.set_xlim([0,15])

    # Save Figure
    fig.savefig("Test_%05d.png" % istep)

    # Clean Figure
    plt.close(fig)
