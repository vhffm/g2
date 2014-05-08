"""
Plot Histograms:

Log Radial | Vertical | Angular Separation
Lin Radial | Vertical | Angular Separation
"""

import matplotlib as mpl; mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Adjust Font Size
mpl.rcParams['xtick.labelsize']  = 'x-small'
mpl.rcParams['ytick.labelsize']  = 'x-small'
mpl.rcParams['axes.labelsize']   = 'small'
mpl.rcParams['axes.titlesize']   = 'small'
mpl.rcParams['legend.fontsize']  = 'x-small'

# Constants
au = 149597871.0 # km
r2d = 360.0 / np.pi / 2.0

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--infile', default='XChaos.npz', \
                    help="Input File Name.")
args = parser.parse_args()

# Load Data
print "// Loading Data"
npz = np.load("%s" % (args.infile,))
x1 = npz["x1"]; x2 = npz["x2"]
y1 = npz["y1"]; y2 = npz["y2"]
z1 = npz["z1"]; z2 = npz["z2"]
# ds = npz["ds"]; tt = npz["tt"][:,0]
ds = npz["ds"]; tt = npz["tt"]

# Status
print "// Found %i Particles, %i Steps" % (ds.shape[1], ds.shape[0])
print "// Computing Separations"

# Relative Heliocentric Distance
r1 = np.sqrt(x1**2. + y1**2. + z1**2.)
r2 = np.sqrt(x2**2. + y2**2. + z2**2.)
dr = r1 - r2

# Cylindrical Radius
rc1 = np.sqrt(x1**2. + y1**2)
rc2 = np.sqrt(x2**2. + y2**2)
drc = rc1 - rc2

# Z Differences
dz = z1 - z2

# 2D Angle
dot_rc1_rc2 = x1 * x2 + y1 * y2
cos_theta = dot_rc1_rc2 / ( rc1 * rc2 )
dt = np.arccos(cos_theta)
dS = dt * ( rc1 + rc2 ) / 2.0

# Set Up Histograms
bins_dr = np.logspace(-12,1,128)
bins_dz = np.logspace(-12,1,128)
bins_dt = np.logspace(-9,1,128)
bins_dr = np.concatenate((-bins_dr[::-1], bins_dr), axis=0)
bins_dz = np.concatenate((-bins_dz[::-1], bins_dz), axis=0)

bins_dr_lin = np.linspace(-0.2,0.2,128)
bins_dz_lin = np.linspace(-0.2,0.2,128)
bins_dt_lin = np.linspace(0,np.pi,128)

# Plot Histograms
for istep in range(ds.shape[0]):
    print "// Rendering Output %05d" % istep
    fig = plt.figure()
    fig.set_size_inches(12,6)

    ax1 = fig.add_subplot(2,3,1)
    ax1.hist(dr[istep,:], bins_dr, \
             weights=np.ones_like(dr[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax1.set_xscale('symlog', linthreshx=1.0e-9)
    ax1.set_title("Radial")
    ax1.set_xticks([-1,-1e-3,-1e-6,-1e-9,1e-9,1e-6,1e-3,1])
    ax1.set_ylabel('Number Fraction')
    ax1.set_xlim([-10,10])
    ax1.set_ylim([0,0.05])

    ax2 = fig.add_subplot(2,3,2)
    ax2.hist(dz[istep,:], bins_dz, \
             weights=np.ones_like(dz[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax2.set_xscale('symlog', linthreshx=1.0e-9)
    ax2.set_title("Vertical")
    ax2.set_xticks([-1,-1e-3,-1e-6,-1e-9,1e-9,1e-6,1e-3,1])
    ax2.set_xlim([-10,10])
    ax2.set_ylim([0,0.05])

    ax3 = fig.add_subplot(2,3,3)
    ax3.hist(dt[istep,:]*r2d, bins_dt*r2d, \
             weights=np.ones_like(dt[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax3.set_xscale('log')
    ax3.set_title("Angle")
    ax3.set_xticks([1e-9,1e-6,1e-3,1])
    #ax3.set_xlim([1.0e-9, 500])
    ax3.set_xlim([1.0e-7, 500])
    ax3.set_ylim([0,0.12])
    
    ax4 = fig.add_subplot(2,3,4)
    ax4.hist(dr[istep,:], bins_dr_lin, \
             weights=np.ones_like(dr[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax4.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
    ax4.set_xlabel('Separation (AU)')
    ax4.set_ylabel('Number Fraction')
    ax4.set_xlim([-0.2,0.2])
    ax4.set_ylim([0,0.10])
    
    ax5 = fig.add_subplot(2,3,5)
    ax5.hist(dz[istep,:], bins_dz_lin, \
             weights=np.ones_like(dz[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax5.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
    ax5.set_xlabel('Separation (AU)')
    ax5.set_xlim([-0.2,0.2])
    ax5.set_ylim([0,0.10])
    
    ax6 = fig.add_subplot(2,3,6)
    ax6.hist(dt[istep,:]*r2d, bins_dt_lin*r2d, \
             weights=np.ones_like(dt[istep,:])/2000.0, \
             facecolor='b', alpha=0.5, edgecolor='none')
    ax6.set_xticks([0, 45, 90, 135, 180])
    ax6.set_xlabel('Separation (Degree)')
    ax6.set_xlim([0, np.pi*r2d])
    ax6.set_ylim([0,0.12])

    fig.suptitle("%.2f yr" % tt[istep])
    fig.savefig('LogHist_%05d.png' % istep)
    plt.close(fig)

print "// Done"
