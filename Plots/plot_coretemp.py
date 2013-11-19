"""
Plot Core Temperature Profile.
"""

import numpy as np
from glob import glob
import matplotlib.pyplot as plt

# Steps
nsteps = []
globs = glob("CoreTemp_*.npz")
globs = sorted(globs)
for g in globs:
    nstep = int(g.split("_")[1].split(".")[0])
    nsteps.append(nstep)
print "// Found %i Snapshots" % len(nsteps)

# Scan Limits
print "// Scanning Limits"
first = True
for istep, nstep in enumerate(nsteps):
    npz = np.load('CoreTemp_%012d.npz' % nstep)
    T = npz["T"]
    if first:
        Tmax = np.nanmax(T); rmax = np.nanmax(npz["r"])
        Tmin = np.nanmin(T); rmin = np.nanmin(npz["r"])
        first = False
    else:
        if np.nanmax(T) > Tmax: Tmax = np.nanmax(T)
        if np.nanmin(T) < Tmin: Tmin = np.nanmin(T)

# Generate Plots
for istep, nstep in enumerate(nsteps):
    print "// Plotting %i/%i" % (istep+1, len(nsteps))
    npz = np.load('CoreTemp_%012d.npz' % nstep)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(npz["r"], npz["T"]-273.0, 'b-d', lw=1.0)
    ax.set_xlim([rmin, rmax])
    ax.set_ylim([(0.9 * Tmin) - 273.0, (1.1 * Tmax) - 273.0])
    ax.set_xlabel('r [m]')
    ax.set_ylabel('T [C]')
    ax.grid(True)
    fig.suptitle('t=%i days / nstep=%012d / r_orbit=%.2f au' % \
                (npz["tout"], npz["nstep"], npz["r_orbit"]))
    fig.savefig('T_%012d.png' % nstep)
    fig.clf()
    plt.close()
