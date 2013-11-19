"""
Compute Core Temperature During Orbit.

Numerical Stability and Convergence Requires dt/dx**2 < 0.5.

Also see:
- http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/keplerOrbit.html
- http://www.ewp.rpi.edu/hartford/~ernesto/S2004/CHT/Notes/s06.pdf
- http://www.ewp.rpi.edu/hartford/~ernesto/S2006/CHT/Notes/ch03.pdf
- https://www.amherst.edu/media/view/103447/original/Diffusion%20part%201%20PDF.pdf
- http://www.math.ubc.ca/~feldman/m267/heatSln.pdf
- http://www.timteatro.net/2010/10/29/performance-python-solving-the-2d-diffusion-equation-with-numpy/
- http://en.wikipedia.org/wiki/Finite_difference_method#Explicit_method
"""

import numpy as np
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--alpha', type=int, default=50, \
                    help="Thermal Diffusivity [mm2/s].")
parser.add_argument('--nsteps', type=int, default=10000, \
                    help="Number of Steps.")
parser.add_argument('--fout', type=int, default=1000, \
                    help="Output Frequency.")
args = parser.parse_args()

# 
# Model Parameters
#
alpha = args.alpha    # Thermal Diffusivity [mm2/s]
nsteps = args.nsteps  # Time Steps
nr = 64               # Grid Zones

#
# Configuration
#
C = 0.5             # CFL Like Stability Factor
fout = args.fout    # Output Profile Frequency

#
# SI
#
alpha = alpha / 1000. / 1000.       # m2/s

#
# Set Grid [m]
#
r, dr = np.linspace(1000,0,nr,endpoint=False,retstep=True)
r = r[::-1]; dr *= -1

#
# Set Initial Profile [K]
# 
T0 = 273.0 + 50.0; Tb = 273.0
T = r * (Tb - T0) / r[-1] + T0

#
# Create Left/Right Ghost Cells For BCs
#
r = np.insert(r,0,0)
r = np.append(r,r[-1]+dr)
T = np.insert(T,0,T[0])
T = np.append(T,T[-1])

#
# Debug
#
print "// Initial Profile"
print T

#
# Evolve Heat Equation
# 
dt = 0.5 * dr**2. * C
print ""
print "// dt=%.2fs" % dt
foutc = fout
for nstep in range(nsteps):
    for ii in range(1,nr):
        x1 = (r[ii]+dr/2.)**2. * (T[ii+1] - T[ii])/dr
        x2 = (r[ii]-dr/2.)**2. * (T[ii] - T[ii-1])/dr
        x3 = (x1 - x2) / r[ii]**2. / dr
        T[ii] = T[ii] + dt * alpha * x3
    # Force BCs
    T[0] = T[1]
    T[-1] = Tb
    T[-2] = T[-1]
    # Output?
    if foutc == nstep:
        foutc += fout
        print ""
        print "// Evolved %i/%i Steps" % (nstep, nsteps)
        print "// Evolved %.2f Days" % (dt*nstep/3600./24.) 
        print "// Current Profile"
        print T

#
# Final Output
#
print ""
print "// Evolved %i Steps" % nsteps
print "// Evolved %.2f Days" % (dt*nsteps/3600./24.) 
print "// Final Profile"
print T
print ""
