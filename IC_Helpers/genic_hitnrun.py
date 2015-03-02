"""
Generate Hit'n'Run ICs. 

a) Calculate impactor orbit corresponding to impact angle, velocity, mass.
b) Iteratively move impactor orbital phase to avoid Earth collision.

For step (b), we rely on Genga to integrate to orbit for a bit.
Make sure Genga is compiled to ignore the lock.dat.
"""

import subprocess as sp
import kepler_helpers as kh
import constants as C
import numpy as np
import pandas as pd

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
# HELPER FUNCTION PARTY
#
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def rotate_xy(xin, yin, theta):
    xout = np.cos(theta) * xin - np.sin(theta) * yin
    yout = np.sin(theta) * xin + np.cos(theta) * yin
    return xout, yout

# mass in kg; default density is 2 g/cc
# earth is 5.5
# returns r in km
def mass2radius(mass, rho=2.0):
    rho /= 1000.0 # kg/cc
    rho *= (100.0*1000.0)**3.0 # kg/km3
    V = mass/rho # km3
    r = (V/np.pi * 0.75)**(1./3.) # km
    return r

def param_file():
    lines = "Time step in days = 6.0\n"
    lines += "Output name = _hitnrun_test\n"
    lines += "Energy output intervall = 100\n"
    lines += "Coordinates output intervall = 100\n"
    lines += "Number of outputs per intervall = 1\n"
    lines += "Integration steps = 100\n"
    lines += "Central Mass = 1.0\n"
    lines += "n1 = 3.0\n"
    lines += "n2 = 0.1\n"
    lines += "Input file = initial.dat\n"
    lines += "Input file Format: << t i m r x y z vx vy vz Sx Sy Sz >>\n"
    lines += "Default rho = 2.0\n"
    lines += "Use Test Particles = 0\n"
    lines += "Restart timestep = 0\n"
    lines += "Minimum number of bodies = 1\n"
    lines += "Order of integrator = 2\n"
    lines += "aeGrid amin = 0\n"
    lines += "aeGrid amax = 20\n"
    lines += "aeGrid emin = 0\n"
    lines += "aeGrid emax = 1\n"
    lines += "aeGrid Na = 10\n"
    lines += "aeGrid Ne = 10\n"
    lines += "aeGrid Start Count = 10000000000000\n"
    lines += "aeGrid name = A\n"
    lines += "Gas dTau_diss = 3000000\n"
    lines += "Gas alpha = 1\n"
    return lines

def ic_file(a1, e1, i1, Omega1, omega1, M1, impactor_mass):
    # Earth Initial Conditions
    m0 = 3.0023628776833749e-06 # solar masses
    x0 = np.array([ -1.7558926311807671e-01, \
                    +9.6752446474730014e-01, \
                    -2.9820874849228900e-05 ]) # AU
    v0 = np.array([ -1.0003188990792635e+00, \
                    -1.8232933966543841e-01, \
                    +6.5366019607806956e-06 ]) # AU/day

    # Earth IC Line
    line0 = "0 2 3.0023628776833749e-06 4.2587504470568303e-05 "
    line0 += "%.16f %.16f %.16f " % (x0[0], x0[1], x0[2])
    line0 += "%.16f %.16f %.16f " % (v0[0], v0[1], v0[2])
    line0 += "0.0 0.0 0.0\n"

    # Compute Cartesian Outputs
    x1, v1 = kh.kep2cart(a1, e1, i1, \
                         Omega1, omega1, M1, \
                         impactor_mass * C.mearth / C.msun)

    # Impactor IC Line
    line1 = "0 10 "
    line1 += "%.16f " % (impactor_mass * C.mearth / C.msun)
    line1 += "%.16f " % (mass2radius(impactor_mass * C.mearth) / C.au2km)
    line1 += "%.16f %.16f %.16f " % (x1[0], x1[1], x1[2])
    line1 += "%.16f %.16f %.16f " % (v1[0], v1[1], v1[2])
    line1 += "0.0 0.0 0.0\n"

    # Cat & Return
    lines = line0 + line1
    return lines

def impactor_ic(alpha, beta_sign, theta, phase_offset, impactor_mass, \
                display_earth):
    # Earth Initial Conditions
    m0 = 3.0023628776833749e-06 # solar masses
    x0 = np.array([ -1.7558926311807671e-01, \
                    +9.6752446474730014e-01, \
                    -2.9820874849228900e-05 ]) # AU
    v0 = np.array([ -1.0003188990792635e+00, \
                    -1.8232933966543841e-01, \
                    +6.5366019607806956e-06 ]) # AU/day

    # Base Orbit, Earth
    if display_earth:
        a0, e0, i0, Omega0, omega0, M0 = kh.cart2kep(x0, v0, m0)
        print "E(a,e,i,Omega,omega,M) = (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)" % (a0, e0, i0, Omega0, omega0, M0)

    # Shift Orbit
    vesc = 11.20   # surface escape velocity (km/s) [present day earth]
    vesc *= C.kms_to_genga   # genga units
    # theta = 35.0 * C.d2r
    # alpha = 1.2 # vrel/vesc
    # scalar multiplier for vector components (can be 1+X or 1-X)
    if beta_sign == "+":
        beta = 1.0 - alpha * vesc / np.sqrt(np.sum(v0**2.0))
    elif beta_sign == "-":
        beta = 1.0 + alpha * vesc / np.sqrt(np.sum(v0**2.0))

    # Rotate, Scale
    vx1, vy1 = rotate_xy(v0[0], v0[1], theta)
    v1 = np.array([vx1, vy1, v0[2]])
    v1 *= beta

    # Compare
    # print np.sqrt(np.sum(v**2.0))*C.genga_to_kms
    # print np.sqrt(np.sum(vnew**2.0))*C.genga_to_kms

    # Compute Impactor Orbit
    a1, e1, i1, Omega1, omega1, M1 = \
        kh.cart2kep(x0, v1, impactor_mass * C.mearth / C.msun)

    # Shift Phase
    M1 += phase_offset
    while M1 > C.twopi:
        M1 -= C.twopi

    # Compute Shifted Impactor Position
    x1, v1 = kh.kep2cart(a1, e1, i1, \
                         Omega1, omega1, M1, \
                         impactor_mass * C.mearth / C.msun)

    # Debug
    dx2 = (x1[0]-x0[0])**2.0
    dy2 = (x1[1]-x0[1])**2.0
    dz2 = (x1[2]-x0[2])**2.0
    ds = np.sqrt(dx2 + dy2 + dz2) * C.au2km
    print "I(a,e,i,Omega,omega,M) = (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)" % (a1, e1, i1, Omega1, omega1, M1)
    print "D(Earth-Impactor) = %.2f km" % (ds,)

    # Return
    return a1, e1, i1, Omega1, omega1, M1

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#
# MAIN CODE
#
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Loop ICs Until Converged
collision = True
niteration = 0
phase_offset = 6500.0 / C.au2km / C.twopi
print ""
while collision or niteration < 10:
    # Book Keeping
    niteration += 1
    print "* Iteration %03d / Phase Offset %0.2e [Rad]" % \
        (niteration, phase_offset)

    # Write Parameters
    lines = param_file()
    with open('param.dat', 'w') as f:
        for line in lines:
            f.write(line)

    # Generate ICs
    # lines = ic_file()
    alpha = 1.2
    beta_sign = "-"
    # beta_sign = "+"
    theta = 35.0
    impactor_mass = 0.009
    if niteration == 1:
        display_earth = True
    else:
        display_earth = False

    # Generate Impactor ICs
    a1, e1, i1, Omega1, omega1, M1 = \
        impactor_ic(alpha, beta_sign, theta, phase_offset, impactor_mass, \
                    display_earth)

    # Write ICs
    lines = ic_file(a1, e1, i1, Omega1, omega1, M1, impactor_mass)
    with open("initial.dat", "w") as f:
        for line in lines:
            f.write(line)

    # Call Genga
    # sp.call(["rm", "lock.dat"], stdout=open("/dev/null", "w"), close_fds=True)
    sp.call(["genga"], stdout=open("/dev/null", "w"), close_fds=True)

    # Check Final Particle Number
    with open("Collisions_hitnrun_test.dat", "r") as f:
        lines = f.readlines()
        # print len(lines)
        if len(lines) == 0:
            collision = False
        else:
            collision = True

    # If Collision, Double Phase Offset
    # Else, Half Phase Offset
    if collision:
        phase_offset *= 1.05
    else:
        phase_offset *= 0.1

    # Space
    print ""

# Load Final Orbit Configuration
# Define CSV File
names_cols = [ "time", "pid", "mass", "radius", \
               "x", "y", "z", \
               "vx", "vy", "vz", \
               "Sx", "Sy", "Sz", \
               "amin", "amax", "emin", "emax", \
               "aecount", "aecountT", "enccount", "test", "X" ]
touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9, 19 ]
types_cols = { "pid": np.int32 }

# Load CSV
df = pd.read_csv("Out_hitnrun_test.dat", \
                 sep=" ", lineterminator="\n", header=None, \
                 names=names_cols, dtype=types_cols, usecols=touse_cols)

# Convert, Sort
df.mass *= C.msun/C.mearth
df.sort(columns=["time"], ascending=True, inplace=True)

# Compute Kepler Elements
x = np.asarray(df[df.pid==10].tail(1).x)[0]
y = np.asarray(df[df.pid==10].tail(1).y)[0]
z = np.asarray(df[df.pid==10].tail(1).z)[0]
vx = np.asarray(df[df.pid==10].tail(1).vx)[0]
vy = np.asarray(df[df.pid==10].tail(1).vy)[0]
vz = np.asarray(df[df.pid==10].tail(1).vz)[0]
m = np.asarray(df[df.pid==10].tail(1).mass)[0]
a1, e1, i1, Omega1, omega1, M1 = kh.cart2kep([x,y,z], [vx,vy,vz], m)

# Final Impactor Orbit
print "* Gravitational-Focusing Cleaned Impactor Orbit"
print "I(a,e,i,Omega,omega,M) = (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)" % \
    (a1, e1, i1, Omega1, omega1, M1)
print ""
