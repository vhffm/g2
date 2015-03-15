"""
Some Initial Conditions. Also Helpers and Wrappers.
"""

import numpy as np
import constants as C
import kepler_helpers as kh
import physics_helpers as ph
import vector_helpers as vh


def Solar2():
    """
    NASA Horizon Query
    Origin = Sun, Body Center
    Epoch = 01-01-2014, 00:00 (Central Time)
    Earth = Earth (No Moon)

    @returns plist - array of IC lines for planets - [various]
    """

    mercury = "0 1 1.6515006786989092e-07 1.6310392545626536e-05 +1.1972692892595370e-01 -4.3426096069639630e-01 -4.6466929699041107e-02 +2.1482009903780699e-02 +8.9119212513564267e-03 -1.2427988748894680e-03 0 0 0"
    venus = "0 2 2.4468352521240763e-06 4.0455121182840900e-05 -4.9881994854843653e-02 +7.1760434138711926e-01 +1.2711854910906791e-02 -2.0246486048184902e-02 -1.5114700566799760e-03 +1.1477577544001651e-03 0 0 0"
    earth = "0 3 3.0023628776833749e-06 4.2587504470568303e-05 -1.7558926311807671e-01 +9.6752446474730014e-01 -2.9820874849228900e-05 -1.7207909676081830e-02 -3.1365065792022962e-03 +1.1244539739591660e-07 0 0 0"
    mars = "0 4 3.2125081695239055e-07 2.2660750299046701e-05 -1.5124387605410869e+00 +6.9681542097585980e-01 +5.1724065828569588e-02 -5.3309348873500006e-03 -1.1514537719341871e-02 -1.1040581507850680e-04 0 0 0"
    jupiter = "0 5 9.5420039213714751e-04 4.6732616936774455e-04 -1.3307968217788140e+00 +5.0187260363223007e+00 +8.9354534335869748e-03 -7.3915532677911646e-03 -1.5771896615651140e-03 +1.7200784821425221e-04 0 0 0"
    saturn = "0 6 2.8570710371524815e-04 3.8925687652332965e-04 -6.8851925827531151e+00 -7.0754782545976083e+00 +3.9713562554594872e-01 +3.6895716420435351e-03 -3.9070640598097067e-03 -7.8752021502553231e-05 0 0 0"
    uranus = "0 7 4.3642853551857632e-05 1.6953449825499186e-04 +1.9645011082335920e+01 +3.9225997140108890e+00 -2.3986827787041001e-01 -8.0447801422034405e-04 +3.6714036954299551e-03 +2.4141556323116709e-05 0 0 0"
    neptune = "0 8 5.1480569101603742e-05 1.6458790379443301e-04 +2.7063448807961450e+01 -1.2891858367498200e+01 -3.5810433631098909e-01 +1.3227724210187141e-03 +2.8500466203526960e-03 -8.9055141344433750e-05 0 0 0"
    pluto = "0 9 6.5808657181639943e-09 7.7608056333903304e-06 +6.2567672814280284e+00 -3.1924879744181990e+01 +1.6063792641261310e+00 +3.1299238024423749e-03 -3.0015467471074471e-05 -8.9263215128274780e-04 0 0 0"

    plist = [ mercury, venus, earth, mars, \
              jupiter, saturn, uranus, neptune, \
              pluto ]

    # Fix velocity units.
    # Fix format for time, particle ID.
    for iplanet, planet in enumerate(plist):
        line = planet.strip().split()
        vx = float(line[7])
        vy = float(line[8])
        vz = float(line[9])
        # Convert velocities (au/day => au/day @ G=Msun=1)
        # Length of a day in units of 1 year = 2 pi
        vx *= 365.25/C.twopi
        vy *= 365.25/C.twopi
        vz *= 365.25/C.twopi
        # t i m r x y z vx vy vz Sx Sy Sz
        # 0 1 2 3 4 5 6  7  8  9 10 11 12
        line_new = "0.0 %05d %s %s " % (int(line[1]), line[2], line[3])
        line_new += "%s %s %s " % (line[4], line[5], line[6])
        line_new += "%+.16e %+.16e %+.16e " % (vx, vy, vz)
        line_new += "0.0 0.0 0.0"
        # Reinsert
        plist[iplanet] = line_new
    
    pnames = [ "mercury", "venus", "earth", "mars", \
               "jupiter", "saturn", "uranus", "neptune", \
               "pluto" ]

    return plist, pnames


def MainFragmentReufer12(sim_name, earth):
    """
    Generates Main Fragment for Collisions in Reufer+ 2012.
    Tabulated Collisions. Calls MainFragmentIC().

    @param sim_name: Simulation Name from Reufer+ 2012 Table - [String]
    @param earth: Genga IC Line for Target (Earth) - []
    @return: Fragment Orbital Elements - []
    """

    # cC03
    if sim_name == "cC03p":
        mass = 0.112426
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.25, "+", 32.5 * C.d2r, 2.25e-6, mass, earth)
    if sim_name == "cC03m":
        mass = 0.112426
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.25, "-", 32.5 * C.d2r, 1.9e-4, mass, earth)

    # fA01
    if sim_name == "fA01p":
        mass = 0.111688
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.30, "+", 30.0 * C.d2r, 3.35e-6, mass, earth)
    if sim_name == "fA01p":
        mass = 0.111688
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.30, "-", 30.0 * C.d2r, 1.98e-4, mass, earth)

    # iA08
    if sim_name == "iA08p":
        mass = 0.111196
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.2, "+", 30.0 * C.d2r, 1.29e-6, mass, earth)
    if sim_name == "iA08m":
        mass = 0.111196
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.2, "-", 30.0 * C.d2r, 1.8e-4, mass, earth)

    # iA14
    if sim_name == "iA14p":
        mass = 0.11328
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.3, "+", 32.5 * C.d2r, 3.35e-6, mass, earth)
    if sim_name == "iA14m":
        mass = 0.11328
        a, e, i, Omega, omega, M = \
        MainFragmentIC(1.3, "-", 32.5 * C.d2r, 2.0e-4, mass, earth)

    # iA27
    if sim_name == "iA27p":
        mass = 0.111811
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.0, "+", 60.0 * C.d2r, 9.35e-7, mass, earth)
    if sim_name == "iA27m":
        mass = 0.111811
        a, e, i, Omega, omega, M = \
            MainFragmentIC(1.0, "-", 60.0 * C.d2r, 1.45e-4, mass, earth)

    # Return
    return a, e, i, Omega, omega, M, mass


def MainFragmentIC(alpha, beta_sign, theta, phase_offset, mass, earth):
    """
    Generates ICs for Main Fragment. After Earth-Collision.

    @param alpha: Impactor Velocity / Escape Velocity - []
    @param beta_sign: Impactor Faster or Slower than Earth? - [+/-]
    @param theta: Impact Angle - [Rad]
    @param phase_offset: Move impactor by this after collision - [Rad]
    @param mass: Impactor Mass - [Earth Masses]
    @param earth: Genga IC Line for Target (Earth) - []

    @return: Fragment Orbital Elements - []
    """

    # Earth Initial Conditions (Genga IC Format)
    earth = earth.strip().split()
    m0 = float(earth[2]) # Solar Masses
    x0 = np.array([float(earth[4]), float(earth[5]), float(earth[6])]) # AU
    v0 = np.array([float(earth[7]), float(earth[8]), float(earth[9])]) # AU/Day

    # Shift Orbit
    vesc = 11.20             # Surface Escape Velocity (km/s) [Today]
    vesc *= C.kms_to_genga   # Genga Units

    # Scalar Multiplier for Vector Components (can be 1+X or 1-X)
    if beta_sign == "-":
        beta = 1.0 - alpha * vesc / np.sqrt(np.sum(v0**2.0))
    elif beta_sign == "+":
        beta = 1.0 + alpha * vesc / np.sqrt(np.sum(v0**2.0))

    # Rotate, Scale
    vx1, vy1 = vh.rotate_xy(v0[0], v0[1], theta)
    v1 = np.array([vx1, vy1, v0[2]])
    v1 *= beta

    # Compute Impactor Orbit
    a1, e1, i1, Omega1, omega1, M1 = \
        kh.cart2kep(x0, v1, mass * C.mearth / C.msun)

    # Shift Phase
    M1 += phase_offset
    while M1 > C.twopi:
        M1 -= C.twopi

    # Return
    return a1, e1, i1, Omega1, omega1, M1


def Cone(x_in, v_in, r_in, nr=32, nalpha=128, aspect=0.1):
    """
    Generates Cone of Test Particles.

    @param x_in: Position Vector (XYZ) of Main Fragment - [AU]
    @param v_in: Velocity Vector (XYZ) of Main Fragment - [AU/Day]
    @param r_in: Radius of Test Particles - [AU]
    @param nr: Radial Samples of Cone - []
    @param nalpha: Aziumthal Samples of Cone - []
    @param aspect: Aspect Ratio (X/Y) of Cones - []

    @return: Genga ICs Lines for Cone Test Particles - []
    """

    # Generate Cone Ranges
    xr, dr = np.linspace(np.pi/64.0, 0.0, nr, endpoint=False, retstep=True)
    xalpha = np.linspace(-np.pi, np.pi, nalpha, endpoint=False)
    xr = xr[::-1]; dr = -dr

    # Allocate Arrays
    vold = np.tile(v_in[np.newaxis,:], (nr*nalpha,1))
    vnew = np.zeros([nr*nalpha,3])

    # Generate Velocity Vectors
    # @todo - Clearly, this code can be vectorized, etc.
    #       - It would probably take an hour of coding for neglible gain.
    #       - So, meh.
    for ir, r in enumerate(xr):
        for ialpha, alpha in enumerate(xalpha):
            # Index Juggling
            ii = ir * nalpha + ialpha
            # XY in Cone (=Theta/Phi)
            theta = r * np.cos(alpha)
            phi = r * np.sin(alpha) * aspect
            # Velocity Normalization
            vnorm = np.sin(np.pi/2.0 - (r / (np.pi/64.0 + dr)))
            # Rotations
            Rz = np.array([[ np.cos(theta), -np.sin(theta), 0.0 ], \
                           [ np.sin(theta),  np.cos(theta), 0.0 ], \
                           [ 0.0, 0.0, 1.0 ]])
            Ry = np.array([[ np.cos(phi), 0.0, np.sin(phi) ], \
                           [ 0.0, 1.0, 0.0 ], \
                           [ - np.sin(phi), 0.0, np.cos(phi) ] ])
            # Rotate, Scale, Replace
            vnew[ii,:] = np.dot(Ry, np.dot(Rz, vold[ii,:])) * vnorm

    # Generate IC Lines
    lines_out = []
    for ii in range(vnew.shape[0]):
        line_new = "0.0 %05d %.16e %.16e " % (ii + 10001, 0.0, r_in)
        line_new += "%+.16e %+.16e %+.16e " % (x_in[0], x_in[1], x_in[2])
        line_new += "%+.16e %+.16e %+.16e " % (vnew[ii,0], \
                                               vnew[ii,1], \
                                               vnew[ii,2])
        line_new += "0.0 0.0 0.0"
        lines_out.append(line_new)

    # Return
    return lines_out
