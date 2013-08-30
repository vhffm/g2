"""
Various Helper Function for Kepler Orbits.

They usually break down for inc=0 and/or ecc=0.
Might need to upgrade to Regular Equations some day
cf. http://server.faia.upm.es/moda/curso1112/kepler.pdf
"""

import numpy as np

def cart2kep(r, v, mass, central_mass=1.0):
    """
    @params
    r - (x,y,z) Cartesian Positions
    v - (vx,vy,vz) Cartesian Velocities
    mass - Particle Mass
    central_mass - Mass of Central Object

    @returns
    a - Semi-Major Axis
    ecc - Eccentricity
    inc - Inclination
    Omega - Longitude of the Ascending Node
    omega - Argument of Periapsis
    M - Mean Anomaly at Epoch

    Cf. http://www.bruce-shapiro.com/pair/ElementConversionRecipes.pdf
    """

    # Gravitational Parameter
    G = 1.0
    mu = G * ( central_mass + mass )

    # Cartesian Unit Vectors
    iHat = np.array([1., 0., 0.])
    jHat = np.array([0., 1., 0.])
    kHat = np.array([0., 0., 1.])

    # Eccentricity
    h = np.cross(r, v)
    evec = 1. / mu * np.cross(v, h) - r / np.linalg.norm(r)
    ecc = np.linalg.norm(evec)

    # Semi Major Axis
    a = np.dot(h,h) / ( mu * ( 1. - ecc**2. ))

    # Inclination
    inc = np.arccos(np.dot(kHat, h) / np.linalg.norm(h))

    # Longitude of the Ascending Node
    n = np.cross(kHat, h)
    if inc == 0.0:
        Omega = 0.0
    else:
        Omega = np.arccos(np.dot(iHat, n) / np.linalg.norm(n))
        if np.dot(n, jHat) < 0:
            Omega = 2. * np.pi - Omega

    # Argument of Perigee
    # For Zero Inclination, Fall Back to 2D Case
    # http://en.wikipedia.org/wiki/Argument_of_periapsis
    if inc == 0.0:
        omega = np.arctan2(evec[1]/ecc,evec[0]/ecc)
    else:
        omega = np.arccos(np.dot(n, evec) / (np.linalg.norm(n) * ecc))
        if np.dot(evec, kHat) < 0:
            omega = 2. * np.pi - omega

    # True Anomaly
    theta = np.arccos(np.dot(evec, r) / (ecc * np.linalg.norm(r)))
    if np.dot(r,v) < 0:
        theta = 2. * np.pi - theta

    # Eccentric Anomaly
    E = np.arccos((ecc + np.cos(theta)) / (1 + ecc * np.cos(theta)))
    if np.pi < theta and theta < 2. * np.pi:
        E = 2. * np.pi - E

    # Mean Anomaly
    M = E - ecc * np.sin(E)

    # Return Set
    return a, ecc, inc, Omega, omega, M

def kep2cart(a, ecc, inc, Omega, omega, M, mass, central_mass=1.0):
    """
    @params
    a - Semi-Major Axis
    ecc - Eccentricity
    inc - Inclination
    Omega - Longitude of the Ascending Node
    omega - Argument of Periapsis
    M - Mean Anomaly at Epoch

    @returns
    r - (x,y,z) Cartesian Positions
    v - (vx,vy,vz) Cartesian Velocities

    Cf. http://www.bruce-shapiro.com/pair/ElementConversionRecipes.pdf
    """

    # Zero Inclination has no Argument of Perigee
    # Set =0 by Convention
    if inc == 0.0:
        Omega = 0.0

    # Mean Anomaly -> Eccentric Anomaly
    E = nr(M, ecc)

    # PQW Unit Vectors
    Px = np.cos(omega) * np.cos(Omega) - \
         np.sin(omega) * np.cos(inc) * np.sin(Omega)
    Py = np.cos(omega) * np.sin(Omega) + \
         np.sin(omega) * np.cos(inc) * np.cos(Omega)
    Pz = np.sin(omega) * np.sin(inc)

    Qx = - np.sin(omega) * np.cos(Omega) - \
           np.cos(omega) * np.cos(inc) * np.sin(Omega)
    Qy = - np.sin(omega) * np.sin(Omega) + \
           np.cos(omega) * np.cos(inc) * np.cos(Omega)
    Qz =   np.sin(inc) * np.cos(omega)

    P = np.array([Px, Py, Pz])
    Q = np.array([Qx, Qy, Qz])

    # Position
    x = a * (np.cos(E) - ecc) * P + \
        a * np.sqrt(1 - ecc**2.) * np.sin(E) * Q

    # Velocity
    G = 1.0
    mu = G * ( central_mass + mass )
    Edot = np.sqrt(mu / a**3.) / ( 1. - ecc * np.cos(E) )
    v = - a * np.sin(E) * Edot * P + \
          a * np.sqrt(1. - ecc**2.) * np.cos(E) * Edot * Q

    # Return
    return x, v

def nr(M, ecc, epsilon_target=1.0e-5):
    """
    Newton-Raphson Iteration to Compute Eccentric Anomaly from Mean Anomaly.
    """

    Ei = M; ii = 1
    # print "Running Newton-Raphson for M=%.2e, e=%.2f" % ( M, ecc )
    while True:
        Ei1 = Ei - ( Ei - ecc * np.sin(Ei) - M ) / (  1 - ecc * np.cos(Ei)  )
        epsilon = np.abs(Ei1 - Ei)
        Ei = Ei1
        # print "Iteration %i, Residual %.2e" % ( ii, epsilon )
        if epsilon < epsilon_target:
            break
        ii += 1
    # print "Found E=%.2f" % Ei1
    return Ei1

def compute_ellipse(a, ecc, inc, Omega, omega):
    """
    Compute XYZ Sequence for a Kepler Ellipse.
    """

    # Eccentric Anomaly
    E = np.linspace(0.0, 2.*np.pi, 512)

    # Allocate
    x = np.zeros_like(E)
    y = np.zeros_like(E)
    z = np.zeros_like(E)

    # PQW Unit Vectors
    Px = np.cos(omega) * np.cos(Omega) - \
         np.sin(omega) * np.cos(inc) * np.sin(Omega)
    Py = np.cos(omega) * np.sin(Omega) + \
         np.sin(omega) * np.cos(inc) * np.cos(Omega)
    Pz = np.sin(omega) * np.sin(inc)

    Qx = - np.sin(omega) * np.cos(Omega) - \
           np.cos(omega) * np.cos(inc) * np.sin(Omega)
    Qy = - np.sin(omega) * np.sin(Omega) + \
           np.cos(omega) * np.cos(inc) * np.cos(Omega)
    Qz =   np.sin(inc) * np.cos(omega)

    P = np.array([Px, Py, Pz])
    Q = np.array([Qx, Qy, Qz])

    # Loop Ellipse
    for iE in range(E.shape[0]):
        r = a * (np.cos(E[iE]) - ecc) * P + \
            a * np.sqrt(1 - ecc**2.) * np.sin(E[iE]) * Q
        x[iE] = r[0]
        y[iE] = r[1]
        z[iE] = r[2]

    # Return
    return x, y, z
