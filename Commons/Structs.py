from kepler_helpers import cart2kep, kep2cart, compute_ellipse
import numpy as np

class Particle():
    """
    Semi-Major Axis               - a           - [AU]
    Eccentricity                  - ecc         - []
    Inclination                   - inc         - [rad]
    Longitude of Ascending Node   - Omega       - [rad]
    Argument of Perigee           - omega       - [rad]
    Orbital Phase (Mean Anomaly)  - M0          - [rad]
    XYZ-Position                  - x, y, z     - [AU]
    XYZ-Velocity                  - vx, vy, vz  - [AU/yr]
                                                - [AU/day/0.01720209895]
    Particle Mass                 - m           - [Msun]
    Particle ID                   - id          - []

    To convert velocities to km/s, do *au2km/24.0/3600.0*0.0172020989.

    NB: For G=M=1, 1 Year = 2 Pi (From Kepler's Third Law).
    """
    def __init__(self):
        # General
        self.id = None
        # Cartesian
        self.x, self.y, self.z = None, None, None
        self.vx, self.vy, self.vz = None, None, None
        # Keplerian
        self.a = None; self.ecc = None; self.inc = None
        self.Omega = None; self.omega = None; self.M0 = None
        # Other
        self.m, self.r = None, None
        # Ellipse
        self.xell = None; self.yell = None; self.zell = None

    def cart2kep(self):
        self.a, self.ecc, self.inc, self.Omega, self.omega, self.M0 = \
            cart2kep(np.array([self.x, self.y, self.z]), \
                     np.array([self.vx, self.vy, self.vz]), self.m, \
                     central_mass=1.0)
    def kep2cart(self):
        x, v = \
            kep2cart(self.a, self.ecc, self.inc, \
                     self.Omega, self.omega, self.M0, self.m, \
                     central_mass=1.0)
        self.x = x[0]; self.y = x[1]; self.z = x[2]
        self.vx = v[0]; self.vy = v[1]; self.vz = v[2]

    def compute_ellipse(self):
        self.xell, self.yell, self.zell = compute_ellipse(self.a, \
                                                          self.ecc, \
                                                          self.inc, \
                                                          self.Omega, \
                                                          self.omega)

class Planet(Particle):
    pass

class Jupiter(Planet):
    pass

class Saturn(Planet):
    pass

class Snapshot():
    def __init__(self):
        self.particles = []
        self.tout = None
        self.nstep = None
        self.nparticles = None
        self.ellipses = None
