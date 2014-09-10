"""
Constants and Conversions.
"""

import numpy as np

# Masses
mjupiter = 1.89e27 # kg
mearth   = 5.97e24 # kg
msun     = 1.99e30 # kg

# Angle Conversions
r2d = 180.0 / np.pi
d2r = np.pi / 180.0

# Distance Conversions
au2km = 149597871.0 # km

# True Constants
G = 6.67384e-11      # m3/kg/s2
G = G / 1000.0**3.0  # km3/kg/s2

# Convenience
twopi = 2.0 * np.pi
