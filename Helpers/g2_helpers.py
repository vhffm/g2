import numpy as np

twopi = 2.0 * np.pi
pihalf = np.pi / 2.0
MEarth = 3.0024584e-6       # Solar Masses

def mkline(x1, y1, x2, y2):
    m = (y2 - y1) / (x2 - x1)
    n = y1 - m * x1
    return m, n
