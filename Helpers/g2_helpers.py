import numpy as np
import scipy.interpolate as spi

twopi = 2.0 * np.pi
pihalf = np.pi / 2.0
MEarth = 3.0024584e-6       # Solar Masses
d2r = twopi / 360.0
r2d = 360. / twopi

def mkline(x1, y1, x2, y2):
    m = (y2 - y1) / (x2 - x1)
    n = y1 - m * x1
    return m, n

def interpolate_to_steps(x, y, xnew):
    """
    Returns Y Interpolated Onto XNew. For X<XNew, Y=0. For X>XNew, Y=Y[-1].
    """
    # Compute, Apply Interpolation Function. Leave Extrapolated Values NaN.
    f = spi.interp1d(x, y, kind='zero', bounds_error=False)
    ynew = f(xnew)
    # We Know Min/Max. Set Manually.
    if np.isnan(ynew[0]) and np.isnan(ynew[-1]):
        nan_left_right = np.where(np.diff(np.isnan(ynew)-1)!=0)[0]
        nan_left, nan_right = nan_left_right
        ynew[:nan_left+1] = 0
        ynew[nan_right+1:] = y[-1]
    # Left Only
    if np.isnan(ynew[0]) and not np.isnan(ynew[-1]):
        nan_left = np.where(np.diff(np.isnan(ynew)-1)!=0)[0][0]
        ynew[:nan_left+1] = 0
    # Right Only
    if not np.isnan(ynew[0]) and np.isnan(ynew[-1]):
        nan_right = np.where(np.diff(np.isnan(ynew)-1)!=0)[0][0]
        ynew[nan_right+1:] = y[-1]
    # Return
    return ynew
