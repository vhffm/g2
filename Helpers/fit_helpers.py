"""
Helper Functions For Fits.
http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html

Conventions:
- Base-10: log = np.log10()
- Base-E :  ln = np.log()
"""

import numpy as np

def fit_pow(x, f):
    """
    Fit Power Law.

    f(x)      = A * x**B.
    log(f(x)) = log(A) + B * log(x)
    """
    fc = np.polyfit(np.log10(x), np.log10(f), 1)
    fc = fc[-1::-1]
    A = 10.0**(fc[0])
    B = fc[1]
    return A, B

def fit_exp(x, f):
    """
    Fit Exponential.

    f(x)     = A * exp(B * x).
    ln(f(x)) = ln(A) + B * x
    """
    fc = np.polyfit(x, np.log(f), 1)
    fc = fc[-1::-1]
    A = np.exp(fc[0])
    B = fc[1]
    return A, B
