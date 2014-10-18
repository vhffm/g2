"""
Helper Functions For Fits.
http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html

Conventions:
- Base-10: log = np.log10()
- Base-E :  ln = np.log()
"""

import numpy as np

# ################################### #
# Polynomial Fits Using Numpy Methods #
# ################################### #
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

# ############################# #
# Self-Rolled Least-Square Fits #
# ############################# #

# Exponential Least Squares Fit
# http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
def lqexp(x,y):
    t1 = np.sum(x**2.0 * y)
    t2 = np.sum(y * np.log(y))
    t3 = np.sum(x * y)
    t4 = np.sum(x * y * np.log(y))
    t5 = np.sum(y)
    t6 = np.sum(x**2.0 * y)
    t7 = np.sum(x * y)**2.0
    a = ( t1 * t2 - t3 * t4 ) / ( t5 * t6 - t7 )
    
    t8 = np.sum(y)
    t9 = np.sum(x * y * np.log(y))
    t10 = np.sum(x*y)
    t11 = np.sum(y * np.log(y))
    t12 = np.sum(y)
    t13 = np.sum(x**2.0 * y)
    t14 = np.sum(x * y)**2.0
    b = ( t8 * t9 - t10 * t11 ) / ( t12 * t13 - t14 )
    
    return np.exp(a), b

# Power Law Least-Squares Fit
# http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
def lqpow(x,y):
    n = len(x)
    t1 = np.sum(np.log(x) * np.log(y))
    t2 = np.sum(np.log(x))
    t3 = np.sum(np.log(y))
    t4 = np.sum(np.log(x)**2.0)
    t5 = np.sum(np.log(x))**2.0
    b = ( n * t1 - t2 * t3 ) / ( n * t4 - t5 )
    
    t6 = np.sum(np.log(y))
    t7 = np.sum(np.log(x))
    a = ( t6 - b * t7 ) / n
    
    return np.exp(a), b

# Linear Least Squares Fit
def lqlin(x,y):
    a, b = np.polyfit(x, y, 1)
    return a, b
