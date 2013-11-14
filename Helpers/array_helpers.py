"""
Helper Functions for Array Processing.
"""

import numpy as np
import sys

def intersect2d(A, B):
    """
    Crazy 2D Intersection.

    @note
    May need to transpose arrays, i.e.
    A = np.tranpose(a).copy()
    B = np.tranpose(a).copy()

    @params
    A = np.array([[1,4],[2,5],[3,6]])
    B = np.array([[1,4],[3,6],[7,8]])

    @returns
    C = np.array([[1,4],[3,6]])

    Source:
    http://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays

    Relevant:
    http://i.imgur.com/xVyoSl.jpg
    """

    nrows, ncols = A.shape
    dtype={'names':['f{}'.format(i) for i in range(ncols)],
           'formats':ncols * [A.dtype]}

    C = np.intersect1d(A.view(dtype), B.view(dtype))

    # This last bit is optional if you're okay with "C" being a structured array...
    C = C.view(A.dtype).reshape(-1, ncols)

    return C
