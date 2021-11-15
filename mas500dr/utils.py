import numpy as np

def inv_median(a):
    """
    Inverse of the median of array a.
    This can be used as the `scale` argument of
    ccdproc.combine when combining flat frames.
    See CCD Data Reduction Guide Sect. 4.3.1
    """
    return 1 / np.median(a)
