import ccdproc as cdp
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
import os
from pathlib import Path
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files
from photutils import detect_sources
from astrowidgets import ImageWidget
import matplotlib.pyplot as plt

def inv_median(a):
    """
    Inverse of the median of array a.
    This can be used as the `scale` argument of
    ccdproc.combine when combining flat frames.
    See CCD Data Reduction Guide Sect. 4.3.1
    """
    return 1 / np.median(a)

def create_ccd_mask(flat1, flat2):
    ## From chapter 6.2.1. Detecting bad pixels with ccdmask
    ## It is required to have 2 flat frames with different count
    print('Creating ccd mask...')
    ratio = flat1.divide(flat2)
    maskr = cdp.ccdmask(ratio)
    ##write the mask to file
    mask_as_ccd = CCDData(data=maskr.astype('uint8'), unit=u.dimensionless_unscaled)
    mask_as_ccd.header['imagetyp'] = 'flat mask'
    mask_as_ccd.write('mask_from_ccdmask.fits')
    return maskr
