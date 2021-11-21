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
import matplotlib.pyplot as plt
from  astroscrappy  import  detect_cosmics

def inv_median(a):
    """
    Inverse of the median of array a.
    This can be used as the `scale` argument of
    ccdproc.combine when combining flat frames.
    See CCD Data Reduction Guide Sect. 4.3.1
    """
    return 1 / np.median(a)

def create_ccd_mask(flat1, flat2):
    """
    From chapter 6.2.1. Detecting bad pixels with ccdmask
    It is required to have 2 flat frames with different count
    """
    print('Creating ccd mask...')
    ratio = flat1.divide(flat2)
    maskr = cdp.ccdmask(ratio)
    ##write the mask to file
    mask_as_ccd = CCDData(data=maskr.astype('uint8'), unit=u.dimensionless_unscaled)
    mask_as_ccd.header['imagetyp'] = 'flat mask'
    mask_as_ccd.write('mask_from_ccdmask.fits')
    return maskr

def cosmic_ray_correction(clean_image):
    """
    This function corrects the science calibrated image from cosmic ray,
    it is a function from astroscrappy, available at https://github.com/astropy/astroscrappy
    the explanation from all the arguments can be found at: 
    https://astroscrappy.readthedocs.io/en/latest/api/astroscrappy.detect_cosmics.html
 
    """
    print('Cleaning from cosmic ray...')
    detect_cosmics(clean_image,inmask=None, sigclip=4.0, sigfrac=0.3, objlim =5.0, 
                   gain =1.15,readnoise =6.5 ,  satlevel=65536, 
                   pssl =0.0,niter=4, sepmed=True ,cleantype ='meanmask' , fsmode='median',
                   psfmodel='gauss' ,  psffwhm =2.5, psfsize=7,
                   psfk=None ,  psfbeta =4.765 ,  verbose=False)
    return mask, _clean


