#!/usr/bin/env python
import os
import warnings
from pathlib import Path
import argparse
import numpy as np
import ccdproc as cdp
from astropy.stats import mad_std
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.wcs import FITSFixedWarning

warnings.simplefilter('ignore', FITSFixedWarning)

class Pipeline:
    def __init__(self, path):
        """
        Class constructor, it creates
        the parent ImageFileCollection object and
        performs basic data organisation

        Parameters
        ----------
        path: pathlib.Path object
        """
        self.path = path
        self.master_products_path = self.path / 'master_calibrations'

        if not os.path.exists(self.master_products_path):
            os.mkdir(self.master_products_path)

        self.parent_ifc = cdp.ImageFileCollection(path)
        self.darks_ifc = self.parent_ifc.filter(imagetyp='Dark Frame')
        self.bias_ifc = self.parent_ifc.filter(imagetyp='Bias Frame')
        self.flats_ifc = self.parent_ifc.filter(imagetyp='Flat Frame')
        self.light_ifc = self.parent_ifc.filter(imagetyp='Light Frame')


    def create_master_bias(self):
        """
        Following https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/02-04-Combine-bias-images-to-make-master.html
        section 2.3.3.1
        """
        bias_ccds = self.bias_ifc.ccds(ccd_kwargs=dict(unit='adu'))
        print('Creating master bias...')
        self.master_bias = cdp.combine(
            bias_ccds,
            method='average',
            sigma_clip=True,
            sigma_clip_low_thresh=5,
            sigma_clip_high_thresh=5,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            mem_limit=350e6
        )
        self.master_bias.meta['combined'] = True
        self.master_bias.write(self.master_products_path / 'master_bias.fits', overwrite=True)

    def create_master_dark(self):
        """
        Combines all available dark frames after master bias subtraction.
        TODO: - combine only darks with equal exposure length
              - identify hot pixels
        """
        dark_ccds = self.darks_ifc.ccds(ccd_kwargs=dict(unit='adu'))
        print('Creating master dark...')
        bias_subtracted_darks = list()
        for dark in dark_ccds:
            bias_subtracted_dark = cdp.subtract_bias(dark, self.master_bias)
            bias_subtracted_darks.append(bias_subtracted_dark)

        self.master_dark = cdp.combine(
            bias_subtracted_darks,
            method='average',
            sigma_clip=True,
            sigma_clip_low_thresh=5,
            sigma_clip_high_thresh=5,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            mem_limit=350e6
        )
        self.master_dark.meta['combined'] = True
        self.master_dark.write(self.master_products_path / 'master_dark.fits', overwrite=True)


        """
        Flat correction: bias substracion, assuming dark current = 0
        """
    
    
        def create_master_flat(self):
            flat_ccds = self.flats_ifc.ccds(ccd_kwargs=dict(unit='adu'))
            print('Creating master flat...')
            bias_subtracted_flats = list()   
            for flat in flat_ccds:
                bias_subtracted_flat = cdp.subtract_bias(flat, self.master_bias)
                bias_subtracted_flats.append(bias_subtracted_flat)
                
            self.master_flat = cdp.combine(
             bias_subtracted_flats,
             method='average',
             sigma_clip=True,
             sigma_clip_low_thresh=5,
             sigma_clip_high_thresh=5,
             sigma_clip_func=np.ma.median,
             sigma_clip_dev_func=mad_std,
             mem_limit=350e6
            )
            self.master_flat.meta['combined'] = True
            self.master_flat.write(self.master_products_path / 'master_flat.fits', overwrite=True)
    

def main():
    pwd = os.getcwd()
    print(f'Current directory: {pwd}')
    pipe = Pipeline(Path(pwd))
    pipe.create_master_bias()
    pipe.create_master_dark()
    pipe.create_master_flat()
    #print('This is script is under development.') #print('Right now it only prints this message.')
    #print('Soon it will serve to automatically reduce data from the MAS500 telescope.')
