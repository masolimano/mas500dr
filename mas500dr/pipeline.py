#!/usr/bin/env python
import os
from pathlib import Path
import argparse
import numpy as np
import ccdproc as cdp
from astropy.stats import mad_std
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits

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
        self.master_bias.write(self.master_products_path / 'master_bias.fits')
#    def create_master_dark(self):
 #       pass

    def create_master_dark(self):
        dark_ccds = self.dark_ifc.ccds(ccd_kwargs=dict(unit='adu'))
        print('Creating master dark...')
        self.master_dark = cdp.combine(
            dark_ccds,
            method='average',
            sigma_clip=True,
            sigma_clip_low_thresh=5,
            sigma_clip_high_thresh=5,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            mem_limit=350e6
        )
        self.master_dark.meta['combined'] = True
        self.master_dark.write(self.master_products_path / 'master_dark.fits')


"""
This function is not similar as the one of the master_bias
"""
    def create_master_dark():
        """
        Following https://ccdproc.readthedocs.io/en/stable/reduction_toolbox.html#subtract-bias-and-dark
        """
        keys = ['file', 'imagetyp', 'object', 'filter', 'exposure']
        ic1 = ImageFileCollection('path/to/your/calibration_folder', keywords=keys) # only keep track of keys

        ic_dark = ImageFileCollection('calibration_dataset', glob_include='Dark*')
        print('Creating master bias...')
        matches_dark = (ic_dark.summary['imagetyp'] == 'Dark Frame') 
        darks = ic_dark.summary['file'][matches_dark]
        master_dark = CCDData(darks, unit=u.adu)

        return master_dark



def main():
    pwd = os.getcwd()
    print(f'Current directory: {pwd}')
    pipe = Pipeline(Path(pwd))
    pipe.create_master_bias()
    #pipe.create_master_dark()
    #print('This is script is under development.')
    #print('Right now it only prints this message.')
    #print('Soon it will serve to automatically reduce data from the MAS500 telescope.')
