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
from astropy.utils.exceptions import AstropyUserWarning

from .utils import inv_median, create_ccd_mask,cosmic_ray_correction

warnings.simplefilter('ignore', FITSFixedWarning)


class Pipeline:
    def __init__(self, path, mem_limit=350e6):
        """
        Class constructor, it creates
        the parent ImageFileCollection object and
        performs basic data organisation

        Parameters
        ----------
        path: pathlib.Path object

        mem_limit: float
            Memory limit for image combination
        """
        self.root = path
        self.mem_limit = mem_limit

        self.raw_path = self.root / 'raw'

        if not self.raw_path.exists():
            self.raw_path.mkdir()

        fits_files_in_root = self.root.glob('*.fit*')
        if len(list(fits_files_in_root)) > 0:
            for filename in fits_files_in_root:
                filename.rename(self.raw_path / filename.name)

        self.master_products_path = self.root / 'master_calibrations'
        self.calibrated_path = self.root / 'calibrated'

        if not self.master_products_path.exists():
            self.master_products_path.mkdir()

        if not self.calibrated_path.exists():
            self.calibrated_path.mkdir()

        self.parent_ifc = cdp.ImageFileCollection(self.raw_path)
        assert len(self.parent_ifc.files) > 0
        self.darks_ifc = self.parent_ifc.filter(imagetyp='Dark Frame')
        self.bias_ifc = self.parent_ifc.filter(imagetyp='Bias Frame')
        self.flats_ifc = self.parent_ifc.filter(imagetyp='Flat Frame')
        self.light_ifc = self.parent_ifc.filter(imagetyp='Light Frame')

        light_pd = self.light_ifc.summary.to_pandas()
        self.grouped_light = light_pd.groupby(by=['filter', 'xbinning', 'readoutm'])


    @staticmethod
    def _avg_combine(frame_list, mem_limit=350e6, **kwargs):
        """
        Wrapper of ccdproc.combine with the keyword values recommended
        in the ccdproc book for use with averaging.

        Parameters
        ----------
        frame_list: generator of CCDData
            Generator that yields the frames to combine

        mem_limit: float
            Memory limit for image combination

        Returns
        ----------
        combined_ccd: CCDData object
        """
        combined_ccd = cdp.combine(
            frame_list,
            method='average',
            sigma_clip=True,
            sigma_clip_low_thresh=5,
            sigma_clip_high_thresh=5,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            mem_limit=mem_limit,
            **kwargs
        )
        combined_ccd.meta['combined'] = True
        return combined_ccd

    @staticmethod
    def _median_combine(frame_list, mem_limit=350e6, **kwargs):
        """
        Wrapper of ccdproc.combine set to median stacking
        Use this when the length of frame_list is not too large.

        Parameters
        ----------
        frame_list: generator of CCDData
            Generator that yields the frames to combine

        mem_limit: float
            Memory limit for image combination

        Returns
        ----------
        combined_ccd: CCDData object
        """
        combined_ccd = cdp.combine(
            frame_list,
            method='median',
            mem_limit=mem_limit,
            **kwargs
        )
        combined_ccd.meta['combined'] = True
        return combined_ccd



    def create_master_bias(self):
        """
        Following https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/02-04-Combine-bias-images-to-make-master.html
        section 2.3.3.1
        """
        bias_ccds = self.bias_ifc.ccds(ccd_kwargs=dict(unit='adu'))
        print('Creating master bias...')
        if len(self.bias_ifc.files) < 3:
            raise RuntimeError('Not enough bias frames to combine')
        elif 3 <= len(self.bias_ifc.files) < 6:
            self.master_bias = self._median_combine(bias_ccds, self.mem_limit)
        else:
            self.master_bias = self._avg_combine(bias_ccds, self.mem_limit)

        binning = self.master_bias.header['XBINNING']
        readout = self.master_bias.header['READOUTM'].replace(' ', '')
        self.master_bias.write(self.master_products_path / f'master_bias_bin{binning}_{readout}.fits', overwrite=True)

    def create_master_dark(self):
        """
        Combines all available dark frames after master bias subtraction.
        TODO: - combine only darks with equal:
                - exposure length
                - binning
                - readout mode
              - identify hot pixels
        """
        available_exptimes = np.unique(self.darks_ifc.summary['exptime'].data)
        self.master_dark = dict()
        for exptime in available_exptimes:
            dark_collection = self.darks_ifc.filter(exptime=exptime, readoutm='1 MHz')
            dark_ccds = dark_collection.ccds(ccd_kwargs=dict(unit='adu'))
            print(f'Creating master dark of t_exp = {exptime:.1f} seconds')

            if len(dark_collection.files) < 3:
                raise RuntimeError('Not enough dark frames to combine')
                continue
            elif 3 <= len(dark_collection.files) < 6:
                self.master_dark[exptime] = self._median_combine(dark_ccds, self.mem_limit)
            else:
                self.master_dark[exptime] = self._avg_combine(dark_ccds, self.mem_limit)


            binning = self.master_dark[exptime].header['XBINNING']
            self.master_dark[exptime].write(self.master_products_path / f'master_dark_{exptime:.0f}s_bin{binning:d}_1MHz.fits',
                                            overwrite=True)




    def create_master_flat(self):
        """
        Flat correction: bias substracion, assuming dark current = 0
        """
        flat_ccds = self.flats_ifc.ccds(ccd_kwargs=dict(unit='adu'))
        print('Creating master flat...')
        bias_subtracted_flats = list()
        for flat in flat_ccds:
            bias_subtracted_flat = cdp.subtract_bias(flat, self.master_bias)
            bias_subtracted_flats.append(bias_subtracted_flat)
        self.master_flat = self._avg_combine(bias_subtracted_flats, self.mem_limit, scale=inv_median)


        # Saving
        filter_ = self.master_flat.header['FILTER']
        self.master_flat.write(self.master_products_path / f'master_flat_{filter_}.fits', overwrite=True)

    def calibrate_science(self):
        """
        Calibrate science exposures. Output is in units of electron
        """

        for raw_science, fname in self.light_ifc.ccds(return_fname=True, ccd_kwargs=dict(unit='adu')):
            exptime = raw_science.header['EXPTIME']
            calib_science = cdp.ccd_process(
                raw_science,
                dark_frame=self.master_dark[exptime],
                master_flat=self.master_flat,
                gain=1.5*u.electron/u.adu,
                readnoise=10*u.electron,
                exposure_key='exptime',
                exposure_unit=u.second,
                gain_corrected=False
            )
            calib_science.write(self.calibrated_path / f'{fname[:-4]}_calibrated.fits', overwrite=True)

    def run(self):
        """
        Work in progress: This function will iterate over all the dataset finding and applying the corresponding
        calibrations
        """
        self.group_iterator = self.grouped_light.__iter__()
        self.calib_ifc  = cdp.ImageFileCollection(self.master_products_path)
        for (flt, binn, rdm), group in self.group_iterator:
            print(flt, binn, rdm)
            if len(self.calib_ifc.files) > 0:
                pre_calibrated_bias = self.calib_ifc.filter(imagetyp='Bias Frame',  xbinning=binn, readoutm=rdm)
                if len(pre_calibrated_bias.files) == 0:
                    print('No processed bias for the current group. Selecting and ')
                    self.bias_ifc = self.parent_ifc.filter(imagetyp='Bias Frame', xbinning=binn, readoutm=rdm)
                    self.create_master_bias()
                    pre_calibrated_bias.refresh()
                elif len(pre_calibrated_bias.files) == 1:
                    print(f'Processed  bias found for the current group. Loading {pre_calibrated_bias.files[0]}')
                    self.master_bias = CCDData.read(pre_calibrated_bias.files[0])
                else:
                    raise RuntimeError("More than one master bias, can't decide")
            else:
                print(f'Empty {self.master_products_path.name} folder')
                self.bias_ifc = self.parent_ifc.filter(imagetyp='Bias Frame', xbinning=binn, readoutm=rdm)
                self.create_master_bias()

            self.calib_ifc.refresh()
            self.darks_ifc = self.parent_ifc.filter(imagetyp='Dark Frame', xbinning=binn, readoutm=rdm)
            pre_calibrated_darks = self.calib_ifc.filter(imagetyp='Dark Frame', combined=True, xbinning=binn)
            if len(pre_calibrated_darks.files) > 0:
                pass


            self.flats_ifc = self.parent_ifc.filter(imagetyp='Flat Frame', xbinning=binn, filter=flt)
            pre_calibrated_flats = self.calib_ifc.filter(imagetyp='Flat frame', filter=flt, xbinning=binn)
            if len(pre_calibrated_flats.files) == 0:
                pass


def main():
    """
    User-facing script to automatically reduce all science
    data in the current directory
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--mem_limit", help="Maximum memory (bytes) to use in image combination", type=float, default=350e6)
    args = parser.parse_args()
    pwd = os.getcwd()
    print(f'Current directory: {pwd}')
    pipe = Pipeline(Path(pwd), mem_limit=args.mem_limit)
    pipe.create_master_bias()
    pipe.create_master_dark()
    pipe.create_master_flat()
    pipe.calibrate_science()
