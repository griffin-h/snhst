from drizzlepac import tweakreg
from glob import glob
import shutil
import os
from astropy.io import fits

from snhst.cosmics_rays import detect_cosmic_rays
from snhst.fits_utils import is_sci_extension


def run_tweakreg(options):
    images = glob(options['input_files'])

    if options['instrument'] != 'wfc3_ir':
        # Run cosmic ray rejection on the input images to make registering the images easier
        for image in images:
            # Make a copy of the original data file to use later
            shutil.copy(image, make_raw_tmp_filename(image))
            # Remove cosmic rays, this operates on the files in place
            detect_cosmic_rays(image, options['crpars'])

    tweakreg.TweakReg(files=options['input_files'], refimage=options['template_image'],
                      interactive=False, writecat=False, clean=True, updatehdr=True,
                      wcsname='TWEAK', reusename=True, rfluxunits='counts', see2dplot=False,
                      separation=0.5, residplot="No plot", runfile='',
                      imagefindcfg={'threshold': options['tweakshifts_threshold'],
                                    'use_sharp_round': True},
                      refimagefindcfg={'threshold': options['tweakshifts_threshold'],
                                       'use_sharp_round': True})

    if options['instrument'] != 'wfc3_ir':
        for image in images:
            # copy the raw data back into the input file with the updated shift
            copy_data(make_raw_tmp_filename(image), image)
            # Clean up
            os.remove(make_raw_tmp_filename(image))


def make_raw_tmp_filename(filename):
    return filename[:-8] + 'rawtmp.fits'


def copy_data(image_data_from, image_data_to):
    hdu_from = fits.open(image_data_from)

    hdu_to = fits.open(image_data_to)

    for i, hdu in enumerate(hdu_from):
        if is_sci_extension(hdu):
            hdu_to[i].data[:, :] = hdu.data[:, :]

    hdu_from.close()

    hdu_to.writeto(image_data_to, overwrite=True)
    hdu_to.close()
