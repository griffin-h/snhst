from drizzlepac import astrodrizzle
from glob import glob
import numpy as np
from astropy.io import fits
from snhst.parameters import get_drizzle_parameters
from snhst.find_offsets import run_tweakreg
from snhst.cosmics_rays import detect_cosmic_rays


def drizzle(instrument, output_name, options):
    # Based on the instrument, get the default parameter values
    # Override any parameters that were provided by the user
    drizzle_options = get_drizzle_parameters(instrument, options)

    # If we want to find shifts
    # find the shifts comparing to the template image (either the first image of the set
    # or the provided template file)
    if drizzle_options['use_tweakshifts']:
        run_tweakreg(drizzle_options)

    # run astrodrizzle on the images
    run_astrodrizzle(output_name, drizzle_options)

    # Clean the cosmic rays from the image
    clean_cosmic_rays(drizzle_options)


def run_astrodrizzle(output_name, drizzle_options):

    num_images = len(glob('input_files'))
    if num_images == 1:
        combine_type = 'mean'
    else:
        combine_type = 'minmed'

    if drizzle_options['use_tweakshifts']:
        wcskey = 'TWEAK'
    else:
        wcskey = None

    # These parameters were either set empirically or are passed in by the user
    astrodrizzle.AstroDrizzle(drizzle_options['input_files'], output=output_name, runfile="",
                              wcskey=wcskey, context=True, group=drizzle_options['group'], build=True,
                              num_cores=drizzle_options['num_cores'], preserve=False,
                              clean=drizzle_options['clean'], skysub=drizzle_options['sky_subtract'],
                              skystat='mode', skylower=0.0, skyupper=None, driz_sep_fillval=-100000,
                              driz_sep_bits=drizzle_options['driz_bits'], driz_sep_wcs=True,
                              driz_sep_refimage=drizzle_options['template_image'],
                              driz_sep_rot=drizzle_options['rotation'],
                              driz_sep_scale=drizzle_options['pixel_scale'],
                              driz_sep_outnx=drizzle_options['nx'], driz_sep_outny=drizzle_options['ny'],
                              driz_sep_ra=drizzle_options['ra'], driz_sep_dec=drizzle_options['dec'],
                              combine_maskpt=0.2, combine_type=combine_type, combine_nsigma='4 3',
                              combine_nlow=0, combine_nhigh=0, combine_lthresh=-10000,
                              combine_hthresh=None, driz_cr=True, driz_cr_snr='3.5 3.0',
                              driz_cr_grow=1, driz_cr_ctegrow=0, driz_cr_scale='1.2 0.7',
                              final_pixfrac=drizzle_options['pixel_fraction'], final_fillval=-50000,
                              final_bits=drizzle_options['driz_bits'], final_units='counts',
                              final_wcs=True, final_refimage=drizzle_options['template_image'],
                              final_rot=drizzle_options['rotation'],
                              final_scale=drizzle_options['pixel_scale'],
                              final_outnx=drizzle_options['nx'], final_outny=drizzle_options['ny'],
                              final_ra=drizzle_options['ra'], final_dec=drizzle_options['dec'])

    # add the sky value back in
    # only add the sky to non flagged pixels (-50k)
    output_filename = glob(output_name + '_dr?.fits')[0]
    sci_hdu = fits.open(output_filename)
    sci_data = sci_hdu['SCI'].data
    if drizzle_options['sky_subtract']:
        mdrizsky = np.min(sci_data[sci_data > -49999.0])
    else:
        mdrizsky = 0.0
    sci_data[sci_data > -49999.0] -= mdrizsky
    no_data = sci_data < -49999.0
    sci_data[no_data] = 0.0

    sci_hdu.writeto(output_filename, clobber=True)


def clean_cosmic_rays(options):
    if options['instrument'] != 'wfc3_ir':
        images = glob('*_dr?.fits')
        for image in images:
            output_image = image[:-8] + 'cr' + image[-6:]
            detect_cosmic_rays(image, options['crpars'], output_image=output_image, masked_value=0.0)
