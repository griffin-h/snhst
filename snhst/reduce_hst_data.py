#!/usr/bin/env python
import os
from glob import glob
import shutil
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import table

from snhst import drizzle, dolphot, wcs, fits_utils
import reproject
import sep

import logging
logging.basicConfig(level=logging.INFO)


def run(center=None, ground_reference=None, options=None):
    if options is None:
        options = {'sep': {'thresh': 5.},
                   'dolphot': {},
                   'dolphot_img': {},
                   'dolphot_sky': {}}

    global topdir
    topdir = os.getcwd()

    # Unzip images
    if glob('*.fits.gz'):
        os.system('gunzip *.fits.gz')

    # Download the reference images
    if not os.path.exists('ref/'):
        download_reference_images()

    # Copy the raw files into the raw directory
    if not os.path.exists('raw/'):
        copy_raw_data()

    # output coords from header of ground reference, if not given
    if center is None and ground_reference is not None:
        center['ra'], center['dec'] = get_target_coordinates(ground_reference)

    raw_images = get_raw_image_filenames()

    # Make sure the position of interest is in the raw frame
    if center is not None:
        remove_images_without_object(center['ra'], center['dec'], raw_images)
        options['ra'], options['dec'] = center['ra'], center['dec']

    # sort the data into instrument/visit/filter
    visit_meta_data = sort_raw_data(raw_images)
    hst_reference = get_default_reference_hst(visit_meta_data)
    instrument_string = go_to_visit(hst_reference)

    # Make the overall HST template image
    if ground_reference is not None:

        drizzled_template_filename, drizzled_template_catalog = match_and_drizzle('match_template', instrument_string, options)
        drizzled_template_header = fits.getheader(drizzled_template_filename, extname='SCI')

        # Get the science hdu of the ground reference (may be compressed)
        reference_hdulist = fits.open(os.path.join(topdir, ground_reference))
        for reference_hdu in reference_hdulist:
            if reference_hdu.header.get('NAXIS'):
                break

        if not os.path.exists('ground_reference.cat'):
            # Reproject ground image onto drizzled HST frame
            data, footprint = reproject.reproject_interp(reference_hdu, drizzled_template_header)
            # Run sep on the ground image to make the ground catalog
            make_sep_catalog(data, options['sep'], 'match_template/ground_reference.cat')

        # Calculate the offset between the frames
        # (both catalogs are projected to the WCS of the drizzled template)
        offsets = wcs.calculate_wcs_offset(drizzled_template_header, drizzled_template_catalog, 'match_template/ground_reference.cat', max_offset=2.)
    else:
        offsets = None

    hst_reference_filename, hst_reference_catalog = match_and_drizzle('hst_template', instrument_string, options, offsets)

    # drizzle each visit separately
    for visit in visit_meta_data:
        if (visit == hst_reference).all():
            continue
        instrument_string = go_to_visit(visit)
        visit_template_filename, visit_template_catalog = match_and_drizzle('match_template', instrument_string, options)
        visit_template_header = fits.getheader(visit_template_filename, extname='SCI')
        offsets = wcs.calculate_wcs_offset(visit_template_header, visit_template_catalog, hst_reference_catalog, max_offset=2.)
        match_and_drizzle('hst_template', instrument_string, options, offsets)

    # drizzle all visits together
    output_images = visit_meta_data.group_by(['instrument', 'detector', 'subarray', 'filter'])
    for output_image in output_images.groups.keys:
        instrument_string = go_to_visit(output_image)
        for image in get_raw_image_filenames('visit*'):
            fits_utils.copy_if_not_exists(image, '.')
        filter_template_filename, filter_template_catalog = match_and_drizzle('match_template', instrument_string, options)
        filter_template_header = fits.getheader(filter_template_filename, extname='SCI')
        offsets = wcs.calculate_wcs_offset(filter_template_header, filter_template_catalog, hst_reference_catalog, max_offset=2.)
        match_and_drizzle('hst_template', instrument_string, options, offsets)

    # TODO: THESE ARE NOT DEFINED YET
    # # Run dolphot again with fake stars
    # dolphot.add_fake_stars()

    # # Copy the final drizzled frames, cosmic ray removed frames, and the dolphot catalog to
    # # a top directory
    # sort_final_data()


def go_to_visit(visit):
    instrument_string = '_'.join(visit[['instrument', 'detector']])
    if visit['subarray']:
        instrument_string += '_sub'
    else:
        instrument_string += '_full'
    visit_path = os.path.join(topdir, visit['instrument'], visit['detector'], visit['filter'])
    if 'visit' in visit.colnames:
        visit_path = os.path.join(visit_path, 'visit{:d}'.format(visit['visit']))
    os.chdir(visit_path)

    return instrument_string


def mkdir_p(path):
    '''Makes a directory if it does not already exist. Equivalent to bash `mkdir -p`.'''
    if not os.path.exists(path):
        os.makedirs(path)


def make_sep_catalog(data, options, output_catalog):
    try:
        bkg = sep.Background(data, bw=32, bh=32, fw=3, fh=3)
    except ValueError:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data, bw=32, bh=32, fw=3, fh=3)

    error = np.sqrt(data)
    sources = sep.extract(data, err=error, deblend_cont=0.005, **options)

    t = table.Table(sources)
    kronrad, krflag = sep.kron_radius(data, sources['x'], sources['y'],
                                      sources['a'], sources['b'],
                                      sources['theta'], 6.0)

    flux, fluxerr, flag = sep.sum_ellipse(data, sources['x'], sources['y'],
                                          sources['a'], sources['b'],
                                          np.pi / 2.0, 2.5 * kronrad,
                                          subpix=1, err=error)

    t['x'] += 1
    t['y'] += 1
    t['mag'] = -2.5 * np.log10(flux)
    t['magerr'] = np.log(10) / 2.5 * fluxerr / flux

    t = t['x', 'y', 'mag', 'magerr']
    t.write(output_catalog, format='ascii.fast_no_header')


def match_and_drizzle(basename, instrument_string, options, offsets=None):
    initial_directory = os.getcwd()

    # Update the raw hst frames with the calculated offset
    mkdir_p(basename)
    for image in get_raw_image_filenames('.'):
        fits_utils.copy_if_not_exists(image, basename)
    os.chdir(basename)
    if offsets is not None:
        wcs.apply_offsets(get_raw_image_filenames('.'), offsets)

    # Drizzle the reference image
    drizzled_templates = glob(basename + '_dr?.fits')
    if not drizzled_templates:
        drizzle.drizzle(instrument_string, basename, options)
        drizzled_templates = glob(basename + '_dr?.fits')
    drizzled_template_filename = os.path.abspath(drizzled_templates[0])

    # Run dolphot on the hst image to build the catalog
    dolphot.dolphot(drizzled_template_filename, get_raw_image_filenames('.'), options)
    output_cat = os.path.abspath(basename + '.cat')
    wcs.parse_dolphot_table('dolphot/dp.out', output_cat)

    os.chdir(initial_directory)

    return drizzled_template_filename, output_cat


def remove_images_without_object(ra, dec, images):
    mkdir_p('unused')
    for image in images:
        if not image_includes_coordinate(image, ra, dec):
            shutil.move(image, 'unused/')


def image_includes_coordinate(image, ra, dec):
    image_hdu = fits.open(image)
    coordinate_in_image = False
    for hdu in image_hdu:
        if hdu.header.get('EXTNAME') == 'SCI':
            image_wcs = WCS(hdu.header)

            # Use zero indexed here for the pixel coordinates
            coordinate_in_pixels = image_wcs.all_world2pix([[ra, dec]], 0)
            coordinate_not_negative = np.all(coordinate_in_pixels > 0)
            coordinate_not_outside = np.all(coordinate_in_pixels < np.array([hdu.header['NAXIS2'],
                                                                             hdu.header['NAXIS1']]))
            if coordinate_not_negative and coordinate_not_outside:
                coordinate_in_image = True
                break

    image_hdu.close()
    return coordinate_in_image


def copy_raw_data():
    mkdir_p('raw')
    fs = glob('*.fits')
    for f in fs:
        fits_utils.copy_if_not_exists(f, 'raw/')


def download_reference_images():
    mkdir_p('ref')
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH'] = os.getcwd() + '/ref/'
    if glob('*flc.fits'):
        os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files *flc.fits')
    if glob('*c0m.fits'):
        os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files *c0m.fits')


def get_filter_name(image):
    if 'c0m.fits' in image:
        f = fits.getval(image, 'FILTNAM1')
        if len(f) == 0:
            f = fits.getval(image, 'FILTNAM2')
    else:
        try:
            f = fits.getval(image, 'FILTER')
        except:
            f = fits.getval(image, 'FILTER1')
            if 'clear' in f.lower():
                f = fits.getval(image, 'FILTER2')
    return f.lower()


def get_raw_image_filenames(raw_directory='raw'):
    for ftype in ['flc', 'flt', 'c0m']:
        images = glob(os.path.join(raw_directory, '*{0}.fits'.format(ftype)))
        if images:
            break
    return images


def get_default_reference_hst(visit_meta_data):
    preferred_filters = ['f555w', 'f625w', 'f606w', 'f622w', 'f814w',
                         'f791w', 'f438w', 'f435w', 'f110w', 'f160w']
    for filt in preferred_filters:
        if filt in visit_meta_data['filter']:
            visit = visit_meta_data[visit_meta_data['filter'] == filt][0]
            break
    else: # choose the central filter
        visit_meta_data.sort('filter')
        visit = visit_meta_data[len(visit_meta_data) // 2]
    return visit


def get_target_coordinates(flt_file):
    hdulist = fits.open(flt_file)
    RA = hdulist[0].header['RA_TARG']
    Dec = hdulist[0].header['DEC_TARG']
    return [RA, Dec]


def get_unique_visits(images):
    exposure_start_mjds = [float(fits.getval(image, 'EXPSTART')) for image in images]
    sorted_start_inds = np.argsort(exposure_start_mjds)
    exposure_start_mjds = np.array(exposure_start_mjds)[sorted_start_inds]
    instruments = np.array([fits.getval(image, 'INSTRUME').lower() for image in images])[sorted_start_inds]
    detectors = np.array([fits.getval(image, 'DETECTOR').lower() if 'flc.fits' in image else "wfpc2"
                          for image in images])[sorted_start_inds]
    # group images by start time & instrument
    unique_exposure_visits = [(x, instruments[i], detectors[i]) for i, x in
                              enumerate(exposure_start_mjds[:-1])
                              if np.abs(x - exposure_start_mjds[i + 1]) > min_visit_separation
                              or detectors[i] != detectors[i + 1]]
    unique_exposure_visits.append((exposure_start_mjds[-1], instruments[-1], detectors[-1]))
    return unique_exposure_visits


def sort_raw_data(images, min_visit_separation=0.2):
    metadata = table.Table()
    metadata['filename'] = images
    metadata['expstart'] = [fits.getval(image, 'EXPSTART') for image in images]
    metadata['instrument'] = np.char.lower([fits.getval(image, 'INSTRUME') for image in images])
    metadata['detector'] = np.char.lower([fits.getval(image, 'DETECTOR') for image in images])
    metadata['filter'] = [get_filter_name(image) for image in images]
    metadata['subarray'] = [fits.getval(image, 'SUBARRAY') for image in images]
    metadata['visit'] = np.ones(len(metadata), int)
    metadata.sort('expstart')
    visit_starts, = np.where(np.diff(metadata['expstart']) > min_visit_separation)
    for visit_num, first_row in enumerate(visit_starts):
        metadata['visit'][first_row + 1:] = visit_num + 2

    for image in metadata:
        folder_name = os.path.join(image['instrument'], image['detector'], image['filter'],
                                   'visit{:d}'.format(image['visit']))
        mkdir_p(folder_name)
        fits_utils.copy_if_not_exists(image['filename'], folder_name)

    visit_list = metadata.group_by(['instrument', 'detector', 'filter', 'visit', 'subarray'])

    return visit_list.groups.keys


if __name__ == '__main__':
    run()
