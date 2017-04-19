#!/usr/bin/env python
import os
from glob import glob
import shutil
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import table

from snhst import drizzle, dolphot, wcs
import reproject
import sep

import logging
logging.basicConfig(level=logging.INFO)


def run(center=None, ground_reference=None, options=None):
    if options is None:
        options = {}

    # Make sure the position of interest is in the raw frame
    remove_images_without_object()

    # Copy the raw files into the raw directory
    copy_raw_data()

    # Download the reference images
    download_reference_images()

    # sort the data into instrument/visit/filter
    sort_raw_data()

    # Make the overall HST template image
    if reference_image is not None:
        options['ra'], options['dec'] = get_center_of_image(reference_image)

        hst_reference = get_default_reference_hst()
        # Drizzle the desired reference HST visit to the desired center
        os.chdir(os.path.join(hst_reference['instrument'], hst_reference['visit'], hst_reference['filter']))
        os.mkdir('match_template')
        os.chdir('match_template')
        shutil.copy('../*fl?.fits', './')
        # Run dolphot on the hst image to build the catalog

        drizzle.drizzle(hst_reference['instrument'], 'match_template', options)
        drizzled_template_filename = glob('match_template_dr?.fits')[0]

        # Resample the ground image to the same WCS as the HST frame
        reference_hdu = fits.open(reference_image)

        data, footprint = reproject.reproject_interp(reference_hdu, fits.getheader(drizzled_template_filename),)
        # Run sep on the ground image to make the ground catalog
        try:
            bkg = sep.Background(data, bw=32, bh=32, fw=3, fh=3)
        except ValueError:
            data = data.byteswap(True).newbyteorder()
            bkg = sep.Background(data, bw=32, bh=32, fw=3, fh=3)

        sources = sep.extract(data, options['sep']['threshold'], minarea=options['sep']['min_area'],
                              err=np.sqrt(data), deblend_cont=0.005)

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
        t['flux'] = flux
        t['fluxerr']
        t['mag'] = -2.5 * np.log10(t['flux'])
        t['magerr'] = np.log(10) / 2.5 * t['fluxerr'] / t['flux']

        t = t['x', 'y', 'mag', 'magerr']
        t.writeto('reference.cat', format='fast_no_header')

        # Run dolphot on the HST image
        dolphot.dolphot(drizzled_template_filename, options['default_images'])

        wcs.parse_dolphot_table('dolphot/dp.out', drizzled_template_filename, 'hst.cat')
        # Calculate the offset between the frames (and apply them)
        offsets = wcs.calculate_wcs_offset(drizzled_template_filename, 'hst.cat',
                                           drizzled_template_filename, 'reference.cat')

        # Update the raw hst frames with the calculated offset
        os.mkdir('../hst_template')
        os.chdir('../hst_template')
        shutil.copy('../*fl?.fits', './')
        wcs.apply_offsets(glob('*fl?.fits', offsets))

        drizzle.drizzle(hst_reference['instrument'], 'hst_template', options)
        shutil.copy('*dr?.fits', topdir)
    else:
        # Drizzle the reference image
        os.chdir(os.path.join(hst_reference['instrument'], hst_reference['visit'], hst_reference['filter']))
        hst_reference = get_default_reference_hst()
        drizzle.drizzle(hst_reference['instrument'], 'hst_template', options)
        dolphot.dolphot(drizzled_template_filename, options['default_images'])
        shutil.copy('hst_template_dr?.fits', topdir)

    instruments, visits, filters = get_visit_info(images)
    # For each instrument
    for instrument in instruments:
        # Drizzle the desired reference visit given the desired pixel scale and image size
        instrument_template = get_instrument_template_visit(instrument, visits[insturment], filters[instrument])
        drizzle_and_match(hst_reference, options)

    for instrument in instruments:
        for visit in visits:
            for f in filters:
                drizzle_and_match(insturment_reference, options)
                # Copy the files into a master dolphot directory
                shutil.copy("*fl?.fits", topdir +'/dolphot/')

    # Run dolphot on the master dolphot directory
    os.chdir(topdir + '/dolphot/')
    dolphot.dolphot(new_hst_template, '*')

    # Convert the dolphot catalog to a fits file
    dolphot.parse_catalog()

    # Run dolphot again with fake stars
    dolphot.add_fake_stars()
    # Copy the final drizzled frames, cosmic ray removed frames, and the dolphot catalog to
    # a top directory
    sort_final_data()


def drizzle_and_match(template_image, options):
    # Run dolphot to make the photometry catalog for each visit/filter
    # Calculate the offset between the drizzled frame and the reference visit
    # update the raw files
    # redrizzle the images with the updated offsets
    pass


def remove_images_without_object(ra, dec):
    if not os.path.exists('unused'):
        os.mkdir('unused')

    images = get_raw_image_filenames(raw_directory='.')
    for image in images:
        if not image_includes_coordinate(image, ra, dec):
            shutil.move(image, 'unused/')


def image_includes_coordinate(image, ra, dec):
    image_hdu = fits.open(image)
    for hdu in image_hdu:
        if hdu.header['EXTNAME'] == 'SCI':
            image_wcs = WCS(hdu.header)

            # Use zero indexed here for the pixel coordinates
            coordinate_in_pixels = image_wcs.all_world2pix(ra, dec, 0)

            coordinate_in_image = False
            coordinate_not_negative = np.all(coordinate_in_pixels > 0)
            coordinate_not_outside = np.all(coordinate_in_pixels >= np.array([hdu.header['NAXIS2'],
                                                                              hdu.header['NAXIS1']]))
            if coordinate_not_negative and coordinate_not_outside:
                coordinate_in_image = True

    image.close()
    return coordinate_in_image


def copy_raw_data():
    if not os.path.exists('raw'):
        os.mkdir('raw')
    fs = glob('*.fits')
    for f in fs:
        shutil.move(f, 'raw/')


def download_reference_images():
    if not os.path.exists('ref'):
        os.mkdir('ref')
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH'] = os.getcwd() + '/ref/'
    os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files raw/*flc.fits')
    os.system('python -m crds.bestrefs --update-bestrefs --sync-references=1 --files raw/*c0m.fits')


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


def get_raw_image_filenames(raw_directory='./raw'):
    images = []
    for ftype in ['flc', 'c0m']:
        images += glob(os.path.join(raw_directory, '*{0}.fits'.format(ftype)))
    return images


def within_half_day(expstart1, expstart2):
    return np.abs(expstart1 - expstart2) > 0.5


def get_unique_visits(images):
    exposure_start_mjds = [float(fits.getval(image, 'EXPSTART')) for image in images]
    sorted_start_inds = np.argsort(exposure_start_mjds)
    exposure_start_mjds = np.array(exposure_start_mjds)[sorted_start_inds]
    instruments = np.array([fits.getval(image, 'INSTRUME').lower() for image in images])[sorted_start_inds]
    detectors = np.array([fits.getval(image, 'DETECTOR').lower() if 'flc.fits' in image else "wfpc2"
                          for image in images])[sorted_start_inds]

    unique_exposure_visits = [(x, instruments[i], detectors[i]) for i, x in
                              enumerate(exposure_start_mjds[:-1])
                              if within_half_day(x, exposure_start_mjds[i + 1]) or detectors[i] != detectors[i + 1]]
    unique_exposure_visits.append((exposure_start_mjds[-1], instruments[-1], detectors[-1]))
    return unique_exposure_visits


def sort_raw_data():
    images = get_raw_image_filenames()

    visits = get_unique_visits(images)

    for visit_start, instrument, detector in visits:
        if not os.path.exists(instrument):
            os.mkdir(instrument)
        if not os.path.exists(os.path.join(instrument, detector)):
            os.mkdir(os.path.join(instrument, detector))

        num_previous_visits = len(glob(os.path.join(instrument, detector, 'visit*')))
        visit_folder = os.path.join(instrument, detector, 'visit%i' % (num_previous_visits + 1))

        if not os.path.exists(visit_folder):
            os.mkdir(visit_folder)

        for f in get_raw_image_filenames():
            same_visit = np.abs(float(fits.getval(f, 'EXPSTART')) - visit_start) < 0.5
            if "flc.fits" in f:
                same_visit &= fits.getval(f, 'DETECTOR').lower() == detector
            else:
                same_visit &= detector == 'wfpc2'
            same_visit &= fits.getval(f, 'INSTRUME').lower() == instrument
            if same_visit:
                folder_name = os.path.join(visit_folder, get_filter_name(f))
                if not os.path.exists(folder_name):
                    os.mkdir(folder_name)
                shutil.copy(f, folder_name + '/')


if __name__ == '__main__':
    run()
