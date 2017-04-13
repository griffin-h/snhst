import os
from snhst import fits_utils
from glob import glob
from astropy.io import fits
from snhst import parameters
from reproject import reproject_interp


def dolphot(template_image, images, options):
    # Copy the raw input files into a new directory
    if not os.path.exists('dolphot'):
        os.mkdir('dolphot')

    for image in images:
        fits_utils.copy_if_not_exists(image, 'dolphot')

    fits_utils.copy_if_not_exists(template_image, 'dolphot')

    os.chdir('dolphot')

    for image in images:
        instrument, detector = fits_utils.get_instrument(image).split('_')

        if needs_to_be_masked(image):
            mask_image(instrument, image)

        if needs_to_split_groups(image):
            split_groups(image)

    overlapping_images = get_overlapping_split_images(template_image, images)

    for image in overlapping_images:
        if needs_to_calc_sky(image):
            calc_sky(image, options)

    if need_to_run_dolphot():
        run_dolphot(template_image, images, options)

    os.chdir('..')


def needs_to_be_masked(image):
    # Masking should remove all of the DQ arrays etc, so make sure that any extensions with data in
    # in them are only SCI extensions. This might not be 100% robust, but should be good enough.
    hdulist = fits.open(image)
    needs_masked = False
    for hdu in hdulist:
        if hdu.data is not None and 'EXTNAME' in hdu.header:
            if hdu.header['EXTNAME'].upper() != 'SCI':
                needs_masked = True
    return needs_masked


def mask_image(instrument, image):
    os.system('{instrument}mask {image}'.format(instrument=instrument, image=image))


def needs_to_calc_sky(image):
    return not os.path.exists(image)


def calc_sky(image, options):
    calcsky_opts = parameters.get_calcsky_parameters(image, options)
    cmd = 'calcsky {image} {rin} {rout} {step} {sigma_low} {sigma_high}")'
    os.system(cmd.format(image=image, rin=calcsky_opts['r_in'], rout=calcsky_opts['r_out'],
                         step=calcsky_opts['step'], sigma_low=calcsky_opts['sigma_low'],
                         sigma_high=calcsky_opts['sigma_high']))


def needs_to_split_groups(image):
    return len(glob(image.replace('.fits', 'chip?.fits'))) == 0


def split_groups(image):
    os.system('splitgroups {filename}'.format(filename=image))


def get_overlapping_split_images(template_image, images):
    overlapping_images = []
    template_hdu = fits.open(template_image)
    for image in images:
        split_images = glob(image.replace('.fits', '.chip?.fits'))
        for split_image in split_images:
            header = fits.getheader(image)
            # remap the template_image onto the image
            _, footprint = reproject_interp(template_hdu, header)
            # If the overlap is at least 10%, consider the image to be overlapping
            if (footprint > 0).sum() >= (0.1 * footprint.size):
                overlapping_images.append(split_image)
    return overlapping_images


def need_to_run_dolphot():
    return not os.path.exists('dp.out')


def write_dolphot_image_parameters(file_object, image, i, options):
    file_object.write('img{i}_file = {file}\n'.format(i=i + 1, file=os.path.splitext(image)[0]))
    for par, value in parameters.get_dolphot_instrument_parameters(image, options).items:
        file_object.write('image{i}_{option} = {value}'.format(i=i + 1, option=par, value=value))


def write_dolphot_master_parameters(file_object, options):
    for par, value in options['dolphot'].items():
        file_object.write('{par} = {value}\n'.format(par=par, value=value))


def run_dolphot(template_image, images, options):
    f = open('dp.params','w')
    write_dolphot_master_parameters(f, options)
    f.write('Nimg = {n}\n'.format(len(images)))
    f.write('img0_file = {drzfile}\n'.format(drzfile=os.path.splitext(template_image)[0]))
    for i, image in enumerate(images):
        write_dolphot_image_parameters(f, image, i, options)
    f.close()
    os.system('dolphot dp.out -pdp.params 2 &> 1 | tee -a dp.log')
