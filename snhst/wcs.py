import numpy as np
from astropy import coordinates, units
from astropy import table
from astropy import wcs
from astropy.io import fits
from scipy import optimize


def calculate_wcs_offset(hst_table, ref_table, n_stars=120, max_offset=0.1):

    # Choose the brightest objects that overlap both images (always assume the reference image
    # and the HST image have the same WCS and are the same shape)
    ref_table.sort('mag')
    ref_table = ref_table[:n_stars]

    hst_table.sort('mag')
    hst_table = hst_table[:n_stars]

    # Match the brighest objects between frames (the closest thing within ~0.1 arcsec)
    hst_skycoords = coordinates.SkyCoord(hst_table['ra'], hst_table['dec'], unit=('deg', 'deg'))
    ref_skycoords = coordinates.SkyCoord(ref_table['ra'], ref_table['dec'], unit=('deg', 'deg'))

    match_id, separations, _ = hst_skycoords.match_to_catalog_sky(ref_skycoords)
    ref_table = ref_table[match_id][separations.arcsec < max_offset]
    hst_table = hst_table[separations.arcsec < max_offset]

    # set the initial guess to be what is already in the header of the HST image
    results = optimize.minimize(wcs_objective, [0, 0], args=(hst_table, ref_table), method='Nelder-Mead')
    return results['x']


def wcs_objective(offset, hst_table, reference_catalog):
    ra_offset, dec_offset = offset
    offset_skycoords = coordinates.SkyCoord(hst_table['ra'] + ra_offset, 
                                            hst_table['dec'] + dec_offset, unit=('deg', 'deg'))
    ref_skycoords = coordinates.SkyCoord(reference_catalog['ra'], reference_catalog['dec'], unit=('deg', 'deg'))
    _, separations, _ = offset_skycoords.match_to_catalog_sky(ref_skycoords)
    # Minimize the mean offset between all of the matched sources
    return np.mean(separations.arcsec)


def read_catalog(catalog_filename):
    return table.Table.read(catalog_filename, format='ascii.fast_no_header', names=['x', 'y', 'mag', 'magerr'])


def apply_offsets(images, offset):
    '''Update the raw hst frames with the calculated offset'''
    for image in images:
        hdulist = fits.open(image)
        for hdu in hdulist:
            hdr = hdu.header
            if hdr.get('EXTNAME') == 'SCI':
                hdr['CRVAL1'] += offset[0]
                hdr['CRVAL2'] += offset[1]
        hdulist.writeto(image, overwrite=True)


def offset_to_match(images, hst_header, hst_table, ref_table,
                    n_stars=120, max_offset=0.1, origin=1):
    # Make WCS object that applies to both catalogs
    final_wcs = wcs.WCS(hst_header)

    # Read in the ref_table x's, y's and convert to ra, dec
    if 'ra' not in ref_table.colnames or 'dec' not in ref_table.colnames:
        ref_table['ra'], ref_table['dec'] = final_wcs.all_pix2world(ref_table['x'], ref_table['y'], origin)
    hst_table['ra'], hst_table['dec'] = final_wcs.all_pix2world(hst_table['x'], hst_table['y'], origin)

    offsets = calculate_wcs_offset(hst_table, ref_table, n_stars, max_offset)

    apply_offsets(images, offsets)