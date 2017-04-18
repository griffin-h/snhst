import numpy as np
from astropy import coordinates, units
from astropy import table
from astropy import wcs
from scipy import optimize

from snhst.dolphot import cut_bad_dolphot_sources


def calculate_wcs_offset(input_image, input_catalog, reference_image, reference_catalog,
                         n_stars=120, max_offset=0.1, origin=1):
    # Read in the reference catalog and the input catalog x's, y's
    catalog = read_catalog(reference_image, reference_catalog)
    # Convert the reference x and y's into ra and dec
    catalog['ra'], catalog['dec'] = wcs.WCS(reference_image).all_pix2world(ref_x, ref_y, origin)

    input_catalog = read_catalog(input_image, input_catalog)
    input_wcs = wcs.WCS(input_image)
    input_catalog['ra'], input_catalog['dec'] = input_wcs.all_pix2world(input_catalog['x'], input_catalog['y'], origin)

    # Choose the brightest objects that overlap both images (always assume the reference image
    # is larger than the input image
    in_input_image = coords_inside_image(input_wcs, catalog)
    catalog =  catalog[in_input_image].sorted('mag', desc=True)[:n_stars]

    input_catalog = input_catalog.sorted('mag', desc=True)[:n_stars]

    # Match the brighest objects between frames (the closest thing within ~0.1 arcsec)
    input_skycoords = coordinates.SkyCoord(input_catalog['ra'], input_catalog['dec'])
    ref_skycoords = coordinates.SkyCoord(catalog['ra'], catalog['dec'])

    match_id, separations, _ = input_skycoords.match_to_catalog_sky(ref_skycoords)
    catalog = catalog[match_id][separations.to(units.arcsec) < max_offset]
    input_catalog = input_catalog[separations.to(units.arcsec) < max_offset]

    # set the initial guess to be what is already in the header of the input catalog
    optimize.minimize(wcs_objective, [0, 0], args=(input_catalog, catalog), method='Nelder-Mead')


def wcs_objective(offset, input_catalog, reference_catalog):
    ra_offset, dec_offset = offset
    offset_skycoords = coordinates.SkyCoord(input_catalog['ra'] + ra_offset, input_catalog['dec'] + dec_offset)
    ref_skycoords = coordinates.SkyCoord(reference_catalog['ra'], reference_catalog['dec'])
    _, separations, _ = offset_skycoords.match_to_catalog_sky(ref_skycoords)
    # Minimize the mean offset between all of the matched sources
    return np.mean(separations.to(units.arcsec))


def coords_inside_image(input_wcs, catalog, origin=1):
    xs, ys = input_wcs.all_world2pix(catalog['ra'], catalog['dec'], origin)
    in_image = np.logical_and(xs > 0, ys > 0)
    in_image = np.logical_and(in_image, xs < input_wcs.naxis1)
    in_image = np.logical_and(in_image, ys < input_wcs.naxis2)
    return in_image


def read_catalog(catalog_filename):
    return ascii.read(catalog_filename, format='fast_no_header', names=['x', 'y', 'mag', 'magerr'])


def parse_dolphot_table(input_catalog, drzfile, output_catalog):
    dolphot_catalog = ascii.read(input_catalog, format='fast_no_header')

    dolphot_catalog = cut_bad_dolphot_sources(dolphot_catalog)
    t = table.Table()
    t['x'] = dolphot_catalog['col3'] + 0.5
    t['y'] = dolphot_catalog['col4'] + 0.5
    t['mag'] = dolphot_catalog['col17']
    t['magerr'] = dolphot_catalog['col19']
    t.writeto(output_catalog, format='fast_no_header', overwrite=True)
