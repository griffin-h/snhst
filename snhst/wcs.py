def read_phot_table(catname, drzfile, mag_column, output_catalog_name, write_radec=False):
    d = ascii.read(catname, format='fast_no_header')
    # Reject bad sources
    d = d[d['col11'] == 1]
    # Sharpness cut
    d = d[np.abs(d['col7']) < 0.3]
    # Crowding cut
    d = d[d['col10'] < 0.5]

    # Calculate the RAs and Decs for the remaining stars
    header = fits.getheader(drzfile)
    wcs = WCS(header)
    ras, decs = wcs.all_pix2world(d['col3'] + 0.5, d['col4'] + 0.5, 1)

    # Select only the stars within 20" of SN 2011fe
    coords = SkyCoord(ras, decs, unit=units.deg)
    offsets = coords.separation(
        SkyCoord('14:03:05.8', '+54:16:25', unit=(units.hourangle, units.deg)))
    d = d[offsets.arcsec <= 20.0]

    ras = ras[offsets.arcsec <= 20.0]
    decs = decs[offsets.arcsec <= 20.0]
    # Take the 120 brightest
    sorted_inds = d.argsort(mag_column)
    d.sort(mag_column)

    d = d[:120]

    if write_radec:
        d['col3'] = ras[sorted_inds][:120]
        d['col4'] = decs[sorted_inds][:120]
    else:
        d['col3'] += 0.5
        d['col4'] += 0.5
    # Write out the coordinates
    d = d['col3', 'col4', mag_column]
    d.write(output_catalog_name, overwrite=True)