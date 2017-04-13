import shutil
import os
from astropy.io import fits


def is_sci_extension(hdu):
    return 'EXTNAME' in hdu.header and hdu.header['EXTNAME'].upper() == 'SCI'


def copy_if_not_exists(filename, path):
    if not os.path.exists(os.path.join(path, filename)):
        shutil.copy(filename, path)


def get_instrument(filename):
    hdu = fits.open(filename, mode='readonly')
    instrument = hdu[0].header['INSTRUME'].lower()
    if instrument.upper == 'WFPC2':
        detector = 'wfpc2'
        subarray = 'full'
    else:
        detector = hdu[0].header['DETECTOR'].lower()
        if hdu[0].header['SUBARRAY'].upper() == 'T':
            subarray = "sub"
        else:
            subarray = "full"
    return "{instrument}_{detector}_{subarray}".format(instrument=instrument, detector=detector,
                                                       subarray=subarray)
