def run_tweakreg(options):
    if find_shifts:
        # This wcskey is created in the header of the flt files.
        # While this works fine on clean flt files, you can not rerun tweakreg
        # on flt files it will fail. Of course it is not just written
        # as a headerlet that is easy to delete.
        # It is written as several WCS keywords in the header.
        # We should add a step to look for the TWEAK wcskey and
        # remove them if they exist.

        wcskey = 'TWEAK'
        # Run tweakreg on the flt frames to get the relative shifts
        # between frames
        imfind_lines = [
            '_task_name_ = imagefindpars# \n',
            'computesig = True# Automatically compute sigma for all inputs?\n',
            'skysigma = 0.0# Standard deviation of background in counts\n',
            'conv_width = 3.5# Convolution kernel width in scale units\n',
            'peakmin = None# Min source peak value\n',
            'peakmax = None# Max source peak value\n',
            'threshold = %g# Threshold in sigma for feature detection\n' % threshold,
            # this seems to work well for crowded fields, still needs to be tested
            # in uncrowded fields
            'nsigma = 1.5# Width of convolution kernel in sigma\n',
            'fluxmin = None# Min good total source flux\n',
            'fluxmax = None# Max good total source flux\n',
            'sharplo = 0.2\n',
            'ratio = 1.0\n',
            'roundhi = 1.0\n',
            'roundlo = -1.0\n',
            'use_sharp_round = True\n',
            'dqbits = None\n',
            'theta = 0.0\n',
            'sharphi = 1.0\n',
            '\n']
        f = open('imfind.cfg', 'w')
        f.writelines(imfind_lines)
        f.close()
        # run la cosmic on the input images to make registering the images
        # easier if the instrument is not wfc3_ir
        for j, img in enumerate(imgs_full):
            if instrument != 'wfc3_ir':
                # Make a copy of the original data file to use later
                os.system('cp -f ' + img + ' ' + imgs[j] + '_raw.fits')

                hdulist = pyfits.open(img)
                for i in range(1, nchip + 1):
                    sci = hdulist['sci', i].data
                    _crmask, crclean = detect_cosmics(sci.copy().astype('<f4'), satlevel=satval,
                                                      gain=img_gain)
                    sci[:, :] = crclean[:, :]

                hdulist.writeto(img, clobber=True)
                hdulist.close()

        tools.teal.unlearn('tweakreg')
        tweaklines = [
            '_task_name_ = tweakreg# \n',
            'input = ' + input_files + '# Input files (name, suffix, or @list)\n',
            'refimage = "%s"# Filename of reference image\n' % template_image,
            'exclusions = ""# Filename for source exclusions catalogs\n',
            'updatewcs = False# Update WCS keywords of input images?\n',
            'writecat = True# Write out source catalogs?\n',
            'clean = True# Remove intermediate files?\n',
            'verbose = False# Print extra messages during processing?\n',
            'interactive = False\n'
            'runfile = tweakreg.log# Filename of processing log\n',
            '\n',
            '[UPDATE HEADER]\n',
            'updatehdr = True# Update headers of input files with shifts?\n',
            'wcsname = TWEAK# Name of updated WCS\n',
            'reusename = True\n'
            '\n',
            '[HEADERLET CREATION]\n',
            'headerlet = False# Create headerlet for solution?\n',
            'attach = True# Create headerlet FITS extension?\n',
            'hdrfile = ""# Filename for headerlet FITS file\n',
            'clobber = False# "Overwrite existing headerlet FITS file?"\n',
            'hdrname = ""# Unique name(HDRNAME) for headerlet\n',
            'author = ""# Author name for creator of headerlet\n',
            'descrip = ""# Short description of headerlet solution\n',
            'catalog = ""# Name of catalog used for headerlet solution\n',
            'history = ""# Name of ASCII file containing history for headerlet\n',
            '\n',
            '[OPTIONAL SHIFTFILE OUTPUT]\n',
            'shiftfile = False# Create output shiftfile?\n',
            'outshifts = shifts.txt# Filename of generated shiftfile\n',
            'outwcs = shifts_wcs.fits# Filename of shiftfile reference WCS\n',
            '\n',
            '[COORDINATE FILE DESCRIPTION]\n',
            'catfile = ""# File containing coordinate filenames for input files\n',
            'xcol = 1# "Column name(s) for X positions"\n',
            'ycol = 2# "Column name(s) for Y positions"\n',
            'fluxcol = ""# "Column name for source flux/mag values"\n',
            'maxflux = None# Maximum flux value for valid objects\n',
            'minflux = None# Minimum flux value for valid objects\n',
            'fluxunits = counts# Units of flux values\n',
            'xyunits = pixels# Units of X/Y positions\n',
            'nbright = None# Number of brightest objects to keep\n',
            '\n',
            '[REFERENCE CATALOG DESCRIPTION]\n',
            'refcat = ""# Filename of reference coordinate catalog\n',
            'refxcol = 1# "Column name(s) for RA"\n',
            'refycol = 2# "Column name(s) for Dec"\n',
            'refxyunits = degrees# Units of sky positions\n',
            'rfluxcol = ""# "Column name for source flux/mag values"\n',
            'rmaxflux = None# Maximum flux value for valid reference objects\n',
            'rminflux = None# Minimum flux value for valid reference objects\n',
            'rfluxunits = mag# Units of flux values\n',
            'refnbright = None# Number of brightest reference objects to keep\n',
            '\n',
            '[OBJECT MATCHING PARAMETERS]\n',
            'minobj = 15# Minimum number of objects acceptable for matching\n',
            'searchrad = 1.0# The search radius for a match\n',
            'searchunits = arcseconds# Units for search radius\n',
            'use2dhist = True# Use 2d histogram to find initial offset?\n',
            'see2dplot = True# See 2d histogram for initial offset?\n',
            # for crowded fields
            'separation = 0.5# Minimum object separation (pixels)\n',
            'tolerance = 1.0# Matching tolerance for xyxymatch(pixels)\n',
            'xoffset = 0.0# Initial guess for X offset(pixels)\n',
            'yoffset = 0.0# Initial guess for Y offset(pixels)\n',
            '\n',
            '[CATALOG FITTING PARAMETERS]\n',
            'fitgeometry = rscale# Fitting geometry\n',
            'residplot = "both"# Plot residuals from fit?\n',
            'nclip = 3 # Number of clipping iterations in fit\n',
            'sigma = 3.0# Clipping limit in sigma units\n',
            '\n']

        f = open('tweak.cfg', 'w')
        f.writelines(tweaklines)
        f.close()

        imagefindcfg = getDefaultConfigObj('imagefind', 'imfind.cfg')
        tweakreg.TweakReg(configobj='tweak.cfg', imagefindcfg=imagefindcfg)
        os.system('rm -rf tweak.cfg')
        os.system('rm -rf imfind.cfg')