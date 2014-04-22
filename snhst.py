'''
Created on Dec 17, 2010

@author: cmccully
'''
import os
from pyraf import iraf
from iraf import stsdas,dither,images,hst_calib,wfc3
import pyfits,acstools,numpy,lacosmicx
import drizzlepac
from drizzlepac import tweakreg, astrodrizzle
from drizzlepac.util import getDefaultConfigObj
from astLib.astWCS import WCS
from stsci import tools
from numpy import cos,sin,array,pi,arctan
from numpy.linalg import norm
import readcol
from pyfits import getval
import resource
from glob import glob
resource.setrlimit(resource.RLIMIT_NOFILE, (10000,-1))

def calccoarsealign(x1s, x2s,s_img, x1t, x2t,t_img):
    #Takes the position of two stars on two different images and calculates the values
    #needed for the shift file. s = source, t = template
    twcs = WCS(t_img)
    swcs = WCS(s_img)
    
    x1s = array(swcs.pix2wcs(x1s[0],x1s[1]))
    x2s  =array(swcs.pix2wcs(x2s[0],x2s[1]))
    x1t = array(twcs.pix2wcs(x1t[0],x1t[1]))
    x2t = array(twcs.pix2wcs(x2t[0],x2t[1]))
    crval_s = array([float(getval(s_img,'CRVAL1')),float(getval(s_img,'CRVAL2'))])
    crval_t = array([float(getval(t_img,'CRVAL1')),float(getval(t_img,'CRVAL2'))])
    
    x1s -= crval_s
    x2s -= crval_s
    x1t -= crval_t
    x2t -= crval_t
    
    dxt = x2t - x1t
    dxs = x2s - x1s
    
    a = norm(dxt)/norm(dxs)
    
    dxtdyt = dxt[0]/dxt[1]
    th = arctan((dxtdyt*dxs[1] - dxs[0])/(dxs[1] + dxtdyt * dxs[0]))
    
    x0 = x1t[0] - a*(cos(th)* x1s[0] + sin(th)*x1s[1])
    y0 = x1t[1] - a*(-sin(th)*x1s[0] + cos(th)*x1s[1])

    #This follows the same convention as the overall wcs params in snhst.drizzle
    return x0,y0,th*180.0/pi,a

def coarsealign(imgs,nchip,drizra,drizdec,rashift,decshift,rot,scale):
    for img in imgs:
        hdulist = pyfits.open(img)
        for i in range(1,nchip*3+1):
            cd1_1 = float(hdulist[i].header['CD1_1'])
            cd1_2 = float(hdulist[i].header['CD1_2'])
            cd2_1 = float(hdulist[i].header['CD2_1'])
            cd2_2 = float(hdulist[i].header['CD2_2'])
            
            cd = array([[cd1_1,cd1_2],[cd2_1,cd2_2]])
            rotmat = array([[scale * cos(rot/180.0 * pi), scale*sin(rot/180.0 * pi)],
                           [-scale*sin(rot/180.0 * pi),scale*cos(rot/180.0 * pi)]])
            
            cdout = rotmat.dot(cd)
            
            hdulist[i].header['CD1_1']= cdout[0][0]
            hdulist[i].header['CD1_2']= cdout[0][1]
            hdulist[i].header['CD2_1']= cdout[1][0]
            hdulist[i].header['CD2_2']= cdout[1][1]
            
            crval1 = float(hdulist[i].header['CRVAL1'])
            crval2 = float(hdulist[i].header['CRVAL2'])
            
            #Translate to the center of the drizzle frame
            rotshift = array([ crval1 - drizra , crval2 - drizdec])
            #rotate and translate back
            rotshift -= rotmat.dot(rotshift) 
            hdulist[i].header['CRVAL1'] = crval1 + rashift - rotshift[0]
            hdulist[i].header['CRVAL2'] = crval2 + decshift - rotshift[1]
        hdulist.writeto(img,clobber=True)
        hdulist.close()
        
def drizzle(output_name,input_files='',ref='',template_image='',
            instrument='wfc3_ir',drizra=0.0,drizdec=0.0,pix_scale=0.0,drizrot=0.0,nx=0,ny=0,
            pix_frac=1.0,acs_cte=False,do_destripe=True,
            clean = True, find_shifts = True,scale = 1.0, rot = 0.0,rashift=0.0,decshift =0.0, threshold= 30.0,num_cores = 8):
    
    this_dir=os.getcwd()+'/'
    if ref=='': ref=this_dir
    #use the multdrizzle parameters for the correct instrument
    #pars=[group,proc_unit,pix_scale,pix_frac,sep_bits,final_bits]
    #nchip = number of chips on the instrument (used to get the right flt or c0m extensions)
    #Do use smaller pixels than native to help differentiate between cosmic rays.
    #Do leave in saturated pixels for cosmics
    if instrument == 'wfc3_uvis':
        [grp,units,sep_bits,final_bits]=['','electrons',0,0]
        os.environ['iref']=ref
        nchip=2
        
        #parameters for la cosmic
        rdnoise=6.5
        img_gain=1.0
        satval=70000.0
        sig_clip=4.0
        sig_frac=0.2
        obj_lim=6.0
        if template_image=='' and nx == 0:
            nx=4850
            ny=4550
            # native_pix_scale=0.04
        if template_image =='' and pix_scale == 0.0: pix_scale=0.04
        if input_files=='': input_files='*_flt.fits'
   
    elif instrument == 'acs':
        [grp,units,sep_bits,final_bits]=['','electrons',0,0]
        os.environ['jref']=ref
        nchip=2

        #parameters for la cosmic
        rdnoise=6.5
        img_gain=1.0
        satval=70000.0
        sig_clip=3.0
        sig_frac=0.1
        obj_lim=5.0
        
        if template_image=='' and nx==0:
            nx=5200
            ny=5200
            # native_pix_scale=0.0495
        if template_image =='' and pix_scale == 0.0: pix_scale=0.05
        if input_files=='': input_files='*_flc.fits'
    
    elif instrument == 'acs_hrc':
        [grp,units,sep_bits,final_bits]=['','electrons',256,256]
        os.environ['jref']=ref
        nchip=1
        
        #parameters for la cosmic
        rdnoise=6.5
        img_gain=1.0
        satval=155000.0
        sig_clip=3.0
        sig_frac=0.1
        obj_lim=5.0
        
        if template_image=='' and nx==0:
            nx=1400
            ny=1400
            # native_pix_scale=0.0495
        if template_image =='' and pix_scale == 0.0: pix_scale=0.025
        if input_files=='': input_files='*_flt.fits'
            
    elif instrument=='wfc3_ir' :
        [grp,units,sep_bits,final_bits]=['','electrons',512+256,512+256]
        os.environ['iref']=ref
        nchip=1
        
        #La Cosmic is not used for WFC3 IR
        if template_image=='' and nx==0:
            nx=1700
            ny=1700
            #native_pix_scale=0.1282
        
        if template_image =='' and pix_scale == 0.0: pix_scale=0.09
        if input_files=='': input_files='*_flt.fits'  
            
    elif instrument == 'stis':
        [grp,units,sep_bits,final_bits]=['','electrons',256,256]
        os.environ['oref']=ref
        nchip=1
        
        #parameters for la cosmic
        rdnoise=6.5
        img_gain=1.0
        satval=30000.0
        sig_clip=3.0
        sig_frac=0.1
        obj_lim=5.0
        if template_image=='' and nx == 0:
            nx=1400
            ny=1400
            # native_pix_scale=0.04
        if template_image =='' and pix_scale == 0.0: pix_scale=0.0507
        if input_files=='':input_files='*_flt.fits'

    elif instrument=='wfpc2_pc' :    
        [grp,units,sep_bits,final_bits]=['1','electrons',8,8]
        os.environ['uref']=ref
        iraf.set(uref=ref)
        nchip=1

        #la cosmic parameters: these have not been optimized for WFPC2, they are just best guesses
        rdnoise=10.0
        img_gain=7.0
        satval=27000.0
        sig_clip=4.0
        sig_frac=0.3
        obj_lim=6.0
        
        if template_image=='' and nx==0:
            nx=1000
            ny=1000
            #native_pix_scale=0.046
        if template_image =='' and pix_scale == 0.0: pix_scale=0.046    
        if input_files=='':input_files='*_c0m.fits'
        
    elif instrument=='wfpc2_all' :    
        [grp,units,sep_bits,final_bits]=['','electrons',8,8]
        os.environ['uref']=ref
        iraf.set(uref=ref)
        nchip=1
        
        #la cosmic parameters: these have not been optimized for WFPC2, they are just best guesses
        rdnoise=10.0
        img_gain=7.0
        satval=27000.0
        sig_clip=4.0
        sig_frac=0.3
        obj_lim=6.0
        
        #you must manually set the pix scale, and nx and ny, or have a wcs_template image
        if input_files=='': input_files='*_c0m.fits'  
           
    elif instrument=='wfpc2_wf' :
        [grp,units,sep_bits,final_bits]=['2,3,4','electrons',8,8]
        os.environ['uref']=ref
        iraf.set(uref=ref)
        nchip=4
        
        #la cosmic parameters: these have not been optimized for WFPC2, they are just best guesses
        rdnoise=10.0
        img_gain=7.0
        satval=27000.0
        sig_clip=4.0
        sig_frac=0.3
        obj_lim=6.0
        
        if template_image=='' and nx==0:
            nx=2400
            ny=2400
            #native_pix_scale=0.0996
        if template_image =='' and pix_scale == 0.0: pix_scale=0.0996
            
        if input_files=='': input_files='*_c0m.fits'  
        
    #convert the input list for multidrizzle into a useable list of images
    imgs = []
    imgs_full=tools.parseinput.parseinput(input_files)[0]
    
    for img in imgs_full:
        #strip off the the extension (_flt or _c0m) and the .fits
        this_img= img[:-9]
        if instrument in ['wfpc2_pc' ,'wfpc2_wf' ,'wfpc2_all']:imgs.append(this_img+'_c0m')
        else:imgs.append(this_img)

    #the ra and the dec are the desired ra and dec for the center of the frame
    if drizra==0.0 and drizdec == 0.0:
        #grab the target ra and dec from the header of the first file
        hdulist=pyfits.open(imgs_full[0])
        #find the midpoint of the cr vals for the different chips
        if instrument in ['wfc3_ir','acs_hrc']:
            hdr=hdulist['SCI'].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra,drizdec = imwcs.getCentreWCSCoords() 

        elif instrument in ['acs','wfc3_uvis']:
            hdr=hdulist['SCI',1].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra,drizdec = imwcs.getCentreWCSCoords() 

            hdr=hdulist['SCI',2].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra +=imwcs.getCentreWCSCoords()[0] 
            drizdec += imwcs.getCentreWCSCoords()[1]
            drizra /= 2.0
            drizdec /= 2.0
        elif instrument=='wfpc2_wf':
            hdr=hdulist['SCI',2].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra,drizdec = imwcs.getCentreWCSCoords() 
            
            hdr=hdulist['SCI',3].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra +=imwcs.getCentreWCSCoords()[0] 
            drizdec += imwcs.getCentreWCSCoords()[1]
            
            hdr=hdulist['SCI',4].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra +=imwcs.getCentreWCSCoords()[0] 
            drizdec += imwcs.getCentreWCSCoords()[1]
            drizra /=   3.0
            drizdec /= 3.0
        elif instrument=='wfpc2_pc':
            hdr=hdulist['SCI',1].header
            imwcs = WCS(hdr,mode='pyfits')
            drizra,drizdec = imwcs.getCentreWCSCoords() 

        hdulist.close()
        
    if template_image!='':
        #Copy the template file here as to not run into file collisions
        os.system('cp -f '+template_image+' ./template.fits')
        template_image = 'template.fits'
        hdulist=pyfits.open('template.fits')
        #Get the WCS parameters for the template image
        twcs = WCS('template.fits')
        hdr=hdulist[0].header
        nx =int(hdr['naxis1'])
        ny =int(hdr['naxis2'])
        
        drizra,drizdec = twcs.getCentreWCSCoords()
        
        pix_scale=3600.0*twcs.getPixelSizeDeg()
        #Add a half pixel to conform to the astrodrizzle convention
        drizra += 0.5*pix_scale/3600.0
        drizdec -= 0.5*pix_scale/3600.0
        drizrot = twcs.getRotationDeg()
        #make a template weight image
        weight_mask=hdulist[0].data.copy()
        weight_hdr=hdr.copy()
        hdulist.close()
        
        #Make a template weight mask for sextractor etc.
        weight_mask[weight_mask != 0.0] = 1.0
        new_hdu=pyfits.PrimaryHDU(weight_mask,weight_hdr)
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto('template_weight.fits',  clobber=True)
        new_hdulist.close()
        
    #if ACS data do a cte correction a destripe the image
    #We no longer do this by default because we default to using the flc files which are already corrected.
    if instrument=='acs' and acs_cte:
        for i,img in enumerate(imgs_full):
            #turn cte into a function use memory properly
            run_cte(img,imgs[i],do_destripe = do_destripe)
    
    wcskey = ''
    if find_shifts:
        #This wcskey is created in the header of the flt files.
        #While this works fine on clean flt files, you can not rerun tweakreg on flt files
        #It will fail. Of course it is not just written as a headerlet that is easy to delete.
        #It is written as several WCS keywords in the header. 
        #We should add a step to look for the TWEAK wcskey and remove them if they exist.
        
        wcskey='TWEAK'
        #Run tweakreg on the flt frames to get the relative shifts between frames
        imfind_lines = ['_task_name_ = imagefindpars# \n',
                        'computesig = True# Automatically compute sigma for all inputs?\n',
                        'skysigma = 0.0# Standard deviation of background in counts\n',
                        'conv_width = 3.5# Convolution kernel width in scale units\n',
                        'peakmin = None# Min source peak value\n',
                        'peakmax = None# Max source peak value\n',
                        'threshold = %g# Threshold in sigma for feature detection\n'% threshold,# this seems to work well for crowded fields, still needs to be tested in uncrowded fields                  
                        'nsigma = 1.5# Width of convolution kernel in sigma\n',
                        'fluxmin = None# Min good total source flux\n',
                        'fluxmax = None# Max good total source flux\n',
                        '\n']
        f = open('imfind.cfg','w')
        f.writelines(imfind_lines)
        f.close()
        #run la cosmic on the input images to make registering the images easier if the instrument is not wfc3_ir
        for j,img in enumerate(imgs_full):
            if instrument!='wfc3_ir':
                #Make a copy of the original data file to use later
                os.system('cp -f '+img+' '+imgs[j]+'_raw.fits')
            
                hdulist=pyfits.open(img)
                for i in range(1,nchip+1):
                    sci=hdulist['sci',i].data               
                    crclean=lacosmicx.run(sci.copy(),satlevel=satval,gain = img_gain)
                    sci[:,:]=crclean[:,:]
            
                hdulist.writeto(img,clobber=True)
                hdulist.close()

        tools.teal.unlearn('tweakreg')
        tweaklines = [  '_task_name_ = tweakreg# \n',
                                'input = '+input_files+'# Input files (name, suffix, or @list)\n',
                                'refimage = "%s"# Filename of reference image\n'%template_image,
                                'exclusions = ""# Filename for source exclusions catalogs\n',
                                'updatewcs = False# Update WCS keywords of input images?\n',
                                'writecat = False# Write out source catalogs?\n',
                                'clean = True# Remove intermediate files?\n',
                                'verbose = False# Print extra messages during processing?\n',
                                'runfile = tweakreg.log# Filename of processing log\n',
                                '\n',
                                '[UPDATE HEADER]\n',
                                'updatehdr = True# Update headers of input files with shifts?\n',
                                'wcsname = TWEAK# Name of updated WCS\n',
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
                                'fluxmax = None# Maximum flux value for valid objects\n',
                                'fluxmin = None# Minimum flux value for valid objects\n',
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
                                'rfluxmax = None# Maximum flux value for valid reference objects\n',
                                'rfluxmin = None# Minimum flux value for valid reference objects\n',
                                'rfluxunits = mag# Units of flux values\n',
                                'refnbright = None# Number of brightest reference objects to keep\n',
                                '\n',
                                '[OBJECT MATCHING PARAMETERS]\n',
                                'minobj = 15# Minimum number of objects acceptable for matching\n',
                                'searchrad = 1.0# The search radius for a match\n',
                                'searchunits = arcseconds# Units for search radius\n',
                                'use2dhist = True# Use 2d histogram to find initial offset?\n',
                                'see2dplot = False# See 2d histogram for initial offset?\n',
                                'separation = 2.0# Minimum object separation (pixels)\n'#for crowded fields
                                'tolerance = 1.0# Matching tolerance for xyxymatch(pixels)\n',
                                'xoffset = 0.0# Initial guess for X offset(pixels)\n',
                                'yoffset = 0.0# Initial guess for Y offset(pixels)\n',
                                '\n',
                                '[CATALOG FITTING PARAMETERS]\n',
                                'fitgeometry = rscale# Fitting geometry\n',
                                'residplot = "No plot"# Plot residuals from fit?\n',
                                'nclip = 20# Number of clipping iterations in fit\n',
                                'sigma = 2.0# Clipping limit in sigma units\n',
                                '\n'
                                ]
        f=open('tweak.cfg','w')
        f.writelines(tweaklines)
        f.close()

        tweakreg.TweakReg( configobj='tweak.cfg', imagefindcfg=getDefaultConfigObj('imagefind','imfind.cfg'))
        os.system('rm -rf tweak.cfg')
        os.system('rm -rf imfind.cfg')
        #Remove the cosmic ray cleaned (and which could have flux issues) images and replace them with the raw files from above
        for j,img in enumerate(imgs_full):
            if instrument!='wfc3_ir':
                hdulist=pyfits.open(img)
                hduraw = pyfits.open(imgs[j]+'_raw.fits')
                for i in range(1,nchip+1):
                    hdulist['sci',i].data[:,:]= hduraw['sci',i].data[:,:]
                hdulist.writeto(img,clobber=True)
                hdulist.close()
                hduraw.close()
                os.system('rm -rf '+img+'_raw.fits')
    
    # Apply overall shifts, rotation and scale to wcs
    if rashift!= 0.0 and decshift != 0.0:coarsealign(imgs_full,nchip,drizra,drizdec,rashift,decshift,rot,scale)
            
    if instrument == 'wfc3_ir': do_drizcr = False
    else: do_drizcr = True

    #run astrodrizzle
    tools.teal.unlearn('astrodrizzle')
    
    num_images=len(glob(input_files))
    if num_images ==1: medstr='mean'
    else: medstr='minmed'
    adlines = ['_task_name_ = astrodrizzle# \n',
                'input = ' +input_files+'# Input files (name, suffix, or @list)\n',
                'output = "'+output_name+'"# Rootname for output drizzled products\n',
                'runfile = ""# File for logging the processing\n', #Don't save a log for now.
                'updatewcs = False# Update the WCS keywords?\n', #This reverts things so it messes up the work tweakreg did if true
                'wcskey = "%s"# WCS version to use in processing\n'%wcskey,
                'proc_unit = '+units+'# Units used during processing\n',
                'coeffs = True# Use header-based distortion coefficients?\n',
                'context = False# Create context image during final drizzle?\n',#Don't create the context right now because we don't use it
                'group = "%s"# Single extension or group to be combined/cleaned\n'%grp,
                'build = False# Create multi-extension output file for final drizzle?\n',
                'crbit = 4096# Bit value for CR ident. in DQ array\n',
                'stepsize = 10# Step size for drizzle coordinate computation\n',
                'resetbits = 4096# Bit values to reset in all input DQ arrays\n' #reset all of the 4096 flags from multidrizzle for the HST pipeline,
                'num_cores = %i# Max CPU cores to use (n<2 disables, None = auto-decide)\n'%num_cores,
                'in_memory = False# Process everything in memory to minimize disk I/O?\n',#In principle this would be nice, but it doesn't update the input flt files with the new DQ mask if set to True
                '\n',
                '[STATE OF INPUT FILES]\n',
                'restore = False# Copy input files FROM archive directory for processing?\n',
                'preserve = False# Copy input files to archive directory, if not already archived?\n', #In principle it might be good to preserve the raw data but we do that manually
                'overwrite = False# Copy input files into archive, overwriting if required?\n',
                'clean = '+str(clean)+'# Delete temporary files after completion?\n',
                '\n',
                '[STEP 1: STATIC MASK]\n',
                'static = True# Create static bad-pixel mask from the data?\n',
                'static_sig = 4.0# "Sigma*rms below mode to clip for static mask"\n',
                '\n',
                '[STEP 2: SKY SUBTRACTION]\n',
                'skysub = True# "Perform sky subtraction?"\n',
                'skywidth = 0.1# "Bin width for sampling sky statistics (in sigma)"\n',
                'skystat = mode# "Sky correction statistics parameter"\n', #We use the mode here instead of the median and use a narrower width
                'skylower = 0.0# "Lower limit of usable data for sky (always in electrons)"\n',
                'skyupper = None# "Upper limit of usable data for sky (always in electrons)"\n',
                'skyclip = 5# "Number of clipping iterations"\n',
                'skylsigma = 4.0# Lower side clipping factor (in sigma)\n',
                'skyusigma = 4.0# Upper side clipping factor (in sigma)\n',
                'skyuser = ""# KEYWORD indicating a sky subtraction value if done by user.\n',
                'skyfile = ""# Name of file with user-computed sky values.\n',
                '\n',
                '[STEP 3: DRIZZLE SEPARATE IMAGES]\n',
                'driz_separate = True# "Drizzle onto separate output images?"\n',
                'driz_sep_kernel = turbo# Shape of kernel function\n',
                'driz_sep_wt_scl = exptime# "Weighting factor for input data image"\n',
                'driz_sep_pixfrac = 1.0# Linear size of drop in input pixels\n',
                'driz_sep_fillval = -100000 # Value to be assigned to undefined output points\n',# arbitrary low value; easy to exclude
                'driz_sep_bits = '+str(sep_bits)+'# Integer mask bit values considered good\n',
                '\n',
                '[STEP 3a: CUSTOM WCS FOR SEPARATE OUTPUTS]\n',
                'driz_sep_wcs = True# "Define custom WCS for separate output images?"\n',
                'driz_sep_refimage = "%s"# Reference image from which to obtain a WCS\n'%template_image,
                'driz_sep_rot = '+str(drizrot)+'# "Position Angle of drizzled image\'s Y-axis w.r.t. North (degrees)"\n',
                'driz_sep_scale = '+str(pix_scale)+'# Absolute size of output pixels in arcsec/pixel\n',
                'driz_sep_outnx = '+str(nx)+'# Size of separate output frame\'s X-axis (pixels)\n',
                'driz_sep_outny = '+str(ny)+'# Size of separate output frame\'s Y-axis (pixels)\n',
                'driz_sep_ra = '+str(drizra)+'# right ascension output frame center in decimal degrees\n',
                'driz_sep_dec = '+str(drizdec)+'# declination output frame center in decimal degrees\n',
                '\n',
                '[STEP 4: CREATE MEDIAN IMAGE]\n',
                'median = True# "Create a median image?"\n',
                'median_newmasks = True# "Create new masks when doing the median?"\n',
                'combine_maskpt = 0.2# "Percentage of weight image value below which it is flagged as a bad pixel."\n',
                'combine_type = %s# "Type of combine operation"\n'%medstr, #Choose minmed for now, this is better if we only have 2 flts
                'combine_nsigma = 4 3# "Significance for accepting minimum instead of median"\n',
                'combine_nlow = 0# "minmax: Number of low pixels to reject"\n',
                'combine_nhigh = 0# "minmax: Number of high pixels to reject"\n',
                'combine_lthresh = -10000# Lower threshold for clipping input pixel values\n',
                'combine_hthresh = None# "Upper threshold for clipping input pixel values"\n',
                'combine_grow = 1# Radius (pixels) for neighbor rejection\n',
                '\n',
                '[STEP 5: BLOT BACK THE MEDIAN IMAGE]\n',
                'blot = True# "Blot the median back to the input frame?"\n',
                'blot_interp = poly5# Interpolant (nearest,linear,poly3,poly5,sinc)\n',
                'blot_sinscl = 1.0# Scale for sinc interpolation kernel\n',
                'blot_addsky = True# "Add sky using MDRIZSKY value from header?"\n', # I think we don't want a sky here so that driz_cr works correctly
                'blot_skyval = 0.0# Custom sky value to be added to blot image\n',
                '\n',
                '["STEP 6: REMOVE COSMIC RAYS WITH DERIV, DRIZ_CR"]\n',
                'driz_cr = '+str(do_drizcr)+'# Perform CR rejection with deriv and driz_cr?\n',
                'driz_cr_corr = False# "Create CR cleaned _crclean file and a _crmask file?"\n',
                'driz_cr_snr = 3.5 3.0# "Driz_cr.SNR parameter"\n',
                'driz_cr_grow = 1# Driz_cr_grow parameter\n',
                'driz_cr_ctegrow = 0# Driz_cr_ctegrow parameter\n',
                'driz_cr_scale = 1.2 0.7# Driz_cr.scale parameter\n',
                '\n',
                '[STEP 7: DRIZZLE FINAL COMBINED IMAGE]\n',
                'driz_combine = True# "Perform final drizzle image combination?"\n',
                'final_wht_type = EXP# Type of weighting for final drizzle\n',
                'final_kernel = square# Shape of kernel function\n',
                'final_wt_scl = exptime# Weighting factor for input data image\n',
                'final_pixfrac = '+str(pix_frac)+'# Linear size of drop in input pixels\n',
                'final_fillval = -50000# "Value to be assigned to undefined output points" flag so we don\'t add the sky back to the bad pixels\n',
                'final_bits = '+str(final_bits)+'# Integer mask bit values considered good\n',
                'final_units = counts # Units for final drizzle image (counts or cps)\n',
                '\n',
                '[STEP 7a: CUSTOM WCS FOR FINAL OUTPUT]\n',
                'final_wcs = True# "Define custom WCS for final output image?"\n',
                'final_refimage = "%s"# Reference image from which to obtain a WCS\n'%template_image,
                'final_rot = '+str(drizrot)+'# "Position Angle of drizzled image\'s Y-axis w.r.t. North (degrees)" \n',
                'final_scale = '+str(pix_scale)+' # Absolute size of output pixels in arcsec/pixel\n',
                'final_outnx = '+str(nx)+'# Size of FINAL output frame X-axis (pixels)\n',
                'final_outny = '+str(ny)+'# Size of FINAL output frame Y-axis (pixels)\n',
                'final_ra = '+str(drizra)+'# right ascension output frame center in decimal degrees\n',
                'final_dec = '+str(drizdec)+'# declination output frame center in decimal degrees\n',
                '\n',
                '[INSTRUMENT PARAMETERS]\n',
                'gain = ""# \n',
                'gnkeyword = ""# \n',
                'rdnoise = ""# \n',
                'rnkeyword = ""# \n',
                'exptime = ""# \n',
                'expkeyword = ""# \n',
                '\n']

    f = open('adriz.cfg','w')
    f.writelines(adlines)
    f.close()
    astrodrizzle.AstroDrizzle( configobj = 'adriz.cfg') 
                 
    
    if instrument =='acs' : drzstr = '_drc'
    else: drzstr = '_drz'    
    if 'flc' in input_files: drzstr='_drc'
    hdulist=pyfits.open(output_name+drzstr+'_sci.fits')
    hdr=hdulist[0].header
    cosmic_hdr=hdr.copy()
    #add the sky value back in 
    #only add the sky to non flagged pixels (-50k)
    sci=hdulist[0].data
    mdrizsky=numpy.min(sci[sci > -49999.0])
    sci[sci > -49999.0 ]-=mdrizsky
    no_data=sci < -49999.0
    sci[no_data]=0.0
    
    sci_cosmic=sci.copy()
    hdulist.writeto(output_name+drzstr+'_sci.fits',clobber=True)
    hdulist.close()
    
    if instrument!='wfc3_ir':
        
        #mask out the missing data
        clean=lacosmicx.run(sci_cosmic,inmask= sci_cosmic==0,outmaskfile=output_name+'_lamask.fits',gain=img_gain, readnoise=rdnoise, sigclip = sig_clip, sigfrac = sig_frac, objlim = obj_lim,satlevel=satval)
        
        new_hdu=pyfits.PrimaryHDU(clean,cosmic_hdr)
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto(output_name+'_lacosmic.fits',  clobber=True)
        new_hdulist.close()
    
    #data products are now finished
    
    #clean up
    os.system('rm -rf adriz.cfg')
    os.system('rm -rf template_weight.fits')
    os.system('rm -rf template.fits')
    
def subtract(image,template_image,output_image,instrument='wfc3_ir',template_instrument='wfc3_ir',imarith=False):    
    #use hotpants with the correct parameters based on the instrument
    if instrument=='wfc3_ir':
        image_gain=1.0
        iuthresh= 75000.0
        iukthresh= 50000.0
    elif instrument=='acs':
        image_gain=1.0
        iuthresh= 70000.0
        iukthresh= 60000.0
    elif instrument=='wfc3_uvis':
        image_gain=1.0
        iuthresh= 70000.0
        iukthresh= 60000.0
    elif instrument=='wfpc2_pc':
        image_gain=7.0
        iuthresh= 27000.0
        iukthresh= 27000.0
    elif instrument=='wfpc2_wf':
        image_gain=7.0
        iuthresh= 27000.0
        iukthresh= 27000.0
    elif instrument=='acs_hrc':
        image_gain=1.0
        iuthresh= 100000.0
        iukthresh= 100000.0
    elif instrument=='wfpc2_all':
        image_gain=7.0
        iuthresh= 10000.0
        iukthresh= 10000.0

        
    if template_instrument=='wfc3_ir':
        template_gain=1.0
        tuthresh=75000.0
        tukthresh=50000.0
    elif template_instrument=='acs':
        template_gain=1.0
        tuthresh=75000.0
        tukthresh=50000.0
    elif template_instrument=='wfc3_uvis':
        template_gain=1.0
        tuthresh=75000.0
        tukthresh=50000.0
    elif template_instrument=='wfpc2_pc':
        template_gain=7.0
        tuthresh=27000.0
        tukthresh=27000.0
    elif template_instrument=='wfpc2_wf':
        template_gain=7.0
        tuthresh=27000.0
        tukthresh=27000.0
    elif template_instrument=='acs_hrc':
        template_gain=1.0
        tuthresh=100000.0
        tukthresh=100000.0
    elif template_instrument=='wfpc2_all':
        template_gain=7.0
        tuthresh=10000.0
        tukthresh=10000.0
    
    os.system('cp -f '+template_image+' template.fits')
    
    if imarith:
        #do an imarith like subtraction
        hdulist=pyfits.open(image)
        im_data=hdulist[0].data.copy()
        im_hdr=hdulist[0].header.copy()
        exptime_im=float(im_hdr['EXPTIME'])
        hdulist.close()
        
        #get the template data
        hdulist=pyfits.open('template.fits')
        temp_data=hdulist[0].data.copy()
        temp_hdr=hdulist[0].header.copy()
        exptime_temp=float(temp_hdr['EXPTIME'])
        hdulist.close()
        
        sub_data=numpy.zeros(im_data.shape)
        #find where there is data for both template and image
        overlap=(im_data != 0.0) & (temp_data != 0.0)  
        sub_data[overlap]=im_data[overlap]-numpy.median(im_data[overlap]) - (temp_data[overlap]-numpy.median(temp_data[overlap])) *exptime_im/exptime_temp 
         
        new_hdu=pyfits.PrimaryHDU(sub_data,im_hdr)
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto(output_image+'_imarithsub.fits',  clobber=True) 
        new_hdulist.close()  
        
        new_hdu_inv=pyfits.PrimaryHDU(sub_data.copy()*-1.0,im_hdr)
        new_hdulist_inv=pyfits.HDUList([new_hdu_inv])
        new_hdulist_inv.writeto(output_image+'_imarithsub_inv.fits',  clobber=True) 
        new_hdulist.close()  

    else:
        #do the hotpants subtraction
        #make mask images
        
        hdulist=pyfits.open(image)
        im_data=hdulist[0].data.copy()
        im_hdr=hdulist[0].header.copy()
        hdulist.close()
        im_mask=numpy.zeros(im_data.shape,dtype='uint16')
        im_mask[im_data==0.0]=1
        im_hdr['BITPIX']=32
        new_hdu=pyfits.PrimaryHDU(im_mask,im_hdr)
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto("im_mask.fits",  clobber=True)
        new_hdulist.close()  
        
        hdulist=pyfits.open("template.fits")
        temp_hdr=hdulist[0].header.copy()
        temp_data=hdulist[0].data.copy()
        hdulist.close()
        temp_mask=numpy.zeros(temp_data.shape,dtype='uint16')
        temp_mask[temp_data==0.0]=1
        temp_hdr['BITPIX']=32
        new_hdu=pyfits.PrimaryHDU(temp_mask,temp_hdr)
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto("temp_mask.fits",  clobber=True)
        new_hdulist.close()  
  
        if instrument=='wfc3_ir':
            cmd='hotpants -inim '+image+' -tmplim template.fits -outim '+output_image+'_sub.fits -tmi temp_mask.fits -imi im_mask.fits -il 1e-29 -tl 1e-29 -n i -ig '+str(image_gain)+' -tg '+str(template_gain)+' -iu '+str(iuthresh)+' -tu '+str(tuthresh)+ ' -iuk '+str(iukthresh)+' -tuk '+str(tukthresh)+' -nsx 25 -nsy 25 -ng 3 8 0.7 6 1.5 4 3.0'
            os.system(cmd)
        else:
            cmd='hotpants -inim '+image+' -tmplim template.fits -outim '+output_image+'_sub.fits -tmi temp_mask.fits -imi im_mask.fits -il 1e-29 -tl 1e-29 -n i -ig '+str(image_gain)+' -tg '+str(template_gain)+' -iu '+str(iuthresh)+' -tu '+str(tuthresh)+ ' -iuk '+str(iukthresh)+' -tuk '+str(tukthresh)+' -nsx 25 -nsy 25'
            os.system(cmd)
            
        hdulist=pyfits.open(output_image+'_sub.fits')
        inv_hdr=hdulist[0].header.copy()
        inv_data=hdulist[0].data.copy()*-1.0
        new_hdu_inv=pyfits.PrimaryHDU(inv_data,inv_hdr)
        new_hdulist_inv=pyfits.HDUList([new_hdu_inv])
        new_hdulist_inv.writeto(output_image+'_sub_inv.fits',  clobber=True) 
        new_hdulist.close()
        #os.system('rm -rf temp_mask.fits')
        #os.system('rm -rf im_mask.fits')
    os.system('rm -rf template.fits')
    
def sextractor(image,thresh=2.5,minarea=4):
    #run sexextractor to find any sources in the image using the correct options
    hdulist=pyfits.open(image)
    weight_mask=hdulist[0].data.copy()
    weight_hdr=hdulist[0].header.copy()
    hdulist.close()
    
    weight_mask[abs(weight_mask) > 1e-10] = 1.0
    new_hdu=pyfits.PrimaryHDU(weight_mask,weight_hdr)
    new_hdulist=pyfits.HDUList([new_hdu])
    new_hdulist.writeto('weight.fits',  clobber=True)
    new_hdulist.close()
    
    #write the addsex file
    lines=[]

    lines.append("# additional SExtractor configuration parameters\n")
    lines.append("WEIGHT_TYPE     MAP_WEIGHT\n")
    lines.append("WEIGHT_IMAGE     weight.fits\n")
    lines.append("WEIGHT_THRESH    0.1\n")
    lines.append("INTERP_TYPE     ALL\n")
    lines.append("WEIGHT_GAIN	y\n")
    
    file = open('addsex.conf', 'w')
    file.writelines(lines)
    file.close()

    os.system('runsex_cm.pl '+image+' '+ str(thresh)+' -minarea '+str(minarea))
    os.system('rm -rf weight.fits')
    x=readcol.readcol(image+'.stars','whitespace','#',1)
    y=readcol.readcol(image+'.stars','whitespace','#',2)
    lines=[]
    for i,thisx in enumerate(x):
        lines.append('circle '+thisx+' '+y[i]+' 10.0 \n')
    f = open(image+'.reg', 'w')
    f.writelines(lines)
    f.close()


def run(output_name,input_files='',ref='',sub_template='',phot_sub_template='',image_instrument='wfc3_ir',templ_instrument='wfc3_ir',wcs_template='',do_cte=False,ra=0.0,dec=0.0,pixel_scale=0.0,this_nx=0,this_ny=0,do_drizzle=True,this_pix_frac=1.0):
  
    if do_drizzle:
        drizzle(output_name, input_files,ref, template_image=wcs_template, instrument=image_instrument,acs_cte=do_cte,drizra=ra,drizdec=dec, pix_scale=pixel_scale,nx=this_nx,ny=this_ny,pix_frac=this_pix_frac)

    if sub_template != '':
        
        if image_instrument !='wfc3_ir':
            subtract(image=output_name+'_lacosmic.fits', template_image=sub_template, output_image=output_name+'_lacosmic', instrument=image_instrument, template_instrument=templ_instrument)
            subtract(image=output_name+'_lacosmic.fits', template_image=sub_template, output_image=output_name+'_lacosmic', instrument=image_instrument, template_instrument=templ_instrument,imarith=True)
            sextractor(image=output_name+'_lacosmic_sub.fits')
            sextractor(image=output_name+'_lacosmic_imarithsub.fits')
            sextractor(image=output_name+'_lacosmic_sub_inv.fits')
            sextractor(image=output_name+'_lacosmic_imarithsub_inv.fits')
        if image_instrument =='acs' : drzstr = '_drc'
        else: drzstr = '_drz'   
        subtract(image=output_name+drzstr+'_sci.fits', template_image=phot_sub_template, output_image=output_name+'_phot', instrument=image_instrument, template_instrument=templ_instrument,imarith=False)
        subtract(image=output_name+drzstr+'_sci.fits', template_image=phot_sub_template, output_image=output_name+'_phot', instrument=image_instrument, template_instrument=templ_instrument,imarith=True)
        sextractor(image=output_name+'_phot_sub.fits')
        sextractor(image=output_name+'_phot_imarithsub.fits')
        sextractor(image=output_name+'_phot_sub_inv.fits')
        sextractor(image=output_name+'_phot_imarithsub_inv.fits')

def cal_wfc3ir(obsdate,ref='/scratch/snmct/ref/wfc3_ir/',mtab='/scratch/snmct/ref/mtab/',crwfc3comp='/scratch/snmct/ref/comp/wfc3/'):
    #get all of the dark frames from 12349 and 12380 proposals from the archive, and put them in the ref/dark directory
    #this in principle could be down programmatically, but do it by hand for now
    
    #use 0.1 cts/s as a threshold for hot pixels
    dark_thresh=0.1
    #16 in the dq header for hot pixel
    dark_bit = 16
    #persistance threshold
    pers_thresh=30000.0
    #persistance dq bit
    pers_bit=1024
    
    #set iref variable
    os.environ['iref']=ref
    os.environ['mtab']=mtab
    os.environ['crwfc3comp']=crwfc3comp
    #get the current directory to come back to at the end
    this_dir=os.getcwd()+'/'
    #we only want raws, so delete flt, ima, and tra files
    os.system('rm -rf '+ref+'dark/*flt.fits')
    os.system('rm -rf '+ref+'dark/*ima.fits')
    os.system('rm -rf '+ref+'dark/*tra.fits')
    os.system('rm -rf '+ref+'dark/*asn.fits')
    
    #move the reference files
    #move the table reference files to mtab
    os.system('find '+ref+'dark/ -mindepth 1 -not -name "i*.fits" -name "*tmc.fits" -exec mv -f {} '+mtab+' \;')
    os.system('find '+ref+'dark/ -mindepth 1 -not -name "i*.fits" -name "*tmg.fits" -exec mv -f {} '+mtab+' \;')
    os.system('find '+ref+'dark/ -mindepth 1 -not -name "i*.fits" -name "*.fits" -exec mv -f {} '+ref+' \;')
    
    #vet out the darks that don't meet our criteria for "guard dark"
    imgs_raw=tools.parseinput.parseinput(ref+'dark/*raw.fits')[0]
    
    darks=[]
   
    for img in imgs_raw:
        print(img)
        primesi=pyfits.getval(img,'PRIMESI',0)
        subarray=pyfits.getval(img,'SUBARRAY',0)
        subtype=pyfits.getval(img,'SUBTYPE',0)
        samp_seq=pyfits.getval(img,'SAMP_SEQ',0)
        expstart=float(pyfits.getval(img,'EXPSTART',0))
        
        if( primesi=='WFC3' and subarray in [False,'F','False']  and 
        subtype=='FULLIMAG' and samp_seq in ['SPARS50','SPARS100','SPARS200'] ) :
            darks.append((img,expstart))
        else:
            os.system('rm -rf '+img)
            os.system('rm -rf '+img[:-9]+'_spt.fits')
            os.system('rm -rf '+img[:-9]+'_trl.fits')

    # find all of the wfc3 ir exposures taken within the given obsdate (+-1 day), this is a brute force way to do this
    #assume we are in the top directory so path will be ./obj/wfc3_ir/raw_filter_date/
    l_obsdate=int(obsdate)
    date_filt = []
    #get object name
    cmd = "find */wfc3_ir -name 'raw_*"+obsdate+"*' -type d"
    for f in os.popen(cmd).readlines():
        #get just the folder name
        str_arr=f.split('/')
        #get the obsdate and the filter, instrument=wfc3_ir, and the start time for each of the raw files and object name
        obj=str_arr[2].split('_')
        obj.append(str_arr[0])
        date_filt.append(obj)
    
    cmd = "find */wfc3_ir -name 'raw_*"+str(l_obsdate+1)+"*' -type d"
    for f in os.popen(cmd).readlines():
        #get just the folder name
        str_arr=f.split('/')
        #get the obsdate and the filter, instrument=wfc3_ir, and the start time for each of the raw files
        obj=str_arr[2].split('_')
        obj.append(str_arr[0])
        date_filt.append(obj)
        
    cmd = "find */wfc3_ir -name 'raw_*"+str(l_obsdate-1)+"*' -type d"
    for f in os.popen(cmd).readlines():
        #get just the folder name
        str_arr=f.split('/')
        #get the obsdate and the filter, instrument=wfc3_ir, and the start time for each of the raw files
        obj=str_arr[2].split('_')
        obj.append(str_arr[0])
        date_filt.append(obj)
        
    #now we have all the obsdates and filter combinations
    #find all of the raw.fits files in these directorys, save: fits file name, filter, obsdate,and expstart
    imgs=[]
    for date_and_filter in date_filt:
        this_obj=date_and_filter[3]
        this_date=date_and_filter[2]
        #strip off the newline character
        this_date=this_date[:-1]
        this_filter=date_and_filter[1]
        os.chdir(this_obj+'/wfc3_ir/raw_'+this_filter+'_'+this_date)
        
        os.system('find '+ref+'dark/ -mindepth 1 -not -name "i*.fits" -name "*tmc.fits" -exec mv -f {} '+mtab+' \;')
        os.system('find '+ref+'dark/ -mindepth 1 -not -name "i*.fits" -name "*tmg.fits" -exec mv -f {} '+mtab+' \;')
        os.system('find . -not -name "i*.fits" -name "*.fits" -exec mv -f {} '+ref+' \;')
        raw_files=tools.parseinput.parseinput('*raw.fits')[0]
        
        #Clear out the old versions of the flt files and other stuff so calwf3 works
        os.system('rm -rf *_flt.fits')
        os.system('rm -rf *_ima.fits')
        os.system('rm -rf *_tra.fits')
        os.system('rm -rf *_asn.fits')
        
        for f in raw_files:
            print(f)
            #Get the start time for each of the exposures and strip off the _raw.fits from the filename
            expstart=float(pyfits.getval(f,'EXPSTART'))
            imgs.append((f[:-9],this_filter,this_date,expstart,this_obj))
        
            #run calwfc3 on the observation data
            pyfits.setval(f,'PHOTCORR',value='OMIT')
            wfc3.calwf3(f)
            
        os.chdir('../../../')
    
    #sort the exposures by start time
    imgs=sorted(imgs, key=lambda img_tup: img_tup[3])
    
    #find the most dark frame recent dark fram ref/dark
    darks=sorted(darks, key=lambda dark_tup: dark_tup[1])
    first_img=imgs[0]
    #Grab the last dark file as default
    this_dark=darks[len(darks)-1]
    for i,dark in enumerate(darks):
        if dark[1] > first_img[3]:
            this_dark=darks[i-1]
            break
    
    #run calwfc3 to make the dark flt file
    os.chdir(ref+'dark/')
    dark_arr=this_dark[0].split('/')
    dark=dark_arr[-1][:-9]
    pyfits.setval( dark+'_raw.fits', 'UNITCORR', value='PERFORM' )
    pyfits.setval( dark+'_raw.fits', 'DARKCORR', value='PERFORM' )
    pyfits.setval( dark+'_raw.fits', 'PHOTCORR', value='OMIT')
    wfc3.calwf3(dark+'_raw.fits')
    
    #any value in the flt dark image > threshold make DQ mask bit of 16
    hdulist=pyfits.open(dark+'_flt.fits')
    dark_dq=hdulist['DQ'].data.copy()
    dark_dq[:,:]=0
    dark_data=hdulist['SCI'].data.copy()
    dark_dq[dark_data > dark_thresh]=dark_bit
    hdulist.close()
    
    os.chdir(this_dir)
    #add the dark dq mask to the first exposure dq mask
    first_img=imgs[0]
    hdulist=pyfits.open(first_img[4]+'/wfc3_ir/raw_'+first_img[1]+'_'+first_img[2]+'/'+first_img[0]+'_flt.fits')
    dq=hdulist['DQ'].data
    dq=dq | dark_dq
    hdulist.writeto(first_img[4]+'/wfc3_ir/raw_'+first_img[1]+'_'+first_img[2]+'/'+first_img[0]+'_flt.fits',clobber=True)
    hdulist.close()
    
    #make a blank dq array
    pers_dq=dark_dq.copy()

    #for all but last exposure
    for i in range(0,len(imgs)-1):
        this_img=imgs[i]
        #flag all of the pixels over the threshold to make this persistance dq
        hdulist=pyfits.open(this_img[4]+'/wfc3_ir/raw_'+this_img[1]+'_'+this_img[2]+'/'+this_img[0]+'_flt.fits')
        this_data=hdulist['SCI'].data
        h=hdulist[0].header
        exptime=float(h['EXPTIME'])
        
        #do a bitwise or with the old persistance dq and this persistance dq
        pers_dq[this_data > (pers_thresh/exptime)]=pers_dq[this_data > (pers_thresh/exptime)] | pers_bit 
        hdulist.close()
        #add the peristance dq and dark dq to the **next** (not this one) exposure
        next_img=imgs[i+1]
        hdulist=pyfits.open(next_img[4]+'/wfc3_ir/raw_'+next_img[1]+'_'+next_img[2]+'/'+next_img[0]+'_flt.fits')
        next_dq=hdulist['DQ'].data
        next_dq=next_dq | pers_dq
        hdulist.writeto(next_img[4]+'/wfc3_ir/raw_'+next_img[1]+'_'+next_img[2]+'/'+next_img[0]+'_flt.fits',clobber=True)
        hdulist.close()
    #endfor

def make_registration_template(output_name,image_to_register='',this_ref='',input_template_files='',image_instrument='wfc3_ir',template_instrument='wfc3_ir',ra=0.0,dec=0.0):
    #image to register is the flt file that you would like to register: either this or ra and dec must be given
    if image_to_register=='' and ra==0.0 and dec==0.0:
        print('You must give the path to the image to register or the ra and dec to center on.')
    else:
        if ra == 0.0 and dec==0.0:
            #get the ra and dec from the image to register target
            hdulist=pyfits.open(image_to_register)
            hdr=hdulist[0].header
            ra = float(hdr['RA_TARG'])
            dec = float(hdr['DEC_TARG'])
            hdulist.close()
        #set nx,ny,pix_scale, and input files based on instrument
        if image_instrument == 'wfc3_uvis':
            this_nx=6000
            this_ny=6000
            # native_pix_scale=0.04
            this_pix_scale=0.04
        
        if template_instrument=='wfc3_uvis' and input_template_files=='':
            input_template_files='*_flt.fits'
        
        if image_instrument == 'acs':
            this_nx=6000
            this_ny=6000
            # native_pix_scale=0.0495
            this_pix_scale=0.0495
            
        if input_template_files=='' and template_instrument=='acs':
            input_template_files='*_flt.fits'
            
        if image_instrument=='wfc3_ir' :
            this_nx=2000
            this_ny=2000
            #native_pix_scale=0.09
            this_pix_scale=0.09
            
        if input_template_files=='' and template_instrument=='wfc3_ir':
            input_template_files='*_flt.fits'  
              
        if image_instrument=='wfpc2_pc' :    
            this_nx=2400
            this_ny=2400
            #native_pix_scale=0.046
            this_pix_scale=0.046
            
        if input_template_files=='':
            input_template_files='*_c0m.fits'  
           
        if image_instrument=='wfpc2_wf' :
            this_nx=2400
            this_ny=2400
            #native_pix_scale=0.0996
            this_pix_scale=0.0996
            
        if input_template_files=='':
            input_template_files='*_c0m.fits'  

        #run drizzle with the parameters from the image instrument instead of the other way around
        drizzle(output_name,input_files=input_template_files,ref=this_ref,template_image='',
                instrument=template_instrument,drizra=ra,drizdec=dec,
               pix_scale=this_pix_scale,nx=this_nx,ny=this_ny,acs_cte=False)

def run_cte(img,img_base,do_destripe=True):
    #put this into a function to make sure all the memory gets released
    hdulist=pyfits.open(img)
    hdr=hdulist[0].header
    hdr.update('PCTEFILE','jref$pctefile_101109.fits')
    hdulist.writeto(img,clobber=True)
    hdulist.close()
    hdr=None
    hdulist=None
    
    acstools.PixCteCorr.CteCorr(img)
    os.system('rm -f '+img)
    os.system('mv '+img_base+'_cte.fits '+img)
    if do_destripe :
        acstools.acs_destripe.clean(img,'dstrp',clobber=True)
        os.system('rm -f '+img)
        os.system('mv '+img_base+'_flt_dstrp.fits '+img)
                   
import sqlcl
from numpy import float64
astrometry_index_dir = '/usr/local/astro64/astrometry/data/'
def makesdsscat(img, ra,dec,filt):
    #Strip of the drz_sci.fits of the filename
    outbase = img.split('_')[0]
    #Query the SDSS database for all objects detected above 3 sigma within +- 5' of the ra and dec
    #Sort by brightness in the given filter 
    import pdb; pdb.set_trace()
    qry = sqlcl.query("select ra, dec, %s from star where ra < %0.8f and  ra > %0.8f and dec <  %0.8f and dec > %0.8f and err_%s < 0.3 order by %s"%(filt,ra+1.0/12.0, ra-1.0/12.0,dec+1.0/12.0, dec-1.0/12.0,filt,filt ),url='http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp')
    qrytxt = qry.read()
    qry.close()
    qrylist = [i.split(',') for i in qrytxt.split('\n')[2:-1]]
    qryarr = array(qrylist,dtype=float64)
    c1 = pyfits.Column('RA',format='1D',array = qryarr.transpose()[0],unit = 'deg')
    c2 = pyfits.Column('DEC',format='1D',array = qryarr.transpose()[1],unit='deg')
    c3 = pyfits.Column('MAG',format='1D',array = qryarr.transpose()[2],unit='mag')
    coldefs=pyfits.ColDefs([c1,c2,c3])
    tbhdu = pyfits.new_table(coldefs)
    outfile = outbase+'.cat.fits'
    if os.path.exists(outfile): os.remove(outfile)
    tbhdu.writeto(outbase+'.cat.fits')
    
    #Build the index file
    indexfile = outbase+'.index.fits'
    if os.path.exists(indexfile):os.remove(indexfile)
    os.system('build-astrometry-index -i %s -o %s -P -2 -n 1000 -E -I %s'%(outfile,indexfile,outbase))
    
    #Write our own backend.cfg file
    lines = []
    lines.append('cpulimit 300\n')

    # In which directories should we search for indices?
    lines.append('add_path %s\n'%os.getcwd())
    lines.append('index %s.index\n'%outbase)
    f = open('backend.cfg','w')
    f.writelines(lines)
    f.close()
    
def runastrometry(img, ra, dec, pix_scale):
    #Strip of the drz_sci.fits of the filename
    outfile = img.split('_')[0]+'_reg.fits'
    #run solve field using the new index file
    os.system('solve-field --backend-config ./backend.cfg --no-plots --ra %0.8f --dec %0.8f --radius 1.0 --downsample 4 --overwrite --scale-units arcsecperpix --scale-low %0.3f --scale-high %0.3f --crpix-center --no-tweak --no-fits2fits -d 1-3000 --new-fits '
             %(ra,dec,pix_scale*0.8, pix_scale*1.2) +outfile +' ' +img)
    
cat_filter={'F250W':'u', 'F435W':'g', 'F439W':'g', 'F450W':'g', 'F475W':'g', 'F550M':'g', 'F555W':'g',
       'F606W':'r', 'F625W':'r', 'F675W':'r', 'F702W':'i', 'F775W':'i', 'F814W':'i'}
def register2sdss(img):
    #get the RA and Dec of the center of the frame
    imwcs = WCS(img)
    ra,dec = imwcs.getCentreWCSCoords()
    
    #Get the pixel scale in arcseconds/pixel
    pix_scale = imwcs.getPixelSizeDeg()/3600.0
    
    #get the filter
    instrument = pyfits.getval(img,'INSTRUME')
    if instrument =='ACS':
        filt  = pyfits.getval(img,'FILTER1')
        if filt =='CLEAR1L': filt =  pyfits.getval(img,'FILTER2')
    else: filt = pyfits.getval('FILTER')
    
    #Figure out which SDSS filter to use for the catalog 
    sdss_filt = cat_filter[filt]
    
    #make the catalog
    makesdsscat(img,ra,dec,sdss_filt)
    
    #run astrometry.net
    runastrometry(img,ra,dec,pix_scale)
    