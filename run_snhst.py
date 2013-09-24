#!/usr/bin/env python
if __name__ == '__main__':

    import sys,os,snhst

    if sys.version_info[1] >= 7 or sys.version_info[0] > 2:
        # Parse the input arguments. You have to use a different package based on the python version: What a pain!
        import argparse
        print('Usage: run_snhst.py target instrument filter obsdate wcs_template template_date template_file template_instrument')
        print('Template instrument defaults to image_instrument.')

        parser = argparse.ArgumentParser(description='Run the pipeline on CLASH data.')
        parser.add_argument('target', metavar='target',help='Name of the target. This will be the folder that the data is processed under.')
        parser.add_argument('instrument', metavar='instrument',help='Name of the instrument. Must be lowercase.')
        parser.add_argument('filter', metavar='filter',help='Filter used for the observations. Must be lowercase')
        parser.add_argument('obsdate', metavar='obsdate',help='Observation date: Format is YYYYMMDD')
        parser.add_argument('wcs_templ', metavar='wcs_template',nargs = '?',help='WCS Template File. The current data will be registered to this image. Full path is required.',default='')
        parser.add_argument('templ_date', metavar='template_date',nargs = '?',help='Date of template epoch. Format is YYYYMMDD. This will automatically locate the correct template files to subtract from the current data. No subtraction will be completed if this is omitted.',default='')
        parser.add_argument('templ_file', metavar='template_file',nargs = '?',help='Template image to subtract from the current data. Full path is required. This overrides the template_date.',default='')
        parser.add_argument('templ_instr', metavar='template_instrument',nargs = '?',help='Template instrument. Default is the image instrument.',default='')
        parser.add_argument('-pixscale', metavar='pixscale', type=float, help='Pixel scale to use for multidrizzle.', default=0.0)
        parser.add_argument('-refdir', metavar='refdir', help='Directory for reference files. The final slash is required.', default='/scratch/snmct/ref/')
        args = parser.parse_args()
        
        target = args.target
        instrument = args.instrument
        filt = args.filter
        obsdate = args.obsdate
        wcs_templ = args.wcs_templ
        template_date = args.templ_date
        template_file = args.templ_file
        template_instrument = args.templ_instr
        pixscale = args.pixscale
        refdir = args.refdir
    else:
        import optparse
        parser = optparse.OptionParser()
        parser.add_option("--pixscale", dest="pixscale",
                          help='Pixel scale to use for multidrizzle.', metavar="pixscale", default=0.0)
        parser.add_option('--refdir', metavar='refdir', help='Directory for reference files. The final slash is required.', default='/scratch/snmct/ref/')
        (options, args) = parser.parse_args()
        
        if len(args) < 4: 
            print('Usage: run_snhst.py target instrument filter obsdate wcs_template template_date template_file template_instrument')
            sys.exit(1)
            #initialize some default variables
        template_date=''
        wcs_templ=''
        template_instrument=''
        template_file=''
        target=args[0]
        instrument=args[1]
        filt=args[2]
        obsdate=args[3]
        if len(args)>=5:
            wcs_templ=args[4]
        if len(args)>=6:
            template_date=args[6]
        if len(args)>=7:
            template_file=args[7]
        if len(args)==8:
            template_instrument=args[8]
        
        pixscale = options.pixscale
        refdir = options.refdir
    
    #Start the main call
    this_dir = os.getcwd()
    os.chdir('raw_'+filt+'_'+obsdate)
    if instrument=='acs':
        data_start='j'
    elif instrument=='wfc3_ir':
        data_start='i'
    elif instrument=='wfc3_uvis':
        data_start='i'

    os.system('find . -name "*.fits" -not -name "'+data_start+'*.fits" -exec mv -f {} '+refdir+instrument+'/ \;')
    os.chdir('../raw_'+filt+'_'+obsdate)
    os.system('rm -rf *_drz.fits')
    os.chdir('../')
    os.system('cp -r raw_'+filt+'_'+obsdate+' '+filt+'_'+obsdate)
    os.system('tar -czvf raw_'+filt+'_'+obsdate+'.tar.gz raw_'+filt+'_'+obsdate)
    os.system('rm -rf raw_'+filt+'_'+obsdate)
    os.chdir(filt+'_'+obsdate)

    if instrument=='wfc3_ir':
        image_ext='_drz_sci.fits'
    else:
        image_ext='_lacosmic.fits'

    if template_date=='':
        sub_templ_file=''
        phot_sub_templ_file=''
    elif template_file=='':
        sub_templ_file=this_dir+'/'+target+'/'+instrument+'/'+filt+'_'+template_date+'/'+filt+'_'+template_date+image_ext
        phot_sub_templ_file=this_dir +'/'+target+'/'+instrument+'/'+filt+'_'+template_date+'/'+filt+'_'+template_date+'_drz_sci.fits'
    else:
        sub_templ_file=template_file
        phot_sub_templ_file=template_file

    if template_instrument=='':
        template_instrument=instrument

    snhst.run(filt+'_'+obsdate,ref=refdir+instrument+'/',sub_template=sub_templ_file,phot_sub_template=phot_sub_templ_file,image_instrument=instrument,templ_instrument=template_instrument,wcs_template=wcs_templ,do_cte=False,pixel_scale = pixscale)
    
    os.chdir('../../')
    os.system('cp '+instrument+'/'+filt+'_'+obsdate+'/*.reg reg/'+instrument+'/')

    if not os.path.exists('./all'):
        os.mkdir('all')

    os.chdir('./all')
    runstr = instrument+'/'+filt+'_'+obsdate+'/'+filt+'_'+obsdate
    os.system('ln -s ../'+runstr+'_drz_sci.fits '+filt+'_'+obsdate+'_drz_sci.fits')    
    os.system('ln -s ../'+runstr+'_drz_weight.fits '+filt+'_'+obsdate+'_drz_weight.fits')

    if instrument!='wfc3_ir':
        os.system('ln -s ../'+runstr+'_lacosmic.fits '+filt+'_'+obsdate+'_lacosmic.fits')    
        os.system('ln -s ../'+runstr+'_lamask.fits '+filt+'_'+obsdate+'_lamask.fits') 
    
    if template_date!='' or template_file!='':
        os.system('ln -s ../'+runstr+'_sub.fits '+filt+'_'+obsdate+'_sub.fits')
        
    os.chdir('../')
    
