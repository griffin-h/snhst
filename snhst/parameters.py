from snhst import fits_utils

global_defaults = {'template_image': None,
                   'ra': None,
                   'dec': None,
                   'rotation': 0.0,
                   'pixel_fraction': 1.0,
                   'clean': True,
                   'num_cores': 8,
                   'sky_subtract': True,
                   'use_tweakshifts': True,
                   'tweakshifts_threshold': 10,
                   'dolphot': {'FitSky': 2,
                               'SkipSky': 2,
                               'RCombine': 1.5,
                               'SkySig': 2.25,
                               'SecondPass': 5,
                               'SigFindMult': 0.85,
                               'MaxIT': 25,
                               'NoiseMult': 0.10,
                               'FSat': 0.999,
                               'ApCor': 1,
                               'RCentroid': 2,
                               'PosStep': 0.25,
                               'dPosMax': 2.5,
                               'SigPSF': 10.0,
                               'PSFres': 1,
                               'Align': 4,
                               'Rotate': 1,
                               'ACSuseCTE': 0,
                               'WFC3useCTE': 0,
                               'WFPC2useCTE': 1,
                               'FlagMask': 7,
                               'SigFind': 2.5,
                               'SigFinal': 3.5,
                               'UseWCS': 1,
                               'AlignIter': 3,
                               'AlignTol': 0.5,
                               'AlignStep': 0.2,
                               'VerboseData': 1,
                               'NegSky': 0,
                               'Force1': 0,
                               'DiagPlotType': 'PNG',
                               'InterpPSFlib': 1}}

instrument_defaults = {'wfc3': {'env_ref': 'iref', 'group': '', 'driz_bits': 0, 'pixel_scale': 0.04,
                                'input_files': '*_flt.fits',
                                'crpars': {'rdnoise': 6.5, 'gain': 1.0, 'saturation': 70000.0,
                                           'sig_clip': 4.0, 'sig_frac': 0.2, 'obj_lim': 6.0}},
                       'acs': {'env_ref': 'jref', 'group': '', 'driz_bits': 0, 'pixel_scale': 0.05,
                               'crpars': {'rdnoise': 6.5,
                                          'gain': 1.0,
                                          'saturation': 70000.0,
                                          'sig_clip': 3.0,
                                          'sig_frac': 0.1,
                                          'obj_lim': 5.0}},
                       'wfpc2': {'env_ref': 'uref', 'group':'', 'driz_bits':0,
                                 'input_files': '*_c0m.fits', 'pixel_scale': 0.046,
                                 'crpars': {'rdnoise': 10.0,
                                            'gain': 7.0,
                                            'saturation': 27000.0,
                                            'sig_clip': 4.0,
                                            'sig_frac': 0.3,
                                            'obj_lim': 6.0}
                                 }}


detector_defaults = {'wfc3_uvis': {'nx': 5200, 'ny': 5200, 'input_files': '*_flc.fits',
                                   'dolphot_sky': {'r_in': 15, 'r_out': 35, 'step': 4,
                                                   'sigma_low': 2.25, 'sigma_high': 2.00},
                                   'dolphot_img': {'apsky': '15 25', 'RAper': 2, 'RChi': 1.5,
                                               'RPSF': 10, 'RSky': '15 35',
                                               'RSky2': '3 6'}},
                     'wfc3_ir': {'driz_bits': 512, 'nx': 1700, 'ny': 1700, 'pixel_scale': 0.09,
                                 'dolphot_sky': {'r_in': 10, 'r_out': 25, 'step': 2,
                                                 'sigma_low': 2.25, 'sigma_high': 2.00},
                                 'dolphot_img': {'apsky': '15 25', 'RAper': 2, 'RChi': 1.5,
                                             'RPSF': 10, 'SkipSky': 1, 'RSky': '8 20',
                                             'RSky2': '3 6'}},
                     'acs_wfc': {'nx': 5200, 'ny': 5200, 'input_files': '*_flc.fits',
                                 'dolphot_sky': {'r_in': 15, 'r_out': 35, 'step': 4,
                                                 'sigma_low': 2.25, 'sigma_high': 2.00},
                                 'dolphot_img': {'apsky': '15 25', 'RAper': 2, 'RChi': 1.5,
                                             'RPSF': 10, 'RSky': '15 35',
                                             'RSky2': '3 6'}
                                 },
                     'wfpc2_wfpc2': {'nx': 5200, 'ny': 5200,
                                     'dolphot_sky': {'r_in': 10, 'r_out': 25, 'step': 2,
                                                     'sigma_low': 2.25, 'sigma_high': 2.00},
                                     'dolphot_img': {'apsky': '15 25', 'RAper': 2, 'RChi': 1.5,
                                                 'RPSF': 10, 'RSky': '15 35',
                                                 'RSky2': '3 6'}}}

subarray_defaults = {'wfc3_uvis_full': {},
                     'wfc3_uvis_sub': {'nx': 1400, 'ny': 1400},
                     # Note that wfc3_ir does not ever do cosmic ray removal
                     # because it automatically flags cosmic rays using up the ramp sampling
                     'wfc3_ir_full': {},

                     'acs_wfc_full': {},
                     'acs_wfc_sub': {'nx': 1400, 'ny': 1400, 'input_files': '*_flt.fits'},

                     'wfpc2_wfpc2_full': {}}


def set_default_parameters(options, defaults):
    for key in defaults:
        if key not in options:
            options[key] = defaults[key]


def get_drizzle_parameters(instrument_string, options):
    options['instrument'] = instrument_string

    set_default_parameters(options, subarray_defaults[instrument_string])

    detector_string = "_".join(instrument_string.split('_')[:2])
    set_default_parameters(options, detector_defaults[detector_string])

    instrument = instrument_string.split('_')[0]
    set_default_parameters(options, instrument_defaults[instrument])

    set_default_parameters(options, global_defaults)

    return options


def get_calcsky_parameters(image, options):
    instrument_string = fits_utils.get_instrument(image)
    detector_string = "_".join(instrument_string.split('_')[:2])
    set_default_parameters(options, detector_defaults[detector_string]['dolphot_sky'])
    return options


def get_dolphot_instrument_parameters(image, options):
    instrument_string = fits_utils.get_instrument(image)
    detector_string = "_".join(instrument_string.split('_')[:2])
    set_default_parameters(options, detector_defaults[detector_string]['dolphot_img'])
    return options
