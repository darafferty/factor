"""
Module that holds all parset-related functions
"""
import sys
import os
import glob
import logging
import ConfigParser
import numpy as np
from factor._logging import set_log_file
from factor.cluster import find_executables, get_compute_nodes

log = logging.getLogger('factor:parset')


def parset_read(parset_file, use_log_file=True):
    """
    Read a Factor-formatted parset file and return dict of parameters

    Parameters
    ----------
    parset_file : str
        Filename of Factor-formated parset file
    use_log_file : bool, optional
        Use a log file as well as outputing to the screen

    Returns
    -------
    parset_dict : dict
        Dict of parset parameters

    """
    if not os.path.isfile(parset_file):
        log.critical("Missing parset file ({}), I don't know what to do :'(".format(parset_file))
        sys.exit(1)

    log.info("Reading parset file: {}".format(parset_file))
    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)

    # Handle global parameters
    parset_dict = get_global_options(parset)

    # Handle calibration parameters
    parset_dict['calibration_specific'].update(get_calibration_options(parset))

    # Handle imaging parameters
    parset_dict['imaging_specific'].update(get_imaging_options(parset))

    # Handle directions-related parameters
    parset_dict['direction_specific'].update(get_directions_options(parset))

    # Handle cluster-specific parameters
    parset_dict['cluster_specific'].update(get_cluster_options(parset))

    # Handle checkfactor parameters
    parset_dict['checkfactor'].update(get_checkfactor_options(parset))

    # Set up working directory. All output will be placed in this directory
    if not os.path.isdir(parset_dict['dir_working']):
        os.mkdir(parset_dict['dir_working'])
    try:
        os.chdir(parset_dict['dir_working'])
        for subdir in ['logs', 'state', 'results', 'regions', 'skymodels']:
            subdir_path = os.path.join(parset_dict['dir_working'], subdir)
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
    except Exception as e:
        log.critical("Cannot use the working dir {0}: {1}".format(parset_dict['dir_working'], e))
        sys.exit(1)
    if use_log_file:
        set_log_file(os.path.join(parset_dict['dir_working'], 'factor.log'))
    log.info("=========================================================\n")
    log.info("Working directory is {0}".format(parset_dict['dir_working']))

    # Get all the MS files in the input directory. These are identified by the
    # extensions 'ms' or 'MS'
    ms_files = []
    for exten in ['MS', 'ms']:
        ms_files += glob.glob(os.path.join(parset_dict['dir_ms'], '*.{}'.format(exten)))
    parset_dict['mss'] = sorted(ms_files)
    if len(parset_dict['mss']) == 0:
        log.error('No MS files found in {}!'.format(parset_dict['dir_ms']))
        sys.exit(1)
    log.info("Input MS directory is {0}".format(parset_dict['dir_ms']))
    log.info("Working on {} input files.".format(len(parset_dict['mss'])))

    # Check for unused sections
    given_sections = parset._sections.keys()
    allowed_sections = ['global', 'calibration', 'imaging', 'directions',
                        'cluster', 'checkfactor']
    for section in given_sections:
        if section not in allowed_sections:
            log.warning('Section "{}" was given in the parset but is not a valid '
                        'section name'.format(section))

    return parset_dict

def get_global_options(parset):
    """
    Handle the global options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all global options

    """
    parset_dict = parset._sections['global'].copy()
    parset_dict.update({'direction_specific': {}, 'calibration_specific': {},
                        'imaging_specific': {}, 'cluster_specific': {},
                        'checkfactor': {}})

    # Size of time chunks in seconds (default = calculate). Generally, the number of
    # chunks should be at least the number of available nodes.
    if 'chunk_size_sec' in parset_dict:
        parset_dict['chunk_size_sec'] = parset.getfloat('global', 'chunk_size_sec')
    else:
        parset_dict['chunk_size_sec'] = None

    # Size of frequency chunks in hz (default = calculate). Generally, the number of
    # chunks should be at least the number of available nodes.
    if 'chunk_size_hz' in parset_dict:
        parset_dict['chunk_size_hz'] = parset.getfloat('global', 'chunk_size_hz')
    else:
        parset_dict['chunk_size_hz'] = None

    # Use Dysco compression (default = False). Enabling this
    # option will result in less storage usage and signifcanctly faster
    # processing. To use this option, you must have the Dysco library in your
    # LD_LIBRARY_PATH. Note: if enabled, Factor will not make symbolic links to the
    # input data, even if they are shorter than chunk_size_sec, but will copy them
    # instead
    if 'use_compression' in parset_dict:
        parset_dict['use_compression'] = parset.getboolean('global', 'use_compression')
    else:
        parset_dict['use_compression'] = False

    # Use interactive mode (default = False). Factor will ask for confirmation of
    # internally derived DDE calibrators and facets
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
    else:
        parset_dict['interactive'] = False

    # Flagging ranges (default = no flagging). A range of times baselines, and
    # frequencies to flag can be specified (see the DPPP documentation for
    # details of syntax) By default, the ranges are AND-ed to produce the final flags,
    # but a set expression can be specified that controls how the selections are
    # combined
    flag_list = []
    if 'flag_abstime' not in parset_dict:
        parset_dict['flag_abstime'] = None
    else:
        flag_list.append('flag_abstime')
    if 'flag_baseline' not in parset_dict:
        parset_dict['flag_baseline'] = None
    else:
        flag_list.append('flag_baseline')
    if 'flag_freqrange' not in parset_dict:
        parset_dict['flag_freqrange'] = None
    else:
        flag_list.append('flag_freqrange')
    if 'flag_expr' not in parset_dict:
        parset_dict['flag_expr'] = ' and '.join(flag_list)
    else:
        for f in flag_list:
            if f not in parset_dict['flag_expr']:
                log.error('Flag selection "{}" was specified but does not '
                    'appear in flag_expr'.format(f))
                sys.exit(1)

    # Check for unused options
    given_options = parset.options('global')
    allowed_options = ['dir_working', 'dir_ms', 'interactive', 'chunk_size_sec',
                       'use_compression', 'flag_abstime', 'flag_baseline', 'flag_freqrange',
                       'flag_expr', 'chunk_size_hz', 'initial_skymodel']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [global] section of the '
                        'parset but is not a valid global option'.format(option))

    return parset_dict

def get_calibration_options(parset):
    """
    Handle the calibration options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all calibration options

    """
    if 'calibration' in parset._sections.keys():
        parset_dict = parset._sections['calibration']
        given_options = parset.options('calibration')
    else:
        parset_dict = {}
        given_options = []

    # Solve for TEC + scalarphase instead of TEC only
    if 'solve_tecandphase' in parset_dict:
        parset_dict['solve_tecandphase'] = parset.getboolean('calibration', 'solve_tecandphase')
    else:
        parset_dict['solve_tecandphase'] = True

    # Maximum number of cycles of the last step of selfcal to perform (default =
    # 10). The last step is looped until the number of cycles reaches this value or
    # until the improvement in dynamic range over the previous image is less than
    # 1.25%. A separate setting can also be used for the target facet only (allowing
    # one to reduce the number for non-target facets)
    if 'max_selfcal_loops' in parset_dict:
        parset_dict['max_selfcal_loops'] = parset.getint('calibration', 'max_selfcal_loops')
    else:
        parset_dict['max_selfcal_loops'] = 10

    # Use baseline-dependent preaveraging to increase the signal-to-noise of the
    # phase-only solve for sources below this flux density (default = 1.0 Jy). When
    # activated, averaging in time is done to exploit the time coherence in the TEC
    # solutions, and averaging in frequency to exploit the frequency coherence of
    # the beam effects
    if 'preaverage_flux_jy' in parset_dict:
        parset_dict['preaverage_flux_jy'] = parset.getfloat('calibration', 'preaverage_flux_jy')
    else:
        parset_dict['preaverage_flux_jy'] = 1.0

    # Use multi-resolution selfcal that starts at 20 arcsec resolution and increases the
    # resolution in stages to the full resolution (default = False). This method may
    # improve convergence, especially when the starting model is poor
    if 'multires_selfcal' in parset_dict:
        parset_dict['multires_selfcal'] = parset.getboolean('calibration', 'multires_selfcal')
    elif 'multiscale_selfcal' in parset_dict:
        log.warning('Option "multiscale_selfcal" is deprecated and should be changed to "multires_selfcal"')
        parset_dict['multires_selfcal'] = parset.getboolean('calibration', 'multiscale_selfcal')
    else:
        parset_dict['multires_selfcal'] = False

    # Minimum uv distance in lambda for calibration (default = 80)
    if 'solve_min_uv_lambda' in parset_dict:
        parset_dict['solve_min_uv_lambda'] = parset.getfloat('calibration', 'solve_min_uv_lambda')
    else:
        parset_dict['solve_min_uv_lambda'] = 80.0

    # Smooth amplitudes with spline fit + 2-D median (default = True, i.e., smooth
    # with a 1-D median only)
    if 'spline_smooth2d' in parset_dict:
        parset_dict['spline_smooth2d'] = parset.getboolean('calibration', 'spline_smooth2d')
    else:
        parset_dict['spline_smooth2d'] = True

    # Include XY and YX correlations during the slow gain solve for sources above
    # this flux density (default = 1000.0). Below this value, only the XX and YY
    # correlations are included. Note that spline_smooth2D must be True to use this
    # option
    if 'solve_all_correlations_flux_jy' in parset_dict:
        parset_dict['solve_all_correlations_flux_jy'] = parset.getfloat('calibration', 'solve_all_correlations_flux_jy')
    else:
        parset_dict['solve_all_correlations_flux_jy'] = 1000.0

    # Solution intervals
    if 'fast_timestep_sec' in parset_dict:
        parset_dict['fast_timestep_sec'] = parset.getfloat('calibration', 'fast_timestep_sec')
    else:
        parset_dict['fast_timestep_sec'] = 8.0
    if 'fast_freqstep_hz' in parset_dict:
        parset_dict['fast_freqstep_hz'] = parset.getfloat('calibration', 'fast_freqstep_hz')
    else:
        parset_dict['fast_freqstep_hz'] = 1e6
    if 'slow_timestep_sec' in parset_dict:
        parset_dict['slow_timestep_sec'] = parset.getfloat('calibration', 'slow_timestep_sec')
    else:
        parset_dict['slow_timestep_sec'] = 600.0
    if 'slow_freqstep_hz' in parset_dict:
        parset_dict['slow_freqstep_hz'] = parset.getfloat('calibration', 'slow_freqstep_hz')
    else:
        parset_dict['slow_freqstep_hz'] = 1e6

    # Check for unused options
    allowed_options = ['max_selfcal_loops', 'preaverage_flux_jy', 'multiscale_selfcal',
                       'multires_selfcal', 'solve_min_uv_lambda', 'spline_smooth2d',
                       'solve_all_correlations_flux_jy', 'solve_tecandphase',
                       'fast_timestep_sec', 'fast_freqstep_hz', 'slow_timestep_sec',
                       'slow_freqstep_hz']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [calibration] section of the '
                        'parset but is not a valid calibration option'.format(option))

    return parset_dict

def get_imaging_options(parset):
    """
    Handle the imaging options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all imaging options

    """
    if 'imaging' in parset._sections.keys():
        parset_dict = parset._sections['imaging']
        given_options = parset.options('imaging')
    else:
        parset_dict = {}
        given_options = []

    # Number of sectors to use in imaging
    if 'nsectors_per_side' in parset_dict:
        parset_dict['nsectors_per_side'] = parset.getint('imaging', 'nsectors_per_side')
    else:
        parset_dict['nsectors_per_side'] = 1

    # Use IDG (image domain gridder) in WSClean. The mode can be cpu, gpu, or hybrid.
    if 'use_idg' in parset_dict:
        parset_dict['use_idg'] = parset.getboolean('imaging', 'use_idg')
    else:
        parset_dict['use_idg'] = True
    if 'idg_mode' not in parset_dict:
        parset_dict['idg_mode'] = 'hybrid'

    # Use baseline-dependent averaging in WSClean (default = True). If enabled,
    # this option can dramatically speed up imaging with WSClean.
    if 'wsclean_bl_averaging' in parset_dict:
        parset_dict['wsclean_bl_averaging'] = parset.getboolean('imaging', 'wsclean_bl_averaging')
    else:
        parset_dict['wsclean_bl_averaging'] = True

    # Max desired peak flux density reduction at center of the facet edges due to
    # bandwidth smearing (at the mean frequency) and time smearing (default = 0.15 =
    # 15% reduction in peak flux). Higher values result in shorter run times but
    # more smearing away from the facet centers. This value only applies to the
    # facet imaging (selfcal always uses a value of 0.15)
    if 'max_peak_smearing' in parset_dict:
        parset_dict['max_peak_smearing'] = parset.getfloat('imaging', 'max_peak_smearing')
    else:
        parset_dict['max_peak_smearing'] = 0.15

    # List of scales in pixels to use when multiscale clean is activated (default =
    # auto). Note that multiscale clean is activated for a direction only when the
    # calibrator or a source in the facet is determined to be larger than 4 arcmin,
    # the facet contains the target (specified below with target_ra and target_dec),
    # or mscale_selfcal_do / mscale_facet_do is set for the direction in the
    # directions file
    if 'selfcal_multiscale_scales_pixel' in parset_dict:
        val_list = parset_dict['selfcal_multiscale_scales_pixel'].strip('[]').split(',')
        str_list = ','.join([v.strip() for v in val_list])
        parset_dict['selfcal_multiscale_scales_pixel'] = str_list
    else:
        parset_dict['selfcal_multiscale_scales_pixel'] = None
    if 'facet_multiscale_scales_pixel' in parset_dict:
        val_list = parset_dict['facet_multiscale_scales_pixel'].strip('[]').split(',')
        str_list = ','.join([v.strip() for v in val_list])
        parset_dict['facet_multiscale_scales_pixel'] = str_list
    else:
        parset_dict['facet_multiscale_scales_pixel'] = None

    # Selfcal imaging parameters: pixel size in arcsec (default = 1.5), Briggs
    # robust parameter (default = -0.5) and minimum uv distance in lambda
    # (default = 80). These settings apply both to selfcal images and to the
    # full facet image used to make the improved facet model that is subtracted
    # from the data
    if 'selfcal_cellsize_arcsec' in parset_dict:
        parset_dict['selfcal_cellsize_arcsec'] = parset.getfloat('imaging', 'selfcal_cellsize_arcsec')
    else:
        parset_dict['selfcal_cellsize_arcsec'] = 1.5
    if 'selfcal_robust' in parset_dict:
        parset_dict['selfcal_robust'] = parset.getfloat('imaging', 'selfcal_robust')
    else:
        parset_dict['selfcal_robust'] = -0.5
    if 'selfcal_min_uv_lambda' in parset_dict:
        parset_dict['selfcal_min_uv_lambda'] = parset.getfloat('imaging', 'selfcal_min_uv_lambda')
    else:
        parset_dict['selfcal_min_uv_lambda'] = 80.0

    # Facet imaging parameters: pixel size in arcsec, Briggs robust parameter, uv
    # taper in arcsec, and minimum uv distance in lambda. These parameters are used
    # only for making full facet images (and not for making improved models). One
    # set of images and one mosaic image will be made for each set of parameters. By
    # default, facets will be imaged using the selfcal imaging parameters above
    len_list = []
    if 'facet_cellsize_arcsec' in parset_dict:
        val_list = parset_dict['facet_cellsize_arcsec'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['facet_cellsize_arcsec'] = val_list
        len_list.append(len(val_list))
    if 'facet_taper_arcsec' in parset_dict:
        val_list = parset_dict['facet_taper_arcsec'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['facet_taper_arcsec'] = val_list
        len_list.append(len(val_list))
    if 'facet_robust' in parset_dict:
        val_list = parset_dict['facet_robust'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['facet_robust'] = val_list
        len_list.append(len(val_list))
    if 'facet_min_uv_lambda' in parset_dict:
        val_list = parset_dict['facet_min_uv_lambda'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['facet_min_uv_lambda'] = val_list
        len_list.append(len(val_list))

    # Check that all the above options have the same number of entries
    if len(set(len_list)) == 0:
        nvals = 1
    elif len(set(len_list)) == 1:
        nvals = len_list[0]
    else:
        log.error('The options facet_cellsize_arcsec, facet_taper_arcsec, facet_robust, and '
            'facet_min_uv_lambda must all have the same number of entires')
        sys.exit(1)

    # Set defaults for any that did not have entries
    if 'facet_cellsize_arcsec' not in parset_dict:
        parset_dict['facet_cellsize_arcsec'] = [parset_dict['selfcal_cellsize_arcsec']] * nvals
    if 'facet_taper_arcsec' not in parset_dict:
        parset_dict['facet_taper_arcsec'] = [0.0] * nvals
    if 'facet_robust' not in parset_dict:
        parset_dict['facet_robust'] = [parset_dict['selfcal_robust']] * nvals
    if 'facet_min_uv_lambda' not in parset_dict:
        parset_dict['facet_min_uv_lambda'] = [80.0] * nvals

    # Padding factor for WSClean images (default = 1.4)
    if 'wsclean_image_padding' in parset_dict:
        parset_dict['wsclean_image_padding'] = parset.getfloat('imaging', 'wsclean_image_padding')
    else:
        parset_dict['wsclean_image_padding'] = 1.4

    # Check for unused options
    allowed_options = ['max_peak_smearing', 'selfcal_cellsize_arcsec', 'selfcal_robust',
                       'selfcal_multiscale_scales_pixel', 'facet_multiscale_scales_pixel',
                       'facet_cellsize_arcsec', 'facet_taper_arcsec', 'facet_robust',
                       'wsclean_image_padding', 'selfcal_min_uv_lambda', 'facet_min_uv_lambda',
                       'selfcal_robust_wsclean', 'wsclean_bl_averaging', 'nsectors_per_side',
                       'fractional_bandwidth_selfcal_facet_image', 'use_idg', 'idg_mode']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [imaging] section of the '
                        'parset but is not a valid imaging option'.format(option))

    return parset_dict

def get_directions_options(parset):
    """
    Handle the directions options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all directions options

    """
    if 'directions' in parset._sections.keys():
        parset_dict = parset._sections['directions']
        given_options = parset.options('directions')
    else:
        parset_dict = {}
        given_options = []

    # Full path to sky model (in makesourcedb format) to be used for calibrator
    # selection and facet-boundary source avoidance (default is to use
    # direction-independent sky model of the highest-frequency band). The sky
    # model must be grouped into patches by source (in PyBDSF, this grouping can be
    # done by setting bbs_patches = 'source' in the write_catalog task)
    if 'faceting_skymodel' not in parset_dict:
        parset_dict['faceting_skymodel'] = None

    # Check whether any sources from the initial subtract sky model fall on facet
    # edges. If any are found, the facet regions are adjusted to avoid them (default
    # is True)
    if 'check_edges' in parset_dict:
        parset_dict['check_edges'] = parset.getboolean('directions', 'check_edges')
    else:
        parset_dict['check_edges'] = True

    # Parameters for selecting directions internally (radius from phase center
    # within which to consider sources as potential calibrators, min flux, max size
    # of a source, and max separation between sources below which they are grouped
    # into one direction; required if no directions_file is given). The number of
    # internally derived directions can be limited to a maximum number of directions
    # if desired with max_num (default = all). Lastly, faceting_radius_deg sets the
    # radius within which facets will be used; outside of this radius, small patches
    # are used that do not appear in the final mosaic.  These parameters will
    # determine the faceting of the field

    # Radius from phase center within which to consider sources as potential
    # calibrators (default = 2 * FWHM of primary beam of highest-frequency band)
    if 'max_radius_deg' in parset_dict:
        parset_dict['max_radius_deg'] = parset.getfloat('directions',
            'max_radius_deg')
    else:
        parset_dict['max_radius_deg'] = None

    # Target flux density in Jy for grouping
    if 'patch_target_flux_jy' in parset_dict:
        parset_dict['patch_target_flux_jy'] = parset.getfloat('directions',
            'patch_target_flux_jy')
    else:
        parset_dict['patch_target_flux_jy'] = 2.5

    # If no directions_file is given, the selection criteria for calibrator sources
    # that follow must be given. For merging of multiple sources into one calibrator
    # group, flux_min_for_merging_Jy (default = 0.1 Jy) and size_max_arcmin set the min
    # flux density and max size of individual sources to be considered for grouping,
    # and separation_max_arcmin sets the max separation between sources below which
    # they are grouped into one calibrator. After grouping, flux_min_Jy sets the
    # min total flux density of a source (or group) to be considered as a DDE
    # calibrator
    if 'flux_min_for_merging_jy' in parset_dict:
        parset_dict['flux_min_for_merging_jy'] = parset.getfloat('directions',
            'flux_min_for_merging_jy')
    else:
        parset_dict['flux_min_for_merging_jy'] = 0.1
    if 'size_max_arcmin' in parset_dict:
        parset_dict['size_max_arcmin'] = parset.getfloat('directions',
            'size_max_arcmin')
    else:
        parset_dict['size_max_arcmin'] = None
    if 'separation_max_arcmin' in parset_dict:
        parset_dict['separation_max_arcmin'] = parset.getfloat('directions',
            'separation_max_arcmin')
    else:
        parset_dict['separation_max_arcmin'] = None
    if 'flux_min_jy' in parset_dict:
        parset_dict['flux_min_jy'] = parset.getfloat('directions',
            'flux_min_jy')
    else:
        parset_dict['flux_min_jy'] = None

    # Number of internally derived directions can be limited to a maximum number
    # of directions if desired with max_num (default = all).
    if 'ndir_max' in parset_dict:
        parset_dict['ndir_max'] = parset.getint('directions', 'ndir_max')
    elif 'max_num' in parset_dict:
        log.warning('Option "max_num" is deprecated and should be changed to "ndir_max"')
        parset_dict['ndir_max'] = parset.getint('directions', 'max_num')
    else:
        parset_dict['ndir_max'] = None

    # Radius within which facets will be used (default = 1.25 * FWHM of primary beam
    # of highest-frequency band); outside of this radius, small patches are used
    # that do not appear in the final mosaic.
    if 'faceting_radius_deg' in parset_dict:
        parset_dict['faceting_radius_deg'] = parset.getfloat('directions',
            'faceting_radius_deg')
    else:
        parset_dict['faceting_radius_deg'] = None

    # Grouping of directions into groups that are selfcal-ed in parallel, defined as
    # grouping:n_total_per_grouping. For example, groupings = 1:5, 4:0 means two
    # groupings are used, with the first 5 directions put into groups of one (i.e.,
    # each direction processed in series) and the rest of the directions divided
    # into groups of 4 (i.e., 4 directions processed in parallel). Default is one at
    # a time (i.e., groupings = 1:0)
    if 'groupings' in parset_dict:
        groupings=[]
        keys = []
        vals = []
        kvs = parset_dict['groupings'].split(',')
        for kv in kvs:
            key, val = kv.split(':')
            keys.append(key.strip())
            vals.append(val.strip())
        for key, val in zip(keys, vals):
            groupings.append({key: int(val)})
        parset_dict['groupings'] = groupings
    else:
        parset_dict['groupings'] = [{'1': 0}]
    log.info("Using the following groupings for directions: {}"
        .format(', '.join(['{0}:{1}'.format(n.keys()[0], n.values()[0])
        for n in parset_dict['groupings']])))

    # If groups are used to process more than one direction in parallel, reordering
    # of the directions in the groups can be done to maximize the flux-weighted
    # separation between directions in each group (default = True)
    if 'allow_reordering' in parset_dict:
        parset_dict['allow_reordering'] = parset.getboolean('directions',
            'allow_reordering')
    else:
        parset_dict['allow_reordering'] = True

    # Total number of directions to selfcal (default = all)
    if 'ndir_selfcal' in parset_dict:
        parset_dict['ndir_selfcal'] = parset.getint('directions', 'ndir_selfcal')
        if parset_dict['ndir_selfcal'] < 1:
            log.error('Total number of directions to selfcal must be 1 or more')
            sys.exit(1)
        log.info("Self calibrating up to %s direction(s)" % (parset_dict['ndir_selfcal']))
    else:
        parset_dict['ndir_selfcal'] = None

    # Total number of directions to process (default = all). If this number is
    # greater than ndir_selfcal, then the remaining directions will not be selfcal-
    # ed but will instead be imaged with the selfcal solutions from the nearest
    # direction for which selfcal succeeded (if a target is specified and
    # target_has_own_facet = True, it will be imaged in this way after ndir_total
    # number of directions are processed)
    if 'ndir_process' in parset_dict:
        parset_dict['ndir_process'] = parset.getint('directions', 'ndir_process')
        if parset_dict['ndir_process'] < 1:
            log.error('Total number of directions to process must be 1 or more')
            sys.exit(1)
        log.info("Processing up to %s direction(s) in total" % (parset_dict['ndir_process']))
    elif 'ndir_total' in parset_dict:
        log.warning('Option "ndir_total" is deprecated and should be changed to "ndir_process"')
        parset_dict['ndir_process'] = parset.getint('directions', 'ndir_total')
        if parset_dict['ndir_process'] < 1:
            log.error('Total number of directions to process must be 1 or more')
            sys.exit(1)
        log.info("Processing up to %s direction(s) in total" % (parset_dict['ndir_process']))
    else:
        parset_dict['ndir_process'] = None

    # A target can be specified to ensure that it falls entirely within a single
    # facet. The values should be those of a circular region that encloses the
    # source and not those of the target itself. Lastly, the target can be placed in
    # a facet of its own. In this case, it will not go through selfcal but will
    # instead use the selfcal solutions of the nearest facet for which selfcal was
    # done
    if 'target_ra' not in parset_dict:
        parset_dict['target_ra'] = None
    if 'target_dec' not in parset_dict:
        parset_dict['target_dec'] = None
    if 'target_radius_arcmin' in parset_dict:
        parset_dict['target_radius_arcmin'] = parset.getfloat('directions',
            'target_radius_arcmin')
    else:
        parset_dict['target_radius_arcmin'] = None
    if 'target_has_own_facet' in parset_dict:
        parset_dict['target_has_own_facet'] = parset.getboolean('directions',
            'target_has_own_facet')
    else:
        parset_dict['target_has_own_facet'] = False

    # Check for unused options
    allowed_options = ['faceting_skymodel', 'directions_file', 'max_radius_deg',
                       'flux_min_for_merging_jy', 'flux_min_jy', 'size_max_arcmin',
                       'separation_max_arcmin', 'max_num', 'ndir_max',
                       'faceting_radius_deg', 'check_edges', 'ndir_total', 'ndir_process',
                       'ndir_selfcal', 'groupings', 'allow_reordering', 'target_ra', 'target_dec',
                       'target_radius_arcmin', 'target_has_own_facet', 'patch_target_flux_jy']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [directions] section of the '
                        'parset but is not a valid directions option'.format(option))

    return parset_dict

def get_cluster_options(parset):
    """
    Handle the compute cluster options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all cluster options

    """
    if 'cluster' in parset._sections.keys():
        parset_dict = parset._sections['cluster']
        given_options = parset.options('cluster')
    else:
        parset_dict = {}
        given_options = []

    # Paths to the LOFAR software
    if 'lofarroot' not in parset_dict:
        if 'LOFARROOT' in os.environ:
            parset_dict['lofarroot'] = os.environ['LOFARROOT']
        else:
            log.critical("The LOFAR root directory cannot be determined. Please "
                         "specify it in the [cluster] section of the parset as lofarroot")
            sys.exit(1)
    if 'lofarpythonpath' not in parset_dict:
        if parset_dict['lofarroot'] in os.environ['PYTHONPATH']:
            pypaths = os.environ['PYTHONPATH'].split(':')
            for pypath in pypaths:
                if parset_dict['lofarroot'] in pypath:
                    parset_dict['lofarpythonpath'] = pypath
                    break
        else:
            log.critical("The LOFAR Python root directory cannot be determined. "
                         "Please specify it in the [cluster] section of the parset as "
                         "lofarpythonpath")
            sys.exit(1)

    # Paths to required executables
    parset_dict = find_executables(parset_dict)

    # Number of CPUs per node to be used.
    if 'ncpu' in parset_dict:
        parset_dict['ncpu'] = parset.getint('cluster', 'ncpu')
    else:
        import multiprocessing
        parset_dict['ncpu'] = multiprocessing.cpu_count()
    log.info("Using up to %s CPU(s) per node" % (parset_dict['ncpu']))

    # Maximum fraction of the total memory per node to use (default = 0.9)
    if 'fmem' in parset_dict:
        parset_dict['fmem'] = parset.getfloat('cluster', 'fmem')
        if parset_dict['fmem'] > 1.0:
            parset_dict['fmem'] = 1.0
    else:
        parset_dict['fmem'] = 0.9
    log.info("Using up to {0}% of the memory per node".format(parset_dict['fmem']*100.0))

    # Maximum number of directions to process in parallel on each node (default =
    # 1). Note that the number of CPUs (set with the ncpu parameter) and the amount
    # of memory available to WSClean (set with the wsclean_fmem parameter) will be
    # divided among the directions on each node
    if 'ndir_per_node' in parset_dict:
        parset_dict['ndir_per_node'] = parset.getint('cluster', 'ndir_per_node')
    else:
        parset_dict['ndir_per_node'] = 1
    log.info("Processing up to %i direction(s) in parallel per node" %
             (parset_dict['ndir_per_node']))

    # Maximum number of io-intensive threads to run per node. If unset, defaults
    # to sqrt of the number of CPUs that will be used (set with the ncpu
    # parameter). Note that this number will be divided among the directions on
    # each node
    if 'nthread_io' in parset_dict:
        parset_dict['nthread_io'] = parset.getint('cluster',
            'nthread_io')
    else:
        parset_dict['nthread_io'] = int(np.ceil(np.sqrt(parset_dict['ncpu'])))
    log.info("Running up to %i IO-intensive job(s) in parallel per node" %
             (parset_dict['nthread_io']))

    # Type of compute cluster. Use cluster_type = localhost for a single machine, cluster_type = PBS to use
    # the PBS / torque reserved nodes, and cluster_type = SLURM to use SLURM
    # reserved ones. If not given, the cluster_type for a single (i.e.,
    # local) node is used
    if 'cluster_type' not in parset_dict:
        parset_dict['cluster_type'] = 'localhost'
    parset_dict['node_list'] = get_compute_nodes(parset_dict['cluster_type'])

    # Full path to a local disk on the nodes for I/O-intensive processing. The path
    # must be the same for all nodes. A selfcal-only path can also be specified to
    # allow certain selfcal data to be cached in memory by setting it to a ram
    # drive (e.g., /dev/shm). By default, dir_local_selfcal is set to dir_local
    if 'dir_local' not in parset_dict:
        parset_dict['dir_local'] = None
    else:
        parset_dict['dir_local'] = parset_dict['dir_local'].rstrip('/')
    if 'dir_local_selfcal' not in parset_dict:
        parset_dict['dir_local_selfcal'] = parset_dict['dir_local']

    # Check for unused options
    allowed_options = ['ncpu', 'fmem', 'fmem', 'ndir_per_node',
                       'clusterdesc_file', 'cluster_type', 'dir_local', 'dir_local_selfcal',
                       'node_list', 'lofarroot', 'lofarpythonpath', 'nthread_io']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [cluster] section of the '
                        'parset but is not a valid cluster option'.format(option))

    return parset_dict

def get_checkfactor_options(parset):
    """
    Handle the options for checkfactor

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all cluster options

    """
    if 'checkfactor' in parset._sections.keys():
        parset_dict = parset._sections['checkfactor']
        given_options = parset.options('checkfactor')
    else:
        parset_dict = {}
        given_options = []

    if 'facet_viewer' not in parset_dict:
        parset_dict['facet_viewer'] = 'casa'

    if 'ds9_limits' not in parset_dict:
        parset_dict['ds9_limits'] = None

    if 'ds9_frames' not in parset_dict:
        parset_dict['ds9_frames'] = 'current'

    if 'image_display' not in parset_dict:
        parset_dict['image_display'] = 'display -geometry 800x600'
    elif parset_dict['image_display'] == 'display':
        parset_dict['image_display'] = 'display -geometry 800x600'

    if 'ds9_load_regions' in parset_dict:
        parset_dict['ds9_load_regions'] = parset.getboolean('checkfactor', 'ds9_load_regions')
    else:
        parset_dict['ds9_load_regions'] = False

    # Check for unused options
    allowed_options = ['facet_viewer', 'ds9_limits', 'ds9_frames', 'image_display',
                       'ds9_load_regions']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [checkfactor] section of the '
                        'parset but is not a valid checkfactor option'.format(option))

    return parset_dict
