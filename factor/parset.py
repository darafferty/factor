"""
Module that holds all parset-related functions
"""
import sys
import os
import glob
import logging
import ConfigParser
from factor._logging import set_log_file
from factor.cluster import find_executables, get_compute_nodes
from astropy.coordinates import Angle

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
        for subdir in ['logs', 'results', 'pipelines', 'regions', 'skymodels', 'images', 'solutions']:
            subdir_path = os.path.join(parset_dict['dir_working'], subdir)
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
    except Exception as e:
        log.critical("Cannot use the working dir {0}: {1}".format(parset_dict['dir_working'], e))
        sys.exit(1)
    if use_log_file:
        set_log_file(os.path.join(parset_dict['dir_working'], 'logs', 'factor.log'))
    log.info("=========================================================\n")
    log.info("Working directory is {}".format(parset_dict['dir_working']))

    # Get all the MS files in the input directory. These are identified by the
    # extensions 'ms' or 'MS'
    ms_files = []
    for exten in ['MS', 'ms']:
        ms_files += glob.glob(os.path.join(parset_dict['dir_ms'], '*.{}'.format(exten)))
    parset_dict['mss'] = sorted(ms_files)
    if len(parset_dict['mss']) == 0:
        log.error('No MS files found in {}!'.format(parset_dict['dir_ms']))
        sys.exit(1)
    log.info("Input MS directory is {}".format(parset_dict['dir_ms']))
    log.info("Working on {} input MS file(s)".format(len(parset_dict['mss'])))

    # Make sure the initial skymodel is present
    if 'initial_skymodel' not in parset_dict:
        log.error('No initial sky model file given. Exiting...')
        sys.exit(1)
    elif not os.path.exists(parset_dict['initial_skymodel']):
        log.error('Initial sky model file "{}" not found. Exiting...'.format(parset_dict['initial_skymodel']))
        sys.exit(1)

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

    # Regroup initial skymodel (default = True)
    if 'regroup_initial_skymodel' in parset_dict:
        parset_dict['regroup_initial_skymodel'] = parset.getboolean('global', 'regroup_initial_skymodel')
    else:
        parset_dict['regroup_initial_skymodel'] = False

    # Define strategy
    if 'strategy' not in parset_dict:
        parset_dict['strategy'] = 'fieldselfcal'

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
    allowed_options = ['dir_working', 'dir_ms', 'chunk_size_sec', 'strategy',
                       'use_compression', 'flag_abstime', 'flag_baseline', 'flag_freqrange',
                       'flag_expr', 'chunk_size_hz', 'initial_skymodel', 'regroup_initial_skymodel']
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
    allowed_options = ['max_selfcal_loops', 'multires_selfcal', 'solve_min_uv_lambda',
                       'solve_tecandphase', 'fast_timestep_sec', 'fast_freqstep_hz',
                       'slow_timestep_sec', 'slow_freqstep_hz']
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

    # Size of area to image when using a grid (default = mean FWHM of the primary beam)
    if 'grid_width_ra_deg' in parset_dict:
        parset_dict['grid_width_ra_deg'] = parset.getfloat('imaging', 'grid_width_ra_deg')
    else:
        parset_dict['grid_width_ra_deg'] = None
    if 'grid_width_dec_deg' in parset_dict:
        parset_dict['grid_width_dec_deg'] = parset.getfloat('imaging', 'grid_width_dec_deg')
    else:
        parset_dict['grid_width_dec_deg'] = None

    # Number of sectors along RA to use in imaging grid (default = 0). The number of sectors in
    # Dec will be determined automatically to ensure the whole area specified with grid_center_ra,
    # grid_center_dec, grid_width_ra_deg, and grid_width_dec_deg is imaged. Set grid_nsectors_ra = 0 to force a
    # single sector for the full area. Multiple sectors are useful for parallelizing the imaging
    # over multiple nodes of a cluster or for computers with limited memory
    if 'grid_nsectors_ra' in parset_dict:
        parset_dict['grid_nsectors_ra'] = parset.getint('imaging', 'grid_nsectors_ra')
    else:
        parset_dict['grid_nsectors_ra'] = 1

    # Center of grid to image (default = phase center of data)
    # grid_center_ra = 14h41m01.884
    # grid_center_dec = +35d30m31.52
    if 'grid_center_ra' in parset_dict:
        parset_dict['grid_center_ra'] = Angle(parset_dict['grid_center_ra']).to('deg').value
    else:
        parset_dict['grid_center_ra'] = None
    if 'grid_center_dec' in parset_dict:
        parset_dict['grid_center_dec'] = Angle(parset_dict['grid_center_dec']).to('deg').value
    else:
        parset_dict['grid_center_dec'] = None

    # Instead of a grid, imaging sectors can be defined individually by specifying
    # their centers and widths. If sectors are specified in this way, they will be
    # used instead of the sector grid. Note that the sectors should not overlap
    # sector_center_ra_list = [14h41m01.884, 14h13m23.234]
    # sector_center_dec_list = [+35d30m31.52, +37d21m56.86]
    # sector_width_ra_deg_list = [0.532, 0.127]
    # sector_width_dec_deg_list = [0.532, 0.127]
    len_list = []
    if 'sector_center_ra_list' in parset_dict:
        val_list = parset_dict['sector_center_ra_list'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [Angle(v).to('deg').value for v in val_list]
        parset_dict['sector_center_ra_list'] = val_list
        len_list.append(len(val_list))
    if 'sector_center_dec_list' in parset_dict:
        val_list = parset_dict['sector_center_dec_list'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [Angle(v).to('deg').value for v in val_list]
        parset_dict['sector_center_dec_list'] = val_list
        len_list.append(len(val_list))
    if 'sector_width_ra_deg_list' in parset_dict:
        val_list = parset_dict['sector_width_ra_deg_list'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['sector_width_ra_deg_list'] = val_list
        len_list.append(len(val_list))
    if 'sector_width_dec_deg_list' in parset_dict:
        val_list = parset_dict['sector_width_dec_deg_list'].strip('[]').split(',')
        if val_list[0] == '':
            val_list = []
        val_list = [float(v) for v in val_list]
        parset_dict['sector_width_dec_deg_list'] = val_list
        len_list.append(len(val_list))

    # Check that all the above options have the same number of entries
    if len(set(len_list)) > 1:
        log.error('The options sector_center_ra_list, sector_center_dec_list, '
            'sector_width_ra_deg_list, and sector_width_dec_deg_list must all '
            'have the same number of entires')
        sys.exit(1)

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
                       'selfcal_multiscale_scales_pixel', 'grid_center_ra', 'grid_center_dec',
                       'grid_width_ra_deg', 'grid_width_dec_deg', 'grid_nsectors_ra',
                       'facet_cellsize_arcsec', 'facet_taper_arcsec', 'facet_robust',
                       'wsclean_image_padding', 'selfcal_min_uv_lambda', 'facet_min_uv_lambda',
                       'selfcal_robust_wsclean', 'wsclean_bl_averaging',
                       'sector_center_ra_list', 'sector_center_dec_list',
                       'sector_width_ra_deg_list', 'sector_width_dec_deg_list',
                       'use_idg', 'idg_mode']
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

    # Radius from phase center within which to consider sources during calibration
    # (default = 2 * FWHM of primary beam of highest-frequency band)
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

    # Check for unused options
    allowed_options = ['max_radius_deg', 'target_ra', 'target_dec',
                       'target_radius_arcmin', 'patch_target_flux_jy']
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

    # Maximum fraction of the total memory per node to use (default = 0.9)
    if 'fmem' in parset_dict:
        parset_dict['fmem'] = parset.getfloat('cluster', 'fmem')
        if parset_dict['fmem'] > 1.0:
            parset_dict['fmem'] = 1.0
    else:
        parset_dict['fmem'] = 0.9

    # Cluster type (default = localhost). Use cluster_type = pbs to use PBS / torque
    # reserved nodes and cluster_type = slurm to use SLURM reserved ones
    if 'cluster_type' not in parset_dict:
        parset_dict['cluster_type'] = 'localhost'
    parset_dict['node_list'] = get_compute_nodes(parset_dict['cluster_type'])

    # Full path to a local disk on the nodes for I/O-intensive processing. The path
    # must be the same for all nodes
    if 'dir_local' not in parset_dict:
        parset_dict['dir_local'] = None
    else:
        parset_dict['dir_local'] = parset_dict['dir_local'].rstrip('/')

    # Check for unused options
    allowed_options = ['ncpu', 'fmem', 'cluster_type', 'dir_local', 'lofarroot',
                       'lofarpythonpath']
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
