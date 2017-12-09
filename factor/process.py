"""
Module that preforms the processing
"""
import sys
import os
import shutil
import numpy as np
import logging
import pickle
import collections
import casacore.tables as pt
from lofarpipe.support.data_map import DataMap
import lsmtool
import factor
import factor.directions
import factor.parset
import factor.cluster
from factor.operations.calibrate import Calibrate
#from factor.operations.outlier_ops import *
#from factor.operations.field_ops import *
#from factor.operations.facet_ops import *
from factor.lib.scheduler import Scheduler
from factor.lib.facet import Facet
from factor.lib.field import Field

log = logging.getLogger('factor')


def run(parset_file, logging_level='info'):
    """
    Processes a dataset using facet calibration

    This function runs the operations in the correct order and handles all the
    bookkeeping for the processing

    Parameters
    ----------
    parset_file : str
        Filename of parset containing processing parameters
    logging_level : str, optional
        One of 'degug', 'info', 'warning' in decreasing order of verbosity

    """
    # Read parset
    parset = factor.parset.parset_read(parset_file)

    # Set up logger
    parset['logging_level'] = logging_level
    factor._logging.set_level(logging_level)

    # Set up scheduler
    scheduler = _set_up_scheduler(parset)

    # Initialize field object
    field = Field(parset)

    # Set up facets for imaging
    facets = _set_up_facets(field)

    # Self calibrate
    max_selfcal_loops = parset['calibration_specific']['max_selfcal_loops']
    do_slowgain_steps = [True] #[False, False, True, True]
    iter = 0
    for do_slowgain in do_slowgain_steps:
        iter += 1

        # Calibrate
        field.do_slowgain_solve = do_slowgain
        op = Calibrate(field, iter)
        scheduler.run(op)

        # Image
        ops = [Image(field, sector, iter) for sector in field.sectors]
        scheduler.run(ops)

        # Combine new sky models and group
        field.skymodel = lsmtool.load(field.combine_skymodels())
        field.group_skymodel()

        # Check for convergence or iteration limit
        converged = field.check_selfcal_convergence()
        if iter >= max_selfcal_loops or converged:
            break

    # Mosaic the final sector images together
    if parset['imaging_specific']['make_mosaic']:
        # Set averaging for primary beam generation
        field.avgpb_freqstep = 20
        field.avgpb_timestep = int(120.0 / bands[0].timepersample)

        # Specify appropriate image, mask, and vertices file
        opname ='facetimage_final'
        field.facet_image_filenames = []
        field.facet_vertices_filenames = []
        for d in facets:
            if not d.is_patch:
                facet_image = DataMap.load(d.facet_image_mapfile[opname])[0].file
                field.facet_image_filenames.append(facet_image)
                field.facet_vertices_filenames.append(d.save_file)

        # Do mosaicking
        op = FieldMosaic(parset, bands, field,
                         parset['imaging_specific']['selfcal_cellsize_arcsec'],
                         parset['imaging_specific']['selfcal_robust'], 0.0,
                         parset['imaging_specific']['selfcal_min_uv_lambda'])
        scheduler.run(op)

    log.info("Factor has finished :)")


def _set_up_scheduler(parset):
    """
    Sets up compute parameters and operation scheduler

    Parameters
    ----------
    parset : dict
        Parset containing processing parameters

    Returns
    -------
    scheduler : Scheduler instance
        The operation scheduler used by the run() function

    """
    log.info('Setting up cluster/node parameters...')

    cluster_parset = parset['cluster_specific']
    if cluster_parset['cluster_type'].lower() == 'pbs':
        log.info('Using cluster setting: "PBS".')
        cluster_parset['clusterdesc'] = factor.cluster.make_pbs_clusterdesc()
        cluster_parset['clustertype'] = 'pbs'
    elif cluster_parset['cluster_type'].lower() == 'slurm':
        log.info('Using cluster setting: "SLURM".')
        cluster_parset['clusterdesc'] = factor.cluster.make_slurm_clusterdesc()
        cluster_parset['clustertype'] = 'slurm'
    elif cluster_parset['cluster_type'].lower() == 'juropa_slurm':
        log.info('Using cluster setting: "JUROPA_slurm" (Single '
                 'genericpipeline using multiple nodes).')
        # slurm_srun on JUROPA uses the local.clusterdesc
        cluster_parset['clusterdesc'] = os.path.join(parset['lofarroot'],
                                                     'share', 'local.clusterdesc')
        cluster_parset['clustertype'] = 'juropa_slurm'
        cluster_parset['node_list'] = ['localhost']
    elif cluster_parset['cluster_type'].lower() == 'mpirun':
        log.info('Using cluster setting: "mpirun".')
        # mpirun uses the local.clusterdesc?
        cluster_parset['clusterdesc'] = os.path.join(parset['lofarroot'],
                                                     'share', 'local.clusterdesc')
        cluster_parset['clustertype'] = 'mpirun'
        cluster_parset['node_list'] = ['localhost']
    else:
        log.info('Using cluster setting: "local" (Single node).')
        cluster_parset['clusterdesc'] = cluster_parset['lofarroot'] + '/share/local.clusterdesc'
        cluster_parset['clustertype'] = 'local'
    if not 'node_list' in cluster_parset:
        cluster_parset['node_list'] = factor.cluster.get_compute_nodes(cluster_parset['clusterdesc'])

    # check ulimit(s)
    try:
        import resource
        nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        if cluster_parset['clustertype'] == 'local' and nof_files_limits[0] < nof_files_limits[1]:
            log.debug('Setting limit for number of open files to: {}.'.format(nof_files_limits[1]))
            resource.setrlimit(resource.RLIMIT_NOFILE,(nof_files_limits[1],nof_files_limits[1]))
            nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        log.debug('Active limit for number of open files is {0}, maximum limit '
                  'is {1}.'.format(nof_files_limits[0],nof_files_limits[1]))
        if nof_files_limits[0] < 2048:
            log.warn('The limit for number of open files is small, this could '
                     'result in a "Too many open files" problem when running factor.')
            log.warn('The active limit can be increased to the maximum for the '
                     'user with: "ulimit -Sn <number>" (bash) or "limit descriptors 1024" (csh).')
    except resource.error:
        log.warn('Cannot check limits for number of open files, what kind of system is this?')

    # Get paths to required executables
    factor.cluster.find_executables(parset)

    # Set up scheduler for operations (pipeline runs)
    ndir_simul = len(cluster_parset['node_list']) * cluster_parset['ndir_per_node']
    if parset['direction_specific']['groupings'] is not None:
        ngroup_max = int(max([int(n.items()[0][0]) for n in
                              parset['direction_specific']['groupings']]))
    else:
        ngroup_max = 1
    if ndir_simul < ngroup_max:
        log.warn('The maximum number of facets that can be proccessed '
                 'simultaneously ({0}) is less than the number of facets in the '
                 'largest group ({1}). For best performance, these values should be '
                 'equal'.format(ndir_simul, ngroup_max))
    scheduler = Scheduler(cluster_parset['genericpipeline_executable'], max_procs=ndir_simul)

    return scheduler


def _set_up_facets(field):
    """
    Sets up facets for imaging

    Parameters
    ----------
    field : Field object
        Field for this run

    Returns
    -------
    facets : List of Direction instances
        All facets to be used by the run() function
    direction_groups : List of lists of Direction instances
        Groups of facets to be selfcal-ed

    """
    dir_parset = field.parset['direction_specific']
    max_radius_deg = dir_parset['max_radius_deg']
    if dir_parset['faceting_skymodel'] is not None:
        log.info("Using {} as sky model for source avoidance and DDE calibrator "
            "selection (if desired)".format(dir_parset['faceting_skymodel']))
        initial_skymodel = lsmtool.load(dir_parset['faceting_skymodel'])
    else:
        initial_skymodel = factor.directions.make_initial_skymodel(field)

    log.info('Setting up facets...')
    facets = _initialize_facets(field, initial_skymodel,
        max_radius_deg=max_radius_deg)

    # Check with user
    if field.parset['interactive']:
        print("\nFacet and DDE calibrator regions saved. Please check that they\n"
            "are OK before continuing. You can edit the facets file and\n"
            "continue; FACTOR will pick up any changes to it. Note: if you\n"
            "choose not to continue and you let FACTOR generate the facets\n"
            "internally, you must delete the FACTOR-made facets file\n"
            "(dir_working/factor_directions.txt) before restarting if you want\n"
            "FACTOR to regenerate it\n")
        prompt = "Continue processing (y/n)? "
        answ = raw_input(prompt)
        while answ.lower() not in  ['y', 'n', 'yes', 'no']:
            answ = raw_input(prompt)
        if answ.lower() in ['n', 'no']:
            log.info('Exiting...')
            sys.exit(0)
        else:
            # Continue processing, but first re-initialize the facets to
            # pick up any changes the user made to the facets file
            facets = _initialize_facets(field.parset, initial_skymodel,
                ref_band, max_radius_deg=max_radius_deg, dry_run=dry_run)

    # Load previously completed operations (if any) and facetsubreset-specific
    # attributes and save the state
    for direction in facets:
        direction.load_state()
        direction.save_state()

    # Select subset of facets to process
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_has_own_facet:
        direction_names = [d.name for d in facets]
        target = facets[direction_names.index('target')]
    if dir_parset['ndir_process'] is not None:
        if dir_parset['ndir_process'] < len(facets):
            facets = facets[:dir_parset['ndir_process']]

            # Make sure target is still included
            direction_names = [d.name for d in facets]
            if target_has_own_facet and 'target' not in direction_names:
                facets.append(target)

    # Set various direction attributes
    for i, facet in enumerate(facets):
        # Set direction sky model
        if initial_skymodel is not None:
            facet.set_skymodel(initial_skymodel.copy())

        # Set field center to that of first band (all bands have the same phase
        # center)
        direction.field_ra = field.ra
        direction.field_dec = field.dec

        # Set initial name of column that contains DATA
        direction.subtracted_data_colname = 'DATA'
        direction.use_compression = True

        # Set any flagging parameters
        direction.flag_abstime = field.parset['flag_abstime']
        direction.flag_baseline = field.parset['flag_baseline']
        direction.flag_freqrange = field.parset['flag_freqrange']
        direction.flag_expr = field.parset['flag_expr']

    # Select facets to selfcal, excluding outliers and target
    if target_has_own_facet:
        # Make sure target is not a DDE calibrator and is at end of facets list
        selfcal_facets = [d for d in facets if d.name != target.name and
                              not d.is_outlier and not d.peel_calibrator]
        facets = [d for d in facets if d.name != target.name] + [target]
    else:
        selfcal_facets = [d for d in facets if not d.is_outlier and
            not d.peel_calibrator]

    if dir_parset['ndir_selfcal'] is not None:
        if dir_parset['ndir_selfcal'] <= len(selfcal_facets):
            selfcal_facets = selfcal_facets[:dir_parset['ndir_selfcal']]

    return facets


def _initialize_facets(field, initial_skymodel, max_radius_deg):
    """
    Read in facets file and initialize resulting facets

    Parameters
    ----------
    field : Field object
        Field for this run

    Returns
    -------
    facets : List of Direction instances
        All facets to be used

    """
    dir_parset = field.parset['direction_specific']
    s = None

    # First check for user-supplied facets file, then for Factor-generated
    # file from a previous run, then for parameters needed to generate it internally
    if 'directions_file' in dir_parset:
        facets = factor.directions.directions_read(dir_parset['directions_file'],
            field.parset['dir_working'])
    elif os.path.exists(os.path.join(field.parset['dir_working'], 'factor_directions.txt')):
        facets = factor.directions.directions_read(os.path.join(field.parset['dir_working'],
            'factor_directions.txt'), field.parset['dir_working'])
    else:
        if dir_parset['flux_min_jy'] is None or \
            dir_parset['size_max_arcmin'] is None or \
            dir_parset['separation_max_arcmin'] is None:
                log.critical('If no facets file is specified, you must '
                    'give values for flux_min_Jy, size_max_arcmin, and '
                    'separation_max_arcmin')
                sys.exit(1)
        else:
            # Make facets from dir-indep sky model of highest-frequency
            # band, as it has the smallest field of view
            log.info("No facets file given. Selecting facets internally...")
            s = initial_skymodel.copy()

            # Filter out sources that lie outside of maximum specific radius from phase
            # center
            if max_radius_deg is None:
                max_radius_deg = field.fwhm_deg # means a diameter of 2 * FWHM
            log.info('Removing sources beyond a radius of {0} degrees)...'.format(max_radius_deg))
            dist = s.getDistance(field.ra, field.dec, byPatch=True)
            s.remove(dist > max_radius_deg, aggregate=True)

            # Generate the facets file
            dir_parset['directions_file'] = factor.directions.make_directions_file_from_skymodel(
                s, dir_parset['flux_min_jy'], dir_parset['size_max_arcmin'],
                dir_parset['separation_max_arcmin'],
                directions_max_num=dir_parset['ndir_max'],
                interactive=field.parset['interactive'],
                flux_min_for_merging_Jy=dir_parset['flux_min_for_merging_jy'])
            facets = factor.directions.directions_read(dir_parset['directions_file'],
                field.parset['dir_working'])

    # Add the target to the facets list if desired
    target_ra = dir_parset['target_ra']
    target_dec = dir_parset['target_dec']
    target_radius_arcmin = dir_parset['target_radius_arcmin']
    target_has_own_facet = dir_parset['target_has_own_facet']
    if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
        # Make target object
        target = Facet('target', target_ra, target_dec,
            factor_working_dir = field.parset['dir_working'])
        if target_has_own_facet:
            target.contains_target = True

            # Check if target is already in facets list because it was
            # selected as a DDE calibrator. If so, remove the duplicate
            nearest, dist = factor.directions.find_nearest(target, facets)
            if dist < dir_parset['target_radius_arcmin']/60.0:
                facets.remove(nearest)

            # Add target to facets list
            facets.append(target)
        else:
            # Find direction that contains target
            nearest, dist = factor.directions.find_nearest(target, facets)
            nearest.contains_target = True
    else:
        if target_has_own_facet:
            log.critical('target_has_own_facet = True, but target RA, Dec, or radius not found in parset')
            sys.exit(1)

    # Load previously completed operations (if any) and facetsubreset-specific
    # attributes and save the state
    for direction in facets:
        direction.load_state()
        direction.save_state()

    # Set calibrator size (must be done before faceting below is done)
    for d in facets:
        d.set_cal_size(field.parset['imaging_specific']['selfcal_cellsize_arcsec'])

    # Create facets and patches
    faceting_radius_deg = dir_parset['faceting_radius_deg']
    if faceting_radius_deg is None:
        faceting_radius_deg = 1.25 * field.fwhm_deg / 2.0
    beam_ratio = 1.0 / np.sin(field.mean_el_rad) # ratio of N-S to E-W beam
    if dir_parset['recalculate_facets']:
        factor.directions.thiessen(facets, field.ra, field.dec,
            faceting_radius_deg, s=s, check_edges=dir_parset['check_edges'],
            target_ra=target_ra, target_dec=target_dec,
            target_radius_arcmin=target_radius_arcmin, beam_ratio=beam_ratio)

    # Make DS9 region files so user can check the facets, etc.
    ds9_facet_reg_file = os.path.join(field.parset['dir_working'], 'regions', 'facets_ds9.reg')
    factor.directions.make_ds9_region_file(facets, ds9_facet_reg_file)
    ds9_calimage_reg_file = os.path.join(field.parset['dir_working'], 'regions', 'calimages_ds9.reg')
    factor.directions.make_ds9_calimage_file(facets, ds9_calimage_reg_file)

    return facets


def _get_image_type_and_name(cellsize_arcsec, taper_arcsec, robust, selfcal_robust,
    min_uv_lambda, parset, opbase='facetimage'):
    """
    Checks input parameters and returns type and associated operation

    Parameters
    ----------
    cellsize_arcsec : float
        Cell size
    taper_arcsec : float
        Taper
    robust : float
        Briggs robust
    min_uv_lambda : float
        Min uv cut
    parset : dict
        Factor parset
    opbase : str, optional
        Basename of operation

    Returns
    -------
    full_res_im, opname : bool, str
        Image type (full resolution or not) and FacetImage operation name

    """
    if (cellsize_arcsec != parset['imaging_specific']['selfcal_cellsize_arcsec'] or
        robust != selfcal_robust or taper_arcsec != 0.0 or
        min_uv_lambda != parset['imaging_specific']['selfcal_min_uv_lambda']):
        opname = '{0}_c{1}r{2}t{3}u{4}'.format(opbase, round(cellsize_arcsec, 1),
                round(robust, 2), round(taper_arcsec, 1), round(min_uv_lambda, 1))
        full_res_im = False
    else:
        opname = opbase
        full_res_im = True

    return full_res_im, opname
