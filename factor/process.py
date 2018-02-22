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
import factor.parset
from factor.operations.calibrate import Calibrate
from factor.operations.image import Image
from factor.lib.scheduler import Scheduler
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

    # Initialize scheduler
    scheduler = Scheduler(parset)

    # Initialize field object
    field = Field(parset)

    # Self calibrate
    max_selfcal_loops = parset['calibration_specific']['max_selfcal_loops']
    do_slowgain_steps = [True] #[False, False, True, True]
    iter = 0
    for do_slowgain in do_slowgain_steps:
        iter += 1

        # Calibrate
        field.do_slowgain_solve = do_slowgain
        scheduler.run(Calibrate(field, iter))

        # Image
        for sector in field.sectors:
            sector.apply_slowgains = do_slowgain
        scheduler.run([Image(field, sector, iter) for sector in field.sectors])

        # Combine new sky models and group
        field.skymodel = lsmtool.load(field.combine_skymodels())
        field.group_skymodel()

        # Check for convergence or iteration limit
        converged = field.check_selfcal_convergence()
        if iter >= max_selfcal_loops or converged:
            break

    # TODO Make final images if needed


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
