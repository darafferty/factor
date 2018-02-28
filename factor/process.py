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

        # Mosaic the sector images together if needed
        field.make_mosaic(iter)

        # Update the sky models
        field.update_skymodels(iter)

        # Check for convergence or iteration limit
        has_converged = field.check_selfcal_convergence()
        if iter >= max_selfcal_loops or has_converged:
            break

    log.info("Factor has finished :)")
