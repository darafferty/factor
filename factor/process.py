"""
Module that preforms the processing
"""
import logging
from factor import _logging
from factor.parset import parset_read
from factor.strategy import set_strategy
from factor.operations.calibrate import Calibrate
from factor.operations.image import Image
from factor.operations.peel import Peel
from factor.operations.predict import Predict
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
    parset = parset_read(parset_file)

    # Set up logger
    parset['logging_level'] = logging_level
    _logging.set_level(logging_level)

    # Initialize scheduler
    scheduler = Scheduler(parset)

    # Initialize field object
    field = Field(parset)

    # Set up strategy
    strategy_steps = set_strategy(field)

    # Run strategy
    for iter, step in enumerate(strategy_steps):
        # Calibrate
        if step['do_calibrate']:
            field.do_slowgain_solve = step['do_slowgain']
            op = Calibrate(field, iter+1)
            scheduler.run(op)

        # Peel outlier sources
        if step['do_peel']:
            op = Peel(field, iter+1)
            scheduler.run(op)

        # Predict and subtract sector models
        if step['do_predict']:
            op = Predict(field, iter+1)
            scheduler.run(op)

        # Image the sectors
        if step['do_image']:
            for sector in field.sectors:
                sector.apply_slowgains = step['do_slowgain']
            ops = [Image(field, sector, iter+1) for sector in field.sectors
                   if sector.name != 'outlier']
            scheduler.run(ops)
            field.make_mosaic(iter+1)

        # Update the sky models
        if step['do_update']:
            field.update_skymodels(iter+1)

        # Check for convergence
        if step['do_check']:
            has_converged = field.check_selfcal_convergence()
            if has_converged:
                break

    log.info("Factor has finished :)")
