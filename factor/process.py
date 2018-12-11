"""
Module that preforms the processing
"""
import logging
from factor import _logging
from factor.parset import parset_read
from factor.strategy import set_strategy
from factor.operations.calibrate import Calibrate
from factor.operations.image import Image
from factor.operations.predict import Predict
from factor.lib.scheduler import Scheduler
from factor.lib.field import Field

log = logging.getLogger('factor')


def run(parset_file, logging_level='info', sectors_to_export=[], export_corrected_data=False):
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

    # Get the processing strategy
    strategy_steps = set_strategy(field)

    # Run the strategy
    for iter, step in enumerate(strategy_steps):
        # Calibrate
        if step['do_calibrate']:
            field.do_slowgain_solve = step['do_slowgain']
            op = Calibrate(field, iter+1)
            scheduler.run(op)

        # Predict and subtract sector models
        if step['do_predict']:
            field.peel_outliers = step['peel_outliers']
            op = Predict(field, iter+1)
            scheduler.run(op)
            if step['peel_outliers']:
                # Update the observations to use the new peeled datasets and remove the
                # outlier sectors (since, once peeled, they are no longer needed)
                for obs in field.observations:
                    obs.ms_filename = obs.ms_field
                    obs.set_calibration_parameters(parset, field.num_patches)
                field.sectors = [sector for sector in field.sectors if not sector.is_outlier]
                field.outlier_sectors = []

        # Image the sectors
        if step['do_image']:
            # Set flag for slow-gain apply
            imaging_sectors = [sector for sector in field.sectors if not sector.is_outlier]
            for sector in imaging_sectors:
                sector.apply_slowgains = step['do_slowgain']

            # Put the sectors using multiscale clean first, as they take the longest to
            # image
            sorted_sectors = [sector for sector in imaging_sectors if sector.multiscale]
            sorted_sectors.extend([sector for sector in imaging_sectors
                                   if not sector.multiscale])
            ops = [Image(field, sector, iter+1) for sector in sorted_sectors]
            scheduler.run(ops)
            field.make_mosaic(iter+1)

        # Check for convergence
        if step['do_check']:
            has_converged = field.check_selfcal_convergence()
            if has_converged:
                break

        # Update the sky models
        if step['do_update']:
            field.update_skymodels(iter+1, step['regroup_model'], step['imaged_sources_only'])

    log.info("Factor has finished :)")
