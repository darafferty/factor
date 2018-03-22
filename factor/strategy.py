"""
Module that holds all strategy-related functions
"""
import os
import logging
import sys

log = logging.getLogger('factor:strategy')


def set_strategy(field):
    """
    Sets up the processing strategy

    Parameters
    ----------
    field : Field object
        Field object

    Returns
    -------
    strategy_steps : list
        List of strategy parameters
    """
    strategy_list = []
    always_do_slowgain = True

    if field.parset['strategy'] == 'fieldselfcal':
        # Selfcal without peeling of non-imaged sources:
        #     - calibration on all sources
        #     - imaging of sectors
        max_selfcal_loops = field.parset['calibration_specific']['max_selfcal_loops']
        for i in range(max_selfcal_loops):
            strategy_list.append({})
            strategy_list[i]['do_calibrate'] = True
            if field.initial_h5parm is not None and i == 0:
                strategy_list[i]['do_calibrate'] = False
            if i < 3 and not always_do_slowgain:
                strategy_list[i]['do_slowgain'] = False
            else:
                strategy_list[i]['do_slowgain'] = True
            strategy_list[i]['peel_outliers'] = False
            nr_outlier_sectors = len(field.outlier_sectors)
            nr_imaging_sectors = len(field.imaging_sectors)
            if nr_imaging_sectors > 1 or nr_outlier_sectors > 1:
                strategy_list[i]['do_predict'] = True
            else:
                strategy_list[i]['do_predict'] = False
            strategy_list[i]['do_image'] = True
            if i == max_selfcal_loops - 1:
                strategy_list[i]['do_update'] = False
            else:
                strategy_list[i]['do_update'] = True
            strategy_list[i]['do_update'] = True
            if i < 1 or i == max_selfcal_loops - 1:
                strategy_list[i]['do_check'] = False
            else:
                strategy_list[i]['do_check'] = True

    elif field.parset['strategy'] == 'sectorselfcal':
        # Selfcal with peeling of non-imaged sources:
        #     - calibration on all sources
        #     - peeling of non-sector sources
        #     - imaging of sectors
        #     - calibration on sector sources only
        max_selfcal_loops = field.parset['calibration_specific']['max_selfcal_loops']
        for i in range(max_selfcal_loops):
            strategy_list.append({})
            strategy_list[i]['do_calibrate'] = True
            if field.initial_h5parm is not None and i == 0:
                strategy_list[i]['do_calibrate'] = False
            strategy_list[i]['do_slowgain'] = True
            if i < 1:
                strategy_list[i]['peel_outliers'] = True
            else:
                strategy_list[i]['peel_outliers'] = False
            nr_outlier_sectors = len(field.outlier_sectors)
            nr_imaging_sectors = len(field.imaging_sectors)
            if nr_imaging_sectors > 1 or nr_outlier_sectors > 1:
                strategy_list[i]['do_predict'] = True
            else:
                strategy_list[i]['do_predict'] = False
            strategy_list[i]['do_image'] = True
            if i == max_selfcal_loops - 1:
                strategy_list[i]['do_update'] = False
            else:
                strategy_list[i]['do_update'] = True
            if i < 1 or i == max_selfcal_loops - 1:
                strategy_list[i]['do_check'] = False
            else:
                strategy_list[i]['do_check'] = True

    elif field.parset['strategy'] == 'image':
        # Image one or more sectors:
        #     - no calibration
        strategy_list.append({})
        strategy_list[0]['do_calibrate'] = False
        strategy_list[0]['peel_outliers'] = False
        strategy_list[0]['do_predict'] = True
        strategy_list[0]['do_image'] = True
        strategy_list[0]['do_update'] = False
        strategy_list[0]['do_check'] = False

    else:
        log.error('Strategy "{}" not understood. Exiting...'.format(field.parset['strategy']))
        sys.exit(1)

    log.info('Using "{}" processing strategy'.format(field.parset['strategy']))
    return strategy_list

