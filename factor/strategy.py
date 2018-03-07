"""
Module that holds all strategy-related functions
"""
import os
import logging
import sys
import re
import numpy as np

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
        # Selfcal of entire field:
        #     - calibration on all sources
        #     - imaging of all sources (excluding outliers)
        max_selfcal_loops = field.parset['calibration_specific']['max_selfcal_loops']
        for i in range(max_selfcal_loops):
            strategy_list.append({'do_calibrate': True})
            if i < 3 and not always_do_slowgain:
                strategy_list.append({'do_slowgain': False})
            else:
                strategy_list.append({'do_slowgain': True})
            strategy_list.append({'do_peel': False})
            if field.sectors[-1].name == 'outlier':
                has_outlier = True
                nr_imaging_sectors = len(field.sectors) - 1
            else:
                has_outlier = False
                nr_imaging_sectors = len(field.sectors)
            if nr_imaging_sectors > 1 or has_outlier:
                strategy_list.append({'do_predict': True})
            strategy_list.append({'do_image': True})
            strategy_list.append({'do_update': True})
            if i < 4:
                strategy_list.append({'do_check': False})
            else:
                strategy_list.append({'do_check': True})

    if field.parset['strategy'] == 'sectorselfcal':
        # Selfcal of target(s) only:
        #     - calibration on all sources
        #     - peeling of non-target sources
        #     - calibration on target(s) only
        #     - imaging of target(s) only
        max_selfcal_loops = field.parset['calibration_specific']['max_selfcal_loops']
        for i in range(max_selfcal_loops):
            strategy_list.append({'do_calibrate': True})
            strategy_list.append({'do_slowgain': True})
            strategy_list.append({'do_peel': True})
            strategy_list.append({'do_predict': False})
            strategy_list.append({'do_image': True})
            strategy_list.append({'do_update': True})
            if i < 1:
                strategy_list.append({'do_check': False})
            else:
                strategy_list.append({'do_check': True})

    if field.parset['strategy'] == 'targetexport':
        # Export of target data:
        #     - no calibration
        #     - peeling of non-target sources
        #     - export of target(s) only
        strategy_list.append({'do_calibrate': False})
        strategy_list.append({'do_peel': True})
        strategy_list.append({'do_predict': False})
        strategy_list.append({'do_image': False})
        strategy_list.append({'do_update': False})
        strategy_list.append({'do_check': False})
        strategy_list.append({'do_export': True})

    return strategy_list

