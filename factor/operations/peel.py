"""
Module that holds the Peel class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.utilities import create_directory

log = logging.getLogger('factor:calibrate')


class Peel(Operation):
    """
    Operation to peel sources
    """
    def __init__(self, field, index):
        name = 'Peel_{0}'.format(index)
        super(Peel, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'peel_pipeline.parset'

        # Set parameters
        outlier_sector = field.sectors[-1]
        sector_skymodel = outlier_sector.predict_skymodel_file
        sector_patches = outlier_sector.patches
        sector_filename = []
        for obs in outlier_sector.observations:
            sector_filename.append(obs.ms_filename)
        obs_filename = []
        for obs in field.observations:
            obs_filename.append(obs.ms_filename)

        self.parms_dict.update({'sector_filename': sector_filename,
                                'sector_skymodel': sector_skymodel,
                                'sector_patches': sector_patches,
                                'obs_filename': obs_filename,
                                'has_outlier': has_outlier})
