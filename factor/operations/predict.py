"""
Module that holds the Predict class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.utilities import create_directory

log = logging.getLogger('factor:calibrate')


class Predict(Operation):
    """
    Operation to predict model data
    """
    def __init__(self, field, index):
        name = 'Predict_{0}'.format(index)
        super(Predict, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'predict_pipeline.parset'

        # Determine if we have any outlier sectors
        nr_outliers = sum([1 for sector in field.sectors if sector.is_outlier])
        nr_imaging_sectors = len(field.sectors) - nr_outliers
        peel_outliers = field.peel_outliers

        # If we have a single sector in addition to outliers, we don't need to predict
        # its model data, just that of the outlier sector(s)
        if nr_imaging_sectors == 1:
            start_sector = 1
        else:
            start_sector = 0
        sector_filename = []
        sector_skymodel = []
        sector_patches = []
        for sector in field.sectors[start_sector:]:
            for obs in sector.observations:
                sector_filename.append(obs.ms_filename)
            sector_skymodel.append(sector.predict_skymodel_file)
            sector_patches.append(sector.patches)
        sector_patches = '[{}]'.format(';'.join(sector_patches))
        obs_filename = []
        for obs in field.observations:
            obs_filename.append(obs.ms_filename)

        self.parms_dict.update({'sector_filename': sector_filename,
                                'sector_skymodel': sector_skymodel,
                                'sector_patches': sector_patches,
                                'obs_filename': obs_filename,
                                'nr_outliers': nr_outliers,
                                'peel_outliers': peel_outliers})
