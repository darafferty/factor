"""
Module that holds the Predict class
"""
import logging
from factor.lib.operation import Operation

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

        # Define extra parameters needed for this operation
        # If we have a single imagine sector, we don't need to predict
        # its model data, just that of any outlier sector(s)
        if len(field.imaging_sectors) == 1:
            start_sector = 1
        else:
            start_sector = 0
        sector_skymodel = []
        sector_filename = []
        sector_model_filename = []
        sector_patches = []
        for sector in field.sectors[start_sector:]:
            sector_skymodel.append(sector.predict_skymodel_file)  # just one per sector
            sector_filename.extend(sector.get_obs_parameters('ms_filename'))
            sector_model_filename.extend(sector.get_obs_parameters('ms_model_filename'))
            sector_patches.extend(sector.get_obs_parameters('patch_names'))
        sector_patches = '[{}]'.format(';'.join(sector_patches))  # convert to ;-separated list
        obs_filename = []
        for obs in field.observations:
            obs_filename.append(obs.ms_filename)
        nr_outliers = len(field.outlier_sectors)
        peel_outliers = field.peel_outliers

        self.parms_dict.update({'sector_filename': sector_filename,
                                'sector_model_filename': sector_model_filename,
                                'sector_skymodel': sector_skymodel,
                                'sector_patches': sector_patches,
                                'obs_filename': obs_filename,
                                'nr_outliers': nr_outliers,
                                'peel_outliers': peel_outliers})
