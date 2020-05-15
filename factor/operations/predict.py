"""
Module that holds the Predict class
"""
import os
import logging
from factor.lib.operation import Operation

log = logging.getLogger('factor:calibrate')


class Predict(Operation):
    """
    Operation to predict model data
    """
    def __init__(self, field, index):
        super(Predict, self).__init__(field, name='predict', index=index)

    def set_parset_parameters(self):
        """
        Define parameters needed for the pipeline parset template
        """
        self.parset_parms = {'factor_pipeline_dir': self.factor_pipeline_dir,
                             'do_slowgain_solve': self.field.do_slowgain_solve}

    def set_input_parameters(self):
        """
        Define the pipeline inputs
        """
        # If we have a single imagine sector, we don't need to predict
        # its model data, just that of any outlier sector(s)
        if len(self.field.imaging_sectors) == 1:
            start_sector = 1
        else:
            start_sector = 0
        sector_skymodel = []
        sector_sourcedb = []
        sector_obs_sourcedb = []
        sector_filename = []
        sector_starttime = []
        sector_ntimes = []
        sector_model_filename = []
        sector_patches = []
        for sector in self.field.sectors[start_sector:]:
            # Set sector-dependent parameters
            sector.set_prediction_parameters()
            sector_skymodel.append(sector.predict_skymodel_file)
            sdb = os.path.splitext(sector.predict_skymodel_file)[0]+'.sourcedb'
            sector_sourcedb.append(sdb)
            sector_obs_sourcedb.extend([sdb]*len(self.field.observations))
            sector_filename.extend(sector.get_obs_parameters('ms_filename'))
            sector_model_filename.extend(sector.get_obs_parameters('ms_model_filename'))
            sector_patches.extend(sector.get_obs_parameters('patch_names'))
            sector_starttime.extend(sector.get_obs_parameters('predict_starttime'))
            sector_ntimes.extend(sector.get_obs_parameters('predict_ntimes'))

        # Set observation-specific parameters
        obs_filename = []
        obs_starttime = []
        obs_infix = []
        obs_solint_sec = []
        obs_solint_hz = []
        for obs in self.field.observations:
            obs_filename.append(obs.ms_filename)
            obs_starttime.append(obs.convert_mjd(obs.starttime))
            obs_infix.append(obs.infix)
            obs_solint_sec.append(obs.parameters['solint_fast_timestep'][0] * obs.timepersample)
            obs_solint_hz.append(obs.parameters['solint_slow_freqstep'][0] * obs.channelwidth)

        nr_outliers = len(self.field.outlier_sectors)
        peel_outliers = self.field.peel_outliers
        min_uv_lambda = self.field.parset['imaging_specific']['min_uv_lambda']
        max_uv_lambda = self.field.parset['imaging_specific']['max_uv_lambda']

        self.input_parms = {'sector_filename': sector_filename,
                            'sector_starttime': sector_starttime,
                            'sector_ntimes': sector_ntimes,
                            'sector_model_filename': sector_model_filename,
                            'sector_skymodel': sector_skymodel,
                            'sector_sourcedb': sector_sourcedb,
                            'sector_obs_sourcedb': sector_obs_sourcedb,
                            'sector_patches': sector_patches,
                            'h5parm': self.field.h5parm_filename,
                            'obs_solint_sec': obs_solint_sec,
                            'obs_solint_hz': obs_solint_hz,
                            'min_uv_lambda': min_uv_lambda,
                            'max_uv_lambda': max_uv_lambda,
                            'obs_filename': obs_filename,
                            'obs_starttime': obs_starttime,
                            'obs_infix': obs_infix,
                            'nr_outliers': nr_outliers,
                            'peel_outliers': peel_outliers}
