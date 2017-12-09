"""
Module that holds all the Calibrate class
"""
import os
import logging
from factor.lib.operation import Operation

log = logging.getLogger('factor:calibrate')


class Calibrate(Operation):
    """
    Operation to calibrate the field
    """
    def __init__(self, field, index):
        name = 'Calibrate_{0}'.format(index)
        super(Calibrate, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'calibrate_pipeline.parset'

        # Define extra parameters needed for this operation
        timechunk_filename = field.get_obs_parameters('timechunk_filename')
        starttime = field.get_obs_parameters('starttime')
        ntimes = field.get_obs_parameters('ntimes')
        freqchunk_filename = field.get_obs_parameters('freqchunk_filename')
        startchan = field.get_obs_parameters('startchan')
        nchan = field.get_obs_parameters('nchan')
        solint_fast_timestep = field.get_obs_parameters('solint_fast_timestep')
        solint_slow_timestep = field.get_obs_parameters('solint_slow_timestep')
        solint_fast_freqstep = field.get_obs_parameters('solint_fast_freqstep')
        solint_slow_freqstep = field.get_obs_parameters('solint_slow_freqstep')
        output_fast_h5parm = [os.path.join(self.pipeline_parset_dir,
                              'fast_phase_{}.h5parm'.format(i))
                              for i in range(field.ntimechunks)]
        output_slow_h5parm = [os.path.join(self.pipeline_parset_dir,
                              'slow_phase_{}.h5parm'.format(i))
                              for i in range(field.nfreqchunks)]

        self.parms_dict.update({'timechunk_filename': timechunk_filename,
                                'freqchunk_filename': freqchunk_filename,
                                'starttime': starttime,
                                'ntimes': ntimes,
                                'startchan': startchan,
                                'nchan': nchan,
                                'solint_fast_timestep': solint_fast_timestep,
                                'solint_slow_timestep': solint_slow_timestep,
                                'solint_fast_freqstep': solint_fast_freqstep,
                                'solint_slow_freqstep': solint_slow_freqstep,
                                'output_fast_h5parm': output_fast_h5parm,
                                'output_slow_h5parm': output_slow_h5parm})

    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to field object for later use
        self.field.fast_h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                      'fast_h5parm.mapfile')
        self.field.slow_h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                      'slow_h5parm.mapfile')
