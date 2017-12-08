"""
Module that holds all the Image class
"""
import os
import logging
from factor.lib.operation import Operation

log = logging.getLogger('factor:image')


class Image(Operation):
    """
    Operation to image the field
    """
    def __init__(self, field, index):
        name = 'Calibrate_{0}'.format(index)
        super(Calibrate, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'calibrate_pipeline.parset'

        # Define extra parameters needed for this operation
        ms_filename = field.get_calibration_parameters('ms_filename')
        starttime = field.get_calibration_parameters('starttime')
        ntimes = field.get_calibration_parameters('ntimes')
        solint_fast_timestep = field.get_calibration_parameters('solint_fast_timestep')
        solint_slow_timestep = field.get_calibration_parameters('solint_slow_timestep')
        solint_fast_freqstep = field.get_calibration_parameters('solint_fast_freqstep')
        solint_slow_freqstep = field.get_calibration_parameters('solint_slow_freqstep')
        output_fast_h5parm = [os.path.join(self.pipeline_parset_dir,
            'fast_phase_{}.h5parm'.format(i)) for i in range(field.nchunks)]
        output_slow_h5parm = [os.path.join(self.pipeline_parset_dir,
            'slow_phase_{}.h5parm'.format(i)) for i in range(field.nchunks)]

        self.parms_dict.update({'ms_filename': ms_filename,
                                'starttime': starttime,
                                'ntimes': ntimes,
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
