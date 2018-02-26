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
        timechunk_filename = field.get_calibration_parameters('timechunk_filename')
        starttime = field.get_calibration_parameters('starttime')
        ntimes = field.get_calibration_parameters('ntimes')
        freqchunk_filename = field.get_calibration_parameters('freqchunk_filename')
        startchan = field.get_calibration_parameters('startchan')
        nchan = field.get_calibration_parameters('nchan')
        solint_fast_timestep = field.get_calibration_parameters('solint_fast_timestep')
        solint_slow_timestep = field.get_calibration_parameters('solint_slow_timestep')
        solint_fast_freqstep = field.get_calibration_parameters('solint_fast_freqstep')
        solint_slow_freqstep = field.get_calibration_parameters('solint_slow_freqstep')
        output_fast_h5parm = [os.path.join(self.pipeline_parset_dir,
                              'fast_phase_{}.h5parm'.format(i))
                              for i in range(field.ntimechunks)]
        output_slow_h5parm = [os.path.join(self.pipeline_parset_dir,
                              'slow_gain_{}.h5parm'.format(i))
                              for i in range(field.nfreqchunks)]

        # If needed (i.e., nsectors > 1), define predict parameters
        sector_filename = []
        sector_skymodel = []
        sector_patches = []
        for sector in field.sectors:
            for obs in sector.observations:
                sector_filename.append(obs.ms_filename)
                sector_skymodel.append(sector.calibration_skymodel_file)
                sector_patches.append(sector.patches)
        sector_patches = '[{}]'.format(';'.join(sector_patches))
        obs_filename = []
        for obs in field.observations:
            obs_filename.append(obs.ms_filename)

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
                                'output_slow_h5parm': output_slow_h5parm,
                                'nsectors': len(field.sectors),
                                'sector_filename': sector_filename,
                                'sector_skymodel': sector_skymodel,
                                'sector_patches': sector_patches,
                                'obs_filename': obs_filename})

    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to field object for later use
        if self.field.do_slowgain_solve:
            self.field.h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                     'combine_all_h5parms_output.mapfile')
        else:
            self.field.h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                     'combine_fast_h5parms_output.mapfile')
