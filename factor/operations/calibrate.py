"""
Module that holds the Calibrate class
"""
import os
import logging
from factor.lib.operation import Operation
from factor.lib import miscellaneous as misc

log = logging.getLogger('factor:calibrate')


class Calibrate(Operation):
    """
    Operation to calibrate
    """
    def __init__(self, field, index):
        super(Calibrate, self).__init__(field, name='calibrate', index=index)

    def set_parset_parameters(self):
        """
        Define parameters needed for the pipeline parset template
        """
        self.parset_parms = {'factor_pipeline_dir': self.factor_pipeline_dir,
                             'mode': self.field.mode,
                             'do_slowgain_solve': self.field.do_slowgain_solve,
                             'use_beam': self.field.use_beam}

    def set_input_parameters(self):
        """
        Define the pipeline inputs
        """
        self.field.set_obs_parameters()
        timechunk_filename = self.field.get_obs_parameters('timechunk_filename')
        starttime = self.field.get_obs_parameters('starttime')
        ntimes = self.field.get_obs_parameters('ntimes')
        slow_starttime = self.field.get_obs_parameters('slow_starttime')
        slow_ntimes = self.field.get_obs_parameters('slow_ntimes')
        freqchunk_filename = self.field.get_obs_parameters('freqchunk_filename')
        startchan = self.field.get_obs_parameters('startchan')
        nchan = self.field.get_obs_parameters('nchan')
        solint_fast_timestep = self.field.get_obs_parameters('solint_fast_timestep')
        solint_fast_timestep_core = self.field.get_obs_parameters('solint_fast_timestep_core')
        solint_fast_timestep_remote = self.field.get_obs_parameters('solint_fast_timestep_remote')
        solint_slow_timestep = self.field.get_obs_parameters('solint_slow_timestep')
        solint_fast_freqstep = self.field.get_obs_parameters('solint_fast_freqstep')
        solint_slow_freqstep = self.field.get_obs_parameters('solint_slow_freqstep')
        output_fast_h5parm = [str(os.path.join(self.pipeline_working_dir,
                              'fast_phase_{}.h5parm'.format(i)))
                              for i in range(self.field.ntimechunks)]
        self.combined_fast_h5parm = os.path.join(self.pipeline_working_dir,
                                                 'fast_phases.h5parm')
        output_fast_core_h5parm = [str(os.path.join(self.pipeline_working_dir,
                                   'fast_phase_core_{}.h5parm'.format(i)))
                                   for i in range(self.field.ntimechunks)]
        output_fast_remote_h5parm = [str(os.path.join(self.pipeline_working_dir,
                                     'fast_phase_remote_{}.h5parm'.format(i)))
                                     for i in range(self.field.ntimechunks)]
        output_slow_h5parm = [str(os.path.join(self.pipeline_working_dir,
                              'slow_gain_{}.h5parm'.format(i)))
                              for i in range(self.field.nfreqchunks)]
        self.combined_slow_h5parm = os.path.join(self.pipeline_working_dir,
                                                 'slow_gains.h5parm')
        baselines_core = self.get_baselines_core()
        antennaconstraint_core = '[[{}]]'.format(','.join(self.get_superterp_stations()))
        antennaconstraint_remote = '[[{}]]'.format(','.join(self.get_core_stations()))
        calibration_skymodel_file = self.field.calibration_skymodel_file
        calibration_sourcedb = str(os.path.join(self.pipeline_working_dir,
                                                'calibration_skymodel.sourcedb'))
        smoothnessconstraint = self.field.smoothnessconstraint
        maxiter = self.field.maxiter
        propagatesolutions = self.field.propagatesolutions
        stepsize = self.field.stepsize
        tolerance = self.field.tolerance
        uvlambdamin = self.field.solve_min_uv_lambda
        sector_bounds_deg = "'{}'".format(self.field.sector_bounds_deg)
        sector_bounds_mid_deg = "'{}'".format(self.field.sector_bounds_mid_deg)
        self.output_aterms_root = str(os.path.join(self.pipeline_working_dir,
                                                   'diagonal_aterms'))
        self.combined_h5parms = str(os.path.join(self.pipeline_working_dir,
                                                 'combined_solutions.h5'))

        self.input_parms = {'timechunk_filename': timechunk_filename,
                            'freqchunk_filename': freqchunk_filename,
                            'starttime': starttime,
                            'ntimes': ntimes,
                            'slow_starttime': slow_starttime,
                            'slow_ntimes': slow_ntimes,
                            'startchan': startchan,
                            'nchan': nchan,
                            'solint_fast_timestep': solint_fast_timestep,
                            'solint_slow_timestep': solint_slow_timestep,
                            'solint_fast_freqstep': solint_fast_freqstep,
                            'solint_slow_freqstep': solint_slow_freqstep,
                            'output_fast_h5parm': output_fast_h5parm,
                            'combined_fast_h5parm': self.combined_fast_h5parm,
                            'output_slow_h5parm': output_slow_h5parm,
                            'combined_slow_h5parm': self.combined_slow_h5parm,
                            'calibration_skymodel_file': calibration_skymodel_file,
                            'calibration_sourcedb': calibration_sourcedb,
                            'smoothnessconstraint': smoothnessconstraint,
                            'maxiter': maxiter,
                            'propagatesolutions': propagatesolutions,
                            'stepsize': stepsize,
                            'tolerance': tolerance,
                            'uvlambdamin': uvlambdamin,
                            'sector_bounds_deg': sector_bounds_deg,
                            'sector_bounds_mid_deg': sector_bounds_mid_deg,
                            'output_aterms_root': self.output_aterms_root,
                            'combined_h5parms': self.combined_h5parms}

    def get_baselines_core(self):
        """
        Returns DPPP string of baseline selection for core calibration

        Returns
        -------
        baselines : str
            Baseline selection string
        """
        cs = self.get_core_stations()
        non_core = [a for a in self.field.stations if a not in cs]

        return '[CR]*&&;!{}'.format(';!'.join(non_core))

    def get_superterp_stations(self):
        """
        Returns list of superterp station names

        Returns
        -------
        stations : list
            Station names
        """
        if self.field.antenna == 'HBA':
            all_st = ['CS002HBA0', 'CS002HBA1', 'CS003HBA0', 'CS003HBA1', 'CS004HBA0', 'CS004HBA1',
                      'CS005HBA0', 'CS005HBA1', 'CS006HBA0', 'CS006HBA1', 'CS007HBA0', 'CS007HBA1']
        elif self.field.antenna == 'LBA':
            all_st = ['CS002LBA', 'CS003LBA', 'CS004LBA', 'CS005LBA', 'CS006LBA', 'CS007LBA']

        return [a for a in all_st if a in self.field.stations]

    def get_core_stations(self):
        """
        Returns list of station names for core calibration

        Returns
        -------
        stations : list
            Station names
        """
        if self.field.antenna == 'HBA':
            all_core = ['CS001HBA0', 'CS002HBA0', 'CS003HBA0', 'CS004HBA0', 'CS005HBA0', 'CS006HBA0',
                        'CS007HBA0', 'CS011HBA0', 'CS013HBA0', 'CS017HBA0', 'CS021HBA0', 'CS024HBA0',
                        'CS026HBA0', 'CS028HBA0', 'CS030HBA0', 'CS031HBA0', 'CS032HBA0', 'CS101HBA0',
                        'CS103HBA0', 'CS201HBA0', 'CS301HBA0', 'CS302HBA0', 'CS401HBA0', 'CS501HBA0',
                        'RS106HBA0', 'RS205HBA0', 'RS305HBA0', 'RS306HBA0', 'RS503HBA0',
                        'CS001HBA1', 'CS002HBA1', 'CS003HBA1', 'CS004HBA1', 'CS005HBA1', 'CS006HBA1',
                        'CS007HBA1', 'CS011HBA1', 'CS013HBA1', 'CS017HBA1', 'CS021HBA1', 'CS024HBA1',
                        'CS026HBA1', 'CS028HBA1', 'CS030HBA1', 'CS031HBA1', 'CS032HBA1', 'CS101HBA1',
                        'CS103HBA1', 'CS201HBA1', 'CS301HBA1', 'CS302HBA1', 'CS401HBA1', 'CS501HBA1',
                        'RS106HBA1', 'RS205HBA1', 'RS305HBA1', 'RS306HBA1', 'RS503HBA1']
        elif self.field.antenna == 'LBA':
            all_core = ['CS001LBA', 'CS002LBA', 'CS003LBA', 'CS004LBA', 'CS005LBA', 'CS006LBA',
                        'CS007LBA', 'CS011LBA', 'CS013LBA', 'CS017LBA', 'CS021LBA', 'CS024LBA',
                        'CS026LBA', 'CS028LBA', 'CS030LBA', 'CS031LBA', 'CS032LBA', 'CS101LBA',
                        'CS103LBA', 'CS201LBA', 'CS301LBA', 'CS302LBA', 'CS401LBA', 'CS501LBA',
                        'RS106LBA', 'RS205LBA', 'RS305LBA', 'RS306LBA', 'RS503LBA']
        return [a for a in all_core if a in self.field.stations]

    def finalize(self):
        """
        Finalize this operation
        """
        # Get the filenames of the aterm images
        with open(self.output_aterms_root+'.txt', 'r') as f:
            self.field.aterm_image_filenames = f.readlines()

        # Save the solutions
        dst_dir = os.path.join(self.parset['dir_working'], 'solutions', 'calibrate_{}'.format(self.index))
        misc.create_directory(dst_dir)
        self.field.h5parm_filename = os.path.join(dst_dir, 'field-solutions.h5')
        if os.path.exists(self.field.h5parm_filename):
            os.remove(self.field.h5parm_filename)
        if self.field.do_slowgain_solve:
            os.system('cp {0} {1}'.format(self.combined_h5parms, self.field.h5parm_filename))
        else:
            os.system('cp {0} {1}'.format(self.combined_fast_h5parm, self.field.h5parm_filename))
