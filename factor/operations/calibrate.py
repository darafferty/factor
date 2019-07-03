"""
Module that holds the Calibrate class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.utilities import create_directory

log = logging.getLogger('factor:calibrate')


class Calibrate(Operation):
    """
    Operation to calibrate
    """
    def __init__(self, field, index):
        name = 'Calibrate_{0}'.format(index)
        super(Calibrate, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'calibrate_pipeline.parset'

        # Define extra parameters needed for this operation
        field.set_obs_parameters()
        timechunk_filename = field.get_obs_parameters('timechunk_filename')
        starttime = field.get_obs_parameters('starttime')
        ntimes = field.get_obs_parameters('ntimes')
        slow_starttime = field.get_obs_parameters('slow_starttime')
        slow_ntimes = field.get_obs_parameters('slow_ntimes')
        freqchunk_filename = field.get_obs_parameters('freqchunk_filename')
        startchan = field.get_obs_parameters('startchan')
        nchan = field.get_obs_parameters('nchan')
        solint_fast_timestep = field.get_obs_parameters('solint_fast_timestep')
        solint_fast_timestep_core = field.get_obs_parameters('solint_fast_timestep_core')
        solint_fast_timestep_remote = field.get_obs_parameters('solint_fast_timestep_remote')
        solint_slow_timestep = field.get_obs_parameters('solint_slow_timestep')
        solint_fast_freqstep = field.get_obs_parameters('solint_fast_freqstep')
        solint_slow_freqstep = field.get_obs_parameters('solint_slow_freqstep')
        output_fast_h5parm = [str(os.path.join(self.pipeline_parset_dir,
                              'fast_phase_{}.h5parm'.format(i)))
                              for i in range(field.ntimechunks)]
        output_fast_core_h5parm = [str(os.path.join(self.pipeline_parset_dir,
                                   'fast_phase_core_{}.h5parm'.format(i)))
                                   for i in range(field.ntimechunks)]
        output_fast_remote_h5parm = [str(os.path.join(self.pipeline_parset_dir,
                                     'fast_phase_remote_{}.h5parm'.format(i)))
                                     for i in range(field.ntimechunks)]
        output_slow_h5parm = [str(os.path.join(self.pipeline_parset_dir,
                              'slow_gain_{}.h5parm'.format(i)))
                              for i in range(field.nfreqchunks)]
        baselines_core = self.get_baselines_core()
        antennaconstraint_core = '[[{}]]'.format(','.join(self.get_superterp_stations()))
        antennaconstraint_remote = '[[{}]]'.format(','.join(self.get_core_stations()))

        self.parms_dict.update({'timechunk_filename': timechunk_filename,
                                'freqchunk_filename': freqchunk_filename,
                                'starttime': starttime,
                                'ntimes': ntimes,
                                'slow_starttime': slow_starttime,
                                'slow_ntimes': slow_ntimes,
                                'startchan': startchan,
                                'nchan': nchan,
                                'solint_fast_timestep': solint_fast_timestep,
                                'solint_fast_timestep_core': solint_fast_timestep_core,
                                'solint_fast_timestep_remote': solint_fast_timestep_remote,
                                'solint_slow_timestep': solint_slow_timestep,
                                'solint_fast_freqstep': solint_fast_freqstep,
                                'solint_slow_freqstep': solint_slow_freqstep,
                                'output_fast_h5parm': output_fast_h5parm,
                                'output_fast_core_h5parm': output_fast_core_h5parm,
                                'output_fast_remote_h5parm': output_fast_remote_h5parm,
                                'output_slow_h5parm': output_slow_h5parm,
                                'baselines_core': baselines_core,
                                'antennaconstraint_core': antennaconstraint_core,
                                'antennaconstraint_remote': antennaconstraint_remote})

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
                        'RS106HBA0', 'RS205HBA0', 'RS305HBA0', 'RS306HBA0', 'RS406HBA0', 'RS407HBA0',
                        'RS503HBA0',
                        'CS001HBA1', 'CS002HBA1', 'CS003HBA1', 'CS004HBA1', 'CS005HBA1', 'CS006HBA1',
                        'CS007HBA1', 'CS011HBA1', 'CS013HBA1', 'CS017HBA1', 'CS021HBA1', 'CS024HBA1',
                        'CS026HBA1', 'CS028HBA1', 'CS030HBA1', 'CS031HBA1', 'CS032HBA1', 'CS101HBA1',
                        'CS103HBA1', 'CS201HBA1', 'CS301HBA1', 'CS302HBA1', 'CS401HBA1', 'CS501HBA1',
                        'RS106HBA1', 'RS205HBA1', 'RS305HBA1', 'RS306HBA1', 'RS406HBA1', 'RS407HBA1',
                        'RS503HBA1']
        elif self.field.antenna == 'LBA':
            all_core = ['CS001LBA', 'CS002LBA', 'CS003LBA', 'CS004LBA', 'CS005LBA', 'CS006LBA',
                        'CS007LBA', 'CS011LBA', 'CS013LBA', 'CS017LBA', 'CS021LBA', 'CS024LBA',
                        'CS026LBA', 'CS028LBA', 'CS030LBA', 'CS031LBA', 'CS032LBA', 'CS101LBA',
                        'CS103LBA', 'CS201LBA', 'CS301LBA', 'CS302LBA', 'CS401LBA', 'CS501LBA',
                        'RS106LBA', 'RS205LBA', 'RS305LBA', 'RS306LBA', 'RS406LBA', 'RS407LBA',
                        'RS503LBA']
        return [a for a in all_core if a in self.field.stations]

    def finalize(self):
        """
        Finalize this operation
        """
        # Save output mapfiles for later use
        if self.field.do_slowgain_solve:
            self.field.h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                     'combine_all_h5parms_output.mapfile')
        else:
            self.field.h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                     'combine_fast_h5parms_output.mapfile')
        self.field.fast_aterms_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                      'fast_aterm_images.mapfile')
        self.field.slow_aterms_mapfile = os.path.join(self.pipeline_mapfile_dir,
                                                      'slow_aterm_images.mapfile')

        # Save the solutions
        dst_dir = os.path.join(self.parset['dir_working'], 'solutions', 'calibrate_{}'.format(self.index))
        create_directory(dst_dir)
        dst = os.path.join(dst_dir, 'field-solutions.h5')
        if os.path.exists(dst):
            os.remove(dst)
        sol_map = DataMap.load(self.field.h5parm_mapfile)
        os.system('cp {0} {1}'.format(sol_map[0].file, dst))

#         dst = os.path.join(dst_dir, 'fast_aterms.fits')
#         if os.path.exists(dst):
#             os.remove(dst)
#         sol_map = DataMap.load(self.field.fast_aterms_mapfile)
#         os.system('cp {0} {1}'.format(sol_map[0].file, dst))
#         dst = os.path.join(dst_dir, 'slow_aterms.fits')
#         if os.path.exists(dst):
#             os.remove(dst)
#         sol_map = DataMap.load(self.field.slow_aterms_mapfile)
#         os.system('cp {0} {1}'.format(sol_map[0].file, dst))
