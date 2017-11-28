"""
Definition of the Observation class
"""
import os
import sys
import shutil
import logging
import casacore.tables as pt
import lofar.parmdb
import numpy as np
import multiprocessing
import itertools


class Observation(object):
    """
    The Observation object contains various MS-related parameters

    Parameters
    ----------
    ms_filename : str
        Filename of the MS file

    """
    def __init__(self, ms_filename, working_dir):

        self.ms_filename = ms_filename
        self.name = os.path.basename(self.ms_filename)
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.working_dir = working_dir
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')

        # Load state (if any)
        if not self.load_state():
            # If state is not current, get observation details
            self.get_info()


    def get_info(self):
        """
        Stores various info from input MS
        """
        # Get time info
        tab = pt.table(self.ms_filename, ack=False)
        self.starttime = np.min(tab.getcol('TIME'))
        self.endtime = np.max(tab.getcol('TIME'))
        self.timepersample = tab.getcell('EXPOSURE',0)
        self.numsamples = len(tab.getcol('TIME'))
        tab.close()

        # Get frequency info
        sw = pt.table(self.ms_filename+'::SPECTRAL_WINDOW', ack=False)
        self.referencefreq = sw.col('REF_FREQUENCY')[0])
        self.startfreq = np.min(sw.col('CHAN_FREQ')[0]))
        self.endfreq = np.max(sw.col('CHAN_FREQ')[0]))
        self.numchannels = sw.col('NUM_CHAN')[0])
        self.channelwidth = sw.col('CHAN_WIDTH')[0][0])
        sw.close()

        # Get pointing info
        obs = pt.table(self.ms_filename+'::FIELD', ack=False)
        self.ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
        if self.ra < 0.:
            self.ra = 360.0 + (self.ra)
        self.dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
        obs.close()

        # Get station diameter
        ant = pt.table(self.ms_filename+'::ANTENNA', ack=False)
        self.diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

        # Find mean elevation and FOV
        tab = pt.table(self.ms_filename, ack=False)
        exiting_colnames = tab.colnames()
        if 'AZEL1' not in exiting_colnames:
            tab.close()
            pt.addDerivedMSCal(self.ms_filename)
            tab = pt.table(self.ms_filename, ack=False)
        el_values = tab.getcol('AZEL1', rowincr=10000)[:, 1]
        tab.close()
        pt.removeDerivedMSCal(self.ms_filename)
        self.mean_el_rad = np.mean(el_values)


    def set_calibration_parameters(self, parset):
        """
        Sets the calibration parameters per calibration time chunk

        Parameters
        ----------
        parset : dict
            Parset with processing parameters

        """
        target_chunksize = parset['chunk_size_sec']
        target_fast_timestep = parset['fast_timestep_sec']
        target_fast_freqstep = parset['fast_freqstep_hz']
        target_slow_timestep = parset['slow_timestep_sec']
        target_slow_freqstep = parset['slow_freqstep_hz']

        self.calibration_parameters = {}

        # Calculate time ranges of calibration chunks
        timepersample = self.timepersample
        numsamples = self.numsamples
        samplesperchunk = int(round(target_chunksize / timepersample))
        chunksize = samplesperchunk * timepersample
        mystarttime = self.starttime
        myendtime = self.endtime
        if (myendtime-mystarttime) > (2.*chunksize):
            nchunks = int((numsamples*timepersample)/chunksize)
        else:
            nchunks = 1
        self.log.debug('Using {1} calibration chunk(s)'.format(nchunks))
        self.calibration_parameters['starttime'] = [mystarttime+(chunksize*i) for i in range(nchunks)]
        self.calibration_parameters['ntimes'] = [samplesperchunk*i for i in range(nchunks)]
        self.calibration_parameters['ntimes'][-1] = 0 # set last entry to extend until end

        # Set solution intervals (same for every calibration chunk)
        freqpersample = self.chan_width_hz
        self.calibration_parameters['solint_fast_timestep'] = [int(round(target_fast_timestep / timepersample))] * nchunks
        self.calibration_parameters['solint_fast_freqstep'] = [self.get_nearest_frequstep(target_fast_freqstep / freqpersample)] * nchunks
        self.calibration_parameters['solint_slow_timestep'] = [int(round(target_slow_timestep / timepersample))] * nchunks
        self.calibration_parameters['solint_slow_freqstep'] = [self.get_nearest_frequstep(target_slow_freqstep / freqpersample)] * nchunks


    def get_nearest_frequstep(self, freqstep):
        """
        Gets the nearest frequstep

        Parameters
        ----------
        freqstep : int
            Target frequency step

        Returns
        -------
        optimum_step : int
            Optimum frequency step nearest to target step

        """
        # first generate a list of possible values for freqstep
        if not hasattr(self, 'freq_divisors'):
            tmp_divisors = []
            for step in range(self.numchannels, 0, -1):
                if (self.numchannels % step) == 0:
                    tmp_divisors.append(step)
            self.freq_divisors = np.array(tmp_divisors)
        idx = np.argmin(np.abs(self.freq_divisors-freqstep))

        return self.freq_divisors[idx]


    def get_md5(self):
        """
        Gets and stores md5 hash for MS file
        """
        self._ms_file_md5 = []

        if not os.path.exists(self.ms_filename):
            self.log.critical('MS file "{0}" not found!'.format(self.ms_filename))
            sys.exit(1)
        filename = '{0}/OBSERVATION/table.f0'.format(self.ms_filename)
        md5hash = hashlib.md5(open(filename, 'rb').read()).hexdigest()
        self._ms_file_md5 = md5hash


    def save_state(self):
        """
        Saves the band state to a file

        """
        import pickle

        with open(self.save_file, 'wb') as f:
            # Remove log object, as it cannot be pickled
            save_dict = self.__dict__.copy()
            save_dict.pop('log')
            pickle.dump(save_dict, f)


    def load_state(self):
        """
        Loads the band state from a file

        Returns
        -------
        success : bool
            True if state was successfully loaded, False if not
        """
        import pickle

        try:
            with open(self.save_file, 'r') as f:
                d = pickle.load(f)
            old_md5 = d['_ms_file_md5']

            # Compare saved and current md5 value
            self.get_md5()
            if old_md5 != self._ms_file_md5:
                # If changes detected, saved state is not current
                return False
            else:
                # If no changes detected, restore state
                self.__dict__.update(d)
                return True
        except:
            return False

