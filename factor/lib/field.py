"""
Definition of the field class and a few related functions
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
import hashlib


class Field(object):
    """
    The Field object contains parameters needed for processing of the field

    Parameters
    ----------
    parset : dict
        Parset with processing parameters

    """
    def __init__(self, parset):

        # Set basic info
        self.name = 'field'
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.parset = parset.copy()
        self.working_dir = self.parset['dir_working']
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')

        # Get info on files
        self.files = self.parset['mss']
        self.numMS = len(self.files)

        # Load state (if any) and check files
        has_state = self.load_state()
        if not has_state:
            self.check_files()
            self.save_state()

        # Set up calibration parameters
        self.set_calibration_parameters()

        self.log.debug("Using {0} files.".format(len(self.files)))


    def check_files(self):
        """
        Checks input MS files and stores various parameters relating to them
        """
        # Make sure all files exist
        self.file_md5s = []
        for MS_id in xrange(self.numMS):
            if not os.path.exists(self.files[MS_id]):
                self.log.critical('MS file "{0}" not found!'.format(self.files[MS_id]))
                sys.exit(1)
            md5hash = hashlib.md5(open('{0}/OBSERVATION/table.f0'.format(self.files[MS_id]),
                'rb').read()).hexdigest()
            self.file_md5s.append(md5hash)

        # Get frequency and time info
        # Calculate time per sample and number of samples
        self.sumsamples = 0
        self.minSamplesPerFile = 4294967295  # If LOFAR lasts that many seconds then I buy you a beer.
        self.starttime = []
        self.endtime = []
        self.timepersample = []
        self.numsamples = []
        self.referencefreq = []
        self.startfreq = []
        self.endfreq = []
        self.numchannels = []
        self.channelwidth = []
        for MSid in xrange(self.numMS):
            tab = pt.table(self.files[MSid], ack=False)
            self.starttime.append(np.min(tab.getcol('TIME')))
            self.endtime.append(np.max(tab.getcol('TIME')))
            self.timepersample.append(tab.getcell('EXPOSURE',0))
            self.numsamples.append(len(tab.getcol('TIME')))
            self.sumsamples += numsamples
            self.minSamplesPerFile = min(self.minSamplesPerFile, numsamples)
            tab.close()
            sw = pt.table(self.files[MSid]+'::SPECTRAL_WINDOW', ack=False)
            self.referencefreq.append(sw.col('REF_FREQUENCY')[0])
            self.startfreq.append(np.min(sw.col('CHAN_FREQ')[0]))
            self.endfreq.append(np.max(sw.col('CHAN_FREQ')[0]))
            self.numchannels.append(sw.col('NUM_CHAN')[0])
            self.channelwidth.append(sw.col('CHAN_WIDTH')[0][0])
            sw.close()

        # Check that all MSs have the same frequency axis
        enforce_uniform_frequency_structure = False
        if enforce_uniform_frequency_structure:
            for MS_id in xrange(1, self.numMS):
                if self.numchannels[0] != self.numchannels[MS_id] \
                        or self.startfreq[0] != self.startfreq[MS_id] \
                        or self.endfreq[0] != self.endfreq[MS_id] \
                        or self.channelwidth[0] != self.channelwidth[MS_id]:
                    self.log.critical('Frequency axis for MS {0} differs from the one for MS {1}! '
                                      'Exiting!'.format(self.files[MS_id], self.files[0]))
                    sys.exit(1)

        # Get pointing info
        obs = pt.table(self.files[0]+'::FIELD', ack=False)
        self.ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
        if self.ra < 0.:
            self.ra = 360.0 + (self.ra)
        self.dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
        obs.close()

        # Check that all MSs have the same pointing
        for MS_id in xrange(1, self.numMS):
            obs = pt.table(self.files[MS_id]+'::FIELD', ack=False)
            ra = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][0]))
            if ra < 0.:
                ra = 360.0 + (ra)
            dec = np.degrees(float(obs.col('REFERENCE_DIR')[0][0][1]))
            obs.close()
            if self.ra != ra or self.dec != dec:
                self.log.critical('Pointing for MS {0} differs from the one for MS {1}! '
                                  'Exiting!'.format(self.files[MS_id], self.files[0]))
                sys.exit(1)

        # Get the station diameter
        ant = pt.table(self.files[0]+'::ANTENNA', ack=False)
        self.diam = float(ant.col('DISH_DIAMETER')[0])
        ant.close()

        # Check that all MSs have the station diameter
        for MS_id in xrange(1, self.numMS):
            ant = pt.table(self.files[MS_id]+'::ANTENNA', ack=False)
            diam = float(ant.col('DISH_DIAMETER')[0])
            ant.close()
            if self.diam != diam:
                self.log.critical('Station diameter for MS {0} differs from the one for MS {1}! '
                                  'Exiting!'.format(self.files[MS_id], self.files[0]))
                sys.exit(1)

        # Find mean elevation and FOV
        for MS_id in xrange(self.numMS):
            # Add (virtual) elevation column to MS
            tab = pt.table(self.files[MS_id], ack=False)
            exiting_colnames = tab.colnames()
            if 'AZEL1' not in exiting_colnames:
                tab.close()
                pt.addDerivedMSCal(self.files[MS_id])
                tab = pt.table(self.files[MS_id], ack=False)

            # Calculate mean elevation
            if MS_id == 0:
                global_el_values = tab.getcol('AZEL1', rowincr=10000)[:, 1]
            else:
                global_el_values = np.hstack( (global_el_values, tab.getcol('AZEL1', rowincr=10000)[:, 1]) )
            tab.close()

            # Remove (virtual) elevation column from MS
            pt.removeDerivedMSCal(self.files[MS_id])
        self.mean_el_rad = np.mean(global_el_values)
        sec_el = 1.0 / np.sin(self.mean_el_rad)
        self.fwhm_deg = 1.1 * ((3.0e8 / np.mean(self.referencefreq)) / self.diam) * 180. / np.pi * sec_el


    def set_calibration_parameters(self):
        """
        Sets the calibration parameters per time chunk

        Parameters
        ----------
        chunksize : int
            Desired size of chunks in sec

        """
        target_chunksize = self.parset['chunk_size_sec']
        target_fast_timestep = self.parset['fast_timestep_sec']
        target_fast_freqstep = self.parset['fast_freqstep_hz']
        target_slow_timestep = self.parset['slow_timestep_sec']
        target_slow_freqstep = self.parset['slow_freqstep_hz']

        self.calibration_parameters = []
        for MS_id in xrange(self.numMS):
            # Create empty dict to hold parameters for this MS file
            self.calibration_parameters.append[{}]

            # Calculate time ranges of calibration chunks
            timepersample = self.timepersample[MS_id]
            numsamples = self.numsamples[MS_id]
            samplesperchunk = int(round(target_chunksize / timepersample))
            chunksize = samplesperchunk * timepersample
            mystarttime = self.starttime[MS_id]
            myendtime = self.endtime[MS_id]
            if (myendtime-mystarttime) > (2.*chunksize):
                nchunks = int((numsamples*timepersample)/chunksize)
            else:
                nchunks = 1
            self.log.debug('Using {1} calibration chunk(s) for {0}'.format(self.files[MS_id], nchunks))
            self.calibration_parameters[MS_id]['starttime'] = [mystarttime+(chunksize*i) for i in range(nchunks)]
            self.calibration_parameters[MS_id]['ntimes'] = [samplesperchunk*i for i in range(nchunks)]
            self.calibration_parameters[MS_id]['ntimes'][-1] = 0 # set last entry to extend until end

            # Set solution intervals per chunk
            freqpersample = self.chan_width_hz
            self.calibration_parameters[MS_id]['solint_time_fast'] = [int(round(target_fast_timestep / timepersample))] * nchunks
            self.calibration_parameters[MS_id]['solint_freq_fast'] = [int(round(target_fast_freqstep / freqpersample))] * nchunks
            self.calibration_parameters[MS_id]['solint_time_slow'] = [int(round(target_slow_timestep / timepersample))] * nchunks
            self.calibration_parameters[MS_id]['solint_freq_slow'] = [int(round(target_slow_freqstep / freqpersample))] * nchunks


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
            for step in range(self.nchan,0,-1):
                if (self.nchan % step) == 0:
                    tmp_divisors.append(step)
            self.freq_divisors = np.array(tmp_divisors)
        idx = np.argmin(np.abs(self.freq_divisors-freqstep))
        return self.freq_divisors[idx]


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
            old_md5s = self.file_md5s[:]
            with open(self.save_file, 'r') as f:
                d = pickle.load(f)
            self.__dict__.update(d)

            # Compare old and new file md5 values. If they differ, return False
            if old_md5s != self.file_md5s:
                return False
            else:
                return True
        except:
            return False
