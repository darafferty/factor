"""
Definition of the Observation class that holds parameters for each measurement
set
"""
import os
import logging
import casacore.tables as pt
import numpy as np
from astropy.time import Time


class Observation(object):
    """
    The Observation object contains various MS-related parameters

    Parameters
    ----------
    ms_filename : str
        Filename of the MS file
    """
    def __init__(self, ms_filename):
        self.ms_filename = ms_filename
        self.name = os.path.basename(self.ms_filename)
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.scan_ms()

    def scan_ms(self):
        """
        Scans input MS and stores info
        """
        # Get time info
        tab = pt.table(self.ms_filename, ack=False)
        self.starttime = np.min(tab.getcol('TIME'))
        self.endtime = np.max(tab.getcol('TIME'))
        self.timepersample = tab.getcell('EXPOSURE', 0)
        self.numsamples = int(np.ceil((self.endtime - self.starttime) / self.timepersample))
        tab.close()

        # Get frequency info
        sw = pt.table(self.ms_filename+'::SPECTRAL_WINDOW', ack=False)
        self.referencefreq = sw.col('REF_FREQUENCY')[0]
        self.startfreq = np.min(sw.col('CHAN_FREQ')[0])
        self.endfreq = np.max(sw.col('CHAN_FREQ')[0])
        self.numchannels = sw.col('NUM_CHAN')[0]
        self.channelwidth = sw.col('CHAN_WIDTH')[0][0]
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
        Sets the calibration parameters

        Parameters
        ----------
        parset : dict
            Parset with processing parameters
        """
        # Get various target parameters
        target_time_chunksize = parset['chunk_size_sec']
        target_freq_chunksize = parset['chunk_size_hz']
        target_fast_timestep = parset['calibration_specific']['fast_timestep_sec']
        target_fast_freqstep = parset['calibration_specific']['fast_freqstep_hz']
        target_slow_timestep = parset['calibration_specific']['slow_timestep_sec']
        target_slow_freqstep = parset['calibration_specific']['slow_freqstep_hz']

        # Calculate time ranges of calibration chunks
        timepersample = self.timepersample
        numsamples = self.numsamples
        samplesperchunk = int(round(target_time_chunksize / timepersample))
        chunksize = samplesperchunk * timepersample
        mystarttime = self.starttime
        myendtime = self.endtime
        if (myendtime-mystarttime) > (2.0*chunksize):
            nchunks = int((numsamples*timepersample)/chunksize)
        else:
            nchunks = 1
        self.ntimechunks = nchunks
        self.log.debug('Using {} time chunk(s) for fast-phase calibration'.format(self.ntimechunks))
        self.calibration_parameters = {}
        self.calibration_parameters['timechunk_filename'] = [self.ms_filename] * self.ntimechunks
        starttimes = [mystarttime+(chunksize * i) for i in range(self.ntimechunks)]
        self.calibration_parameters['starttime'] = [self.convert_mjd(t) for t in starttimes]
        self.calibration_parameters['ntimes'] = [samplesperchunk] * self.ntimechunks
        self.calibration_parameters['ntimes'][-1] = 0  # set last entry to extend until end

        # Calculate frequency ranges of calibration chunks
        channelwidth = self.channelwidth
        numchannels = self.numchannels
        channelsperchunk = int(round(target_freq_chunksize / channelwidth))
        chunksize = channelsperchunk * channelwidth
        mystartfreq = self.startfreq
        myendfreq = self.endfreq
        if (myendfreq-mystartfreq) > (2.0*chunksize):
            nchunks = int((numchannels*channelwidth)/chunksize)
        else:
            nchunks = 1
        self.nfreqchunks
        self.log.debug('Using {1} frequency chunk(s) for slow-gain calibration'.format(self.nfreqchunks))
        self.calibration_parameters['freqchunk_filename'] = [self.ms_filename] * self.nfreqchunks
        self.calibration_parameters['startchan'] = [channelsperchunk * i for i in range(nchunks)]
        self.calibration_parameters['nchan'] = [channelsperchunk] * nchunks
        self.calibration_parameters['nchan'][-1] = 0  # set last entry to extend until end

        # Set solution intervals (same for every calibration chunk)
        freqpersample = self.channelwidth
        self.calibration_parameters['solint_fast_timestep'] = max(1, [int(round(target_fast_timestep /
                                                                      timepersample))] * self.nchunks)
        self.calibration_parameters['solint_fast_freqstep'] = max(1, [self.get_nearest_frequstep(target_fast_freqstep /
                                                                      freqpersample)] * self.nchunks)
        self.calibration_parameters['solint_slow_timestep'] = max(1, [int(round(target_slow_timestep /
                                                                      timepersample))] * self.nchunks)
        self.calibration_parameters['solint_slow_freqstep'] = max(1, [self.get_nearest_frequstep(target_slow_freqstep /
                                                                      freqpersample)] * self.nchunks)

    def convert_mjd(self, mjd_sec):
        """
        Converts MJD to casacore MVTime

        Parameters
        ----------
        mjd_sec : float
            MJD time in seconds

        Returns
        -------
        mvtime : str
            Casacore MVTime string
        """
        t = Time(mjd_sec / 3600 / 24, format='mjd', scale='utc')
        date, hour = t.iso.split(' ')
        year, month, day = date.split('-')
        d = t.datetime
        month = d.ctime().split(' ')[1]

        return '{0}{1}{2}/{3}'.format(day, month, year, hour)

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
        # Generate a list of possible values for freqstep
        if not hasattr(self, 'freq_divisors'):
            tmp_divisors = []
            for step in range(self.numchannels, 0, -1):
                if (self.numchannels % step) == 0:
                    tmp_divisors.append(step)
            self.freq_divisors = np.array(tmp_divisors)

        # Find nearest
        idx = np.argmin(np.abs(self.freq_divisors - freqstep))

        return self.freq_divisors[idx]
