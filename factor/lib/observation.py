"""
Definition of the Observation class that holds parameters for each measurement
set
"""
import os
import logging
import casacore.tables as pt
import numpy as np
from astropy.time import Time
from factor.directions import get_target_timewidth, get_target_bandwidth
from factor.cluster import get_time_chunksize, get_frequency_chunksize

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
        if 'HBA' in ant.col('NAME')[0]:
            self.antenna = 'HBA'
        elif 'LBA' in ant.col('NAME')[0]:
            self.antenna = 'LBA'
        else:
            self.log.warning('Antenna type not recognized (only LBA and HBA data '
                             'are supported at this time)')
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
        # Get the target solution intervals
        target_fast_timestep = parset['calibration_specific']['fast_timestep_sec']
        target_fast_freqstep = parset['calibration_specific']['fast_freqstep_hz']
        target_slow_timestep = parset['calibration_specific']['slow_timestep_sec']
        target_slow_freqstep = parset['calibration_specific']['slow_freqstep_hz']

        # Find solution intervals for fast-phase solve
        timepersample = self.timepersample
        channelwidth = self.channelwidth
        solint_fast_timestep = max(1, int(round(target_fast_timestep / timepersample)))
        solint_fast_freqstep = max(1, self.get_nearest_frequstep(target_fast_freqstep / channelwidth))

        # Calculate time ranges of calibration chunks for fast-phase solve. Try
        # to ensure that the number of samples per chunk is an even multiple of
        # the solution interval
        if parset['chunk_size_sec'] is None:
            target_time_chunksize = get_time_chunksize(parset['cluster_specific'], self.timepersample,
                                                       self.numsamples, solint_fast_timestep)
        else:
            target_time_chunksize = parset['chunk_size_sec']
        samplesperchunk = int(round(target_time_chunksize / timepersample))
        chunksize = samplesperchunk * timepersample
        mystarttime = self.starttime
        myendtime = self.endtime
        if (myendtime-mystarttime) > (2*chunksize):
            nchunks = int(round(float(self.numsamples)*timepersample)/chunksize)
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

        # Find solution intervals for slow-gain solve
        solint_slow_timestep = max(1, int(round(target_slow_timestep / timepersample)))
        solint_slow_freqstep = max(1, self.get_nearest_frequstep(target_slow_freqstep / channelwidth))

        # Calculate frequency ranges of calibration chunks for slow-gain solve. Try
        # to ensure that the number of samples per chunk is an even multiple of
        # the solution interval
        numchannels = self.numchannels
        if parset['chunk_size_hz'] is None:
            target_freq_chunksize = get_freq_chunksize(parset['cluster_specific'], channelwidth,
                                                       solint_slow_freqstep, solint_slow_timestep)
        else:
            target_freq_chunksize = parset['chunk_size_hz']
        chunksize = channelsperchunk * channelwidth
        mystartfreq = self.startfreq
        myendfreq = self.endfreq
        if (myendfreq-mystartfreq) > (2.0*chunksize):
            nchunks = int((numchannels*channelwidth)/chunksize)
        else:
            nchunks = 1
        self.nfreqchunks = nchunks
        self.log.debug('Using {} frequency chunk(s) for slow-gain calibration'.format(self.nfreqchunks))
        self.calibration_parameters['freqchunk_filename'] = [self.ms_filename] * self.nfreqchunks
        self.calibration_parameters['startchan'] = [channelsperchunk * i for i in range(nchunks)]
        self.calibration_parameters['nchan'] = [channelsperchunk] * nchunks
        self.calibration_parameters['nchan'][-1] = 0  # set last entry to extend until end

        # Set solution intervals (same for every calibration chunk)
        self.calibration_parameters['solint_fast_timestep'] = [solint_fast_timestep] * self.ntimechunks
        self.calibration_parameters['solint_fast_freqstep'] = [solint_fast_freqstep] * self.ntimechunks
        self.calibration_parameters['solint_slow_timestep'] = [solint_slow_timestep] * self.nfreqchunks
        self.calibration_parameters['solint_slow_freqstep'] = [solint_slow_freqstep] * self.nfreqchunks

    def set_imaging_parameters(self, cellsize_arcsec, max_peak_smearing,
                               width_ra, width_dec):
        """
        Sets the imaging parameters

        Parameters
        ----------
        field : Field object
            Field object
        cellsize_arcsec : float
            Pixel size in arcsec for imaging
        width_ra : float
            Width in RA of image in degrees
        width_dec : float
            Width in Dec of image in degrees
        """
        mean_freq_mhz = self.referencefreq
        peak_smearing_factor = np.sqrt(1.0 - max_peak_smearing)
        chan_width_hz = self.channelwidth
        nchan = self.numchannels
        timestep_sec = self.timepersample

        # Get target time and frequency averaging steps
        delta_theta_deg = max(width_ra, width_dec) / 2.0
        resolution_deg = 3.0 * cellsize_arcsec / 3600.0  # assume normal sampling of restoring beam
        target_timewidth_sec = min(120.0, get_target_timewidth(delta_theta_deg,
                                   resolution_deg, peak_smearing_factor))
        target_bandwidth_mhz = min(2.0, get_target_bandwidth(mean_freq_mhz,
                                   delta_theta_deg, resolution_deg, peak_smearing_factor))
        self.log.debug('Target timewidth for imaging is {} s'.format(target_timewidth_sec))
        self.log.debug('Target bandwidth for imaging is {} MHz'.format(target_bandwidth_mhz))

        # Find averaging steps for above target values
        self.imaging_parameters = {}
        image_freqstep = max(1, min(int(round(target_bandwidth_mhz * 1e6 / chan_width_hz)), nchan))
        self.imaging_parameters['image_freqstep'] = self.get_nearest_frequstep(image_freqstep)
        self.imaging_parameters['image_timestep'] = max(1, int(round(target_timewidth_sec / timestep_sec)))
        self.log.debug('Using averaging steps of {0} channels and {1} time slots '
                       'for imaging'.format(self.imaging_parameters['image_freqstep'],
                                            self.imaging_parameters['image_timestep']))

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
