"""
Definition of the Field class
"""
import os
import sys
import logging
import numpy as np
import lsmtool
from factor.lib.observation import Observation


class Field(object):
    """
    The Field object stores parameters needed for processing of the field

    Parameters
    ----------
    parset : dict
        Parset with processing parameters
    """
    def __init__(self, parset):
        # Save basic info
        self.name = 'field'
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.parset = parset.copy()
        self.working_dir = self.parset['dir_working']
        self.ms_filenames = self.parset['mss']
        self.numMS = len(self.ms_filenames)
        self.data_colname = 'DATA'
        self.solve_min_uv_lambda = self.parset['calibration_specific']['solve_min_uv_lambda']
        self.solve_tecandphase = self.parset['calibration_specific']['solve_tecandphase']

        # Scan MS files to get observation info
        self.scan_observations()

        # Make calibration sky model by grouping the initial sky model
        self.make_calibration_skymodel(self.parset['initial_skymodel'])

    def scan_observations(self):
        """
        Checks input MS and stores the associated Observation objects
        """
        self.observations = []
        for ms_filename in self.ms_filenames:
            self.observations.append(Observation(ms_filename))

        # Check that all observations have the same frequency axis
        # NOTE: this may not be necessary and is disabled for now
        obs0 = self.observations[0]
        enforce_uniform_frequency_structure = False
        if enforce_uniform_frequency_structure:
            for obs in self.observations:
                if (obs0.numchannels != obs.numchannels or
                        obs0.startfreq != obs.startfreq or
                        obs0.endfreq != obs.endfreq or
                        obs0.channelwidth != obs.channelwidth):
                    self.log.critical('Frequency axis for MS {0} differs from the one for MS {1}! '
                                      'Exiting!'.format(self.obs.ms_filename, self.obs0.ms_filename))
                    sys.exit(1)

        # Check that all MSs have the same pointing
        self.ra = obs0.ra
        self.dec = obs0.dec
        for obs in self.observations:
            if self.ra != obs.ra or self.dec != obs.dec:
                self.log.critical('Pointing for MS {0} differs from the one for MS {1}! '
                                  'Exiting!'.format(self.obs.ms_filename, self.obs0.ms_filename))
                sys.exit(1)

        # Check that all observations have the same station diameter
        self.diam = obs0.diam
        for obs in self.observations:
            if self.diam != obs.diam:
                self.log.critical('Station diameter for MS {0} differs from the one for MS {1}! '
                                  'Exiting!'.format(self.obs.ms_filename, self.obs0.ms_filename))
                sys.exit(1)

        # Find mean elevation and FOV over all observations
        el_rad_list = []
        ref_freq_list = []
        for obs in self.observations:
            el_rad_list.append(obs.mean_el_rad)
            ref_freq_list.append(obs.referencefreq)
        sec_el = 1.0 / np.sin(np.mean(el_rad_list))
        self.mean_el_rad = np.mean(el_rad_list)
        self.fwhm_deg = 1.1 * ((3.0e8 / np.mean(ref_freq_list)) /
                               self.diam) * 180. / np.pi * sec_el

        # Set calibration parameters for each observation
        ntimechunks = 0
        nfreqchunks = 0
        for obs in self.observations:
            obs.set_calibration_parameters(self.parset)
            ntimechunks += obs.ntimechunks
            nfreqchunks += obs.nfreqchunks
        self.ntimechunks = ntimechunks
        self.nfreqchunks = nfreqchunks

    def get_obs_parameters(self, parameter):
        """
        Returns list of calibration parameters for all observations

        Parameters
        ----------
        parameter : str
            Name of calibration parameter to return

        Returns
        -------
        parameters : list
            List of parameters, with one entry for each time or frequency chunk
            of each observation
        """
        return sum([obs.calibration_parameters[parameter] for obs in self.observations], [])

    def make_calibration_skymodel(self, skymodel_filename):
        """
        Groups the sky model into groups of target flux and saves to disk

        Parameters
        ----------
        skymodel_filename : str
            Filename of input makesourcedb sky model file
        """
        skymodel = lsmtool.load(skymodel_filename)
        flux = self.parset['direction_specific']['patch_target_flux_jy']
        self.log.info('Grouping sky model to form calibration patches...')
        skymodel.group('threshold', FWHM='60.0 arcsec')
        skymodel.group(algorithm='tessellate', targetFlux=flux, method='mid', byPatch=True)

        self.skymodel_file = os.path.join(self.working_dir, 'skymodels', 'calibration_skymodel.txt')
        skymodel.write(self.skymodel_file, clobber=True)
        self.log.info('Using {0} calibration patches of ~ {1} Jy each'.format(len(skymodel), flux))
