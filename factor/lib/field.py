"""
Definition of the Field class
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
    The Field object stores parameters needed for processing of the field

    Parameters
    ----------
    parset : dict
        Parset with processing parameters

    """
    def __init__(self, parset):

        # Initialize basic info
        self.name = 'field'
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.parset = parset.copy()
        self.working_dir = self.parset['dir_working']
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.ms_filenames = self.parset['mss']
        self.numMS = len(self.ms_files)
        self.patch_target_flux_jy = None

        # Scan MS files to get observation info
        self.scan_observations()

        # Load initial sky model and group
        self.skymodel = lsmtool.load(parset['initial_skymodel'])
        self.patch_target_flux_jy = parset[direction_specific]['patch_target_flux_jy']
        self.group_skymodel()

        # Set calibration parameters
        self.set_calibration_parameters()

        # Save state
        self.save_state()


    def scan_observations(self):
        """
        Checks input MS and stores the associated Observation objects
        """
        self.observations = []
        for ms_filename in self.ms_filenames:
            self.observations.append(Observation(ms_filename, self.working_dir))

        # Check that all observations have the same frequency axis
        # NOTE: this may not be necessary and is disable for now
        obs0 = self.observations[0]
        enforce_uniform_frequency_structure = False
        if enforce_uniform_frequency_structure:
            for obs in self.observations:
                if (obs0.numchannels != obs.numchannels
                    or obs0.startfreq != obs.startfreq
                    or obs0.endfreq != obs.endfreq
                    or obs0.channelwidth != obs.channelwidth:
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
        sec_el = 1.0 / np.sin(np.mean(self.mean_el_rad))
        self.fwhm_deg = 1.1 * ((3.0e8 / np.mean(ref_freq_list)) / self.diam) * 180. / np.pi * sec_el


    def group_skymodel(self):
        """
        Groups the sky model into groups of target flux
        """
        flux = self.patch_target_flux_jy
    	self.skymodel.group('threshold', FWHM='60.0 arcsec’)
	    self.skymodel.group(algorithm='tessellate', targetFlux=flux, method=‘mid’, byPatch=True)


    def save_state(self):
        """
        Saves the state to a file

        """
        import pickle

        with open(self.save_file, 'wb') as f:
            # Remove log object, as it cannot be pickled
            save_dict = self.__dict__.copy()
            save_dict.pop('log')
            pickle.dump(save_dict, f)


    def load_state(self):
        """
        Loads the state from a file

        Returns
        -------
        success : bool
            True if state was successfully loaded, False if not
        """
        import pickle

        try:
            with open(self.save_file, 'r') as f:
                d = pickle.load(f)
            self.__dict__.update(d)
            return True
        except:
            return False
