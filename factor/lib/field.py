"""
Definition of the Field class
"""
import os
import sys
import logging
import numpy as np
import lsmtool
from factor.lib.observation import Observation
from factor.lib.sector import Sector


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

        # Set up imaging sectors
        self.makeWCS()
        self.define_sectors()

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
        self.fwhm_ra_deg = self.fwhm_deg / sec_el
        self.fwhm_dec_deg = self.fwhm_deg

        # Set calibration parameters for each observation
        ntimechunks = 0
        nfreqchunks = 0
        for obs in self.observations:
            obs.set_calibration_parameters(self.parset)
            ntimechunks += obs.ntimechunks
            nfreqchunks += obs.nfreqchunks
        self.ntimechunks = ntimechunks
        self.nfreqchunks = nfreqchunks

    def get_calibration_parameters(self, parameter):
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
        self.log.info('Using {0} calibration patches of ~ {1} Jy each'.format(len(skymodel.getPatchNames()), flux))

    def define_sectors(self):
        """
        Defines the imaging sectors
        """
        nsectors_ra = self.parset['imaging_specific']['num_sectors']
        if nsectors_ra == 1 and len(self.parset['cluster_specific']['node_list']) == 1:
            nsectors_dec = 1
            width_ra = self.fwhm_ra_deg / nsectors_ra
            width_dec = self.fwhm_dec_deg / nsectors_dec
            center_x, center_y = self.radec2xy([self.ra], [self.dec])
            x = np.array([[center_x]])
            y = np.array([[center_y]])
        else:
            nsectors_dec = int(np.ceil(nsectors_ra / np.sin(self.mean_el_rad)))
            width_ra = self.fwhm_ra_deg / nsectors_ra
            width_dec = self.fwhm_dec_deg / nsectors_dec
            width_x = width_ra / abs(self.wcs.wcs.cdelt[0])
            width_y = width_dec / abs(self.wcs.wcs.cdelt[1])
            center_x, center_y = self.radec2xy([self.ra], [self.dec])
            min_x = center_x - nsectors_ra / 2.0 * width_x
            max_x = center_x + nsectors_ra / 2.0 * width_x
            min_y = center_y - nsectors_dec / 2.0 * width_y
            max_y = center_y + nsectors_dec / 2.0 * width_y
            x = np.linspace(min_x, max_x, width_x)
            y = np.linspace(min_y, max_y, width_y)
            x, y = np.meshgrid(x, y)

        # Define sectors
        self.sectors = []
        for i in range(nsectors_ra):
            for i in range(nsectors_dec):
                name = 'sector_{0}_{1}'.format(i, j)
                ra, dec = self.xy2radec([x[i, j]], [y[i, j]])
                self.sectors.append(Sector(name, ra, dec, width_ra, width_dec, field))

        for this_sector in self.sectors:
            # For each sector, check for intersection with other sectors
            other_sectors = self.sectors[:].remove(this_sector)
            for other_sector in other_sectors:
                if this_sector.poly.contains(other_sector.poly.centroid):
                    # If centroid is outside, difference the polys
                    this_sector.poly = this_sector.poly.difference(other_sector.poly)
                else:
                    # If point is inside, union the polys
                    this_sector.poly = this_sector.poly.union(other_sector.poly)

            # Make sector region and vertices files
            this_sector.make_vertices_file()
            this_sector.make_region_file(os.path.join(self.working_dir, 'regions',
                                                      '{}_region.txt'.format(self.name)))

    def get_imaging_parameters(self, parameter):
        """
        Returns list of imaging parameters for all sectors

        Parameters
        ----------
        parameter : str
            Name of parameter to return

        Returns
        -------
        parameters : list
            List of parameters, with one entry for each observation
        """
        return sum([sector.get_imaging_parameters(parameter) for sector in self.sectors], [])

    def radec2xy(self, RA, Dec):
        """
        Returns x, y for input RA, Dec

        Parameters
        ----------
        RA : list
            List of RA values in degrees
        Dec : list
            List of Dec values in degrees

        Returns
        -------
        x, y : list, list
            Lists of x and y pixel values corresponding to the input RA and Dec
            values
        """
        x = []
        y = []

        for ra_deg, dec_deg in zip(RA, Dec):
            ra_dec = np.array([[ra_deg, dec_deg]])
            x.append(self.wcs.wcs_world2pix(ra_dec, 0)[0][0])
            y.append(self.wcs.wcs_world2pix(ra_dec, 0)[0][1])

        return x, y

    def xy2radec(self, x, y):
        """
        Returns input RA, Dec for input x, y

        Parameters
        ----------
        x : list
            List of x values in pixels
        y : list
            List of y values in pixels

        Returns
        -------
        RA, Dec : list, list
            Lists of RA and Dec values corresponding to the input x and y pixel
            values
        """
        RA = []
        Dec = []

        for xp, yp in zip(x, y):
            x_y = np.array([[xp, yp]])
            RA.append(self.wcs.wcs_pix2world(x_y, 0)[0][0])
            Dec.append(self.wcs.wcs_pix2world(x_y, 0)[0][1])

        return RA, Dec

    def makeWCS(self):
        """
        Makes simple WCS object

        Returns
        -------
        w : astropy.wcs.WCS object
            A simple TAN-projection WCS object for specified reference position

        """
        from astropy.wcs import WCS

        wcs_pixel_scale = 0.066667 # degrees/pixel
        w = WCS(naxis=2)
        w.wcs.crpix = [1000, 1000]
        w.wcs.cdelt = np.array([-wcs_pixel_scale, wcs_pixel_scale])
        w.wcs.crval = [self.ra, self.dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.set_pv([(2, 1, 45.0)])
        self.wcs = w

