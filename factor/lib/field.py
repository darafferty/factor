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
from factor.scripts import blank_image, mosaic_images


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
        self.flag_abstime = self.parset['flag_abstime']
        self.flag_baseline = self.parset['flag_baseline']
        self.flag_freqrange = self.parset['flag_freqrange']
        self.flag_expr = self.parset['flag_expr']

        # Scan MS files to get observation info
        self.scan_observations()

        # Make calibration and source sky models by grouping the initial sky model
        self.make_skymodels(self.parset['initial_skymodel'], self.parset['regroup_initial_skymodel'])

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

    def make_source_skymodel(self, skymodel_filename, regroup=True):
        """
        Groups the sky model into sources

        Parameters
        ----------
        skymodel_filename : str
            Filename of input makesourcedb sky model file
        """
        skymodel = lsmtool.load(skymodel_filename)

        # Group by thresholding
        skymodel.group('threshold', FWHM='60.0 arcsec')
        self.source_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'source_skymodel.txt')
        skymodel.write(self.source_skymodel_file, clobber=True)

    def make_skymodels(self, skymodel, regroup=True):
        """
        Groups a sky model into source and calibration patches

        Parameters
        ----------
        skymodel : str or LSMTool skymodel object
            Filename of input makesourcedb sky model file
        regroup : bool, optional
            If False, the calibration sky model is not regrouped to the target flux.
            Instead, the existing calibration groups are used
        """
        if type(skymodel) is str:
            skymodel = lsmtool.load(skymodel)

        if regroup:
            flux = self.parset['direction_specific']['patch_target_flux_jy']
            self.log.info('Grouping sky model to form calibration patches...')
            source_skymodel = skymodel
        else:
            source_skymodel = skymodel.copy()

        # Group by thresholding and write out the source sky model
        source_skymodel.group('threshold', FWHM='60.0 arcsec')
        self.source_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'source_skymodel.txt')
        source_skymodel.write(self.source_skymodel_file, clobber=True)

        # Now tesselate to get patches of the target flux and write out calibration sky model
        if regroup:
            self.log.info('Using {0} calibration patches of ~ {1} Jy each'.format(len(skymodel.getPatchNames()), flux))
            skymodel.group(algorithm='tessellate', targetFlux=flux, method='mid', byPatch=True)
        self.calibration_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'calibration_skymodel.txt')
        skymodel.write(self.calibration_skymodel_file, clobber=True)

    def update_skymodels(self, iter):
        """
        Updates the source and calibration sky models from the output sector sky model(s)
        """
        sector_skymodels = [sector.get_output_skymodel_filename() for sector in self.sectors]
        skymodel = lsmtool.load(sector_skymodels[0])
        sector_skymodels.pop(0)
        for s2 in sector_skymodels:
            skymodel.concatenate(s2)
        skymodel.write(self.calibration_skymodel_file, clobber=True)
        self.make_skymodels(skymodel)

    def define_sectors(self):
        """
        Defines the imaging sectors
        """
        nsectors_ra = self.parset['imaging_specific']['nsectors_per_side']
        if nsectors_ra == 1 and len(self.parset['cluster_specific']['node_list']) == 1:
            # For a single machine, use a single sector
            nsectors_dec = 1
            width_ra = self.fwhm_ra_deg
            width_dec = self.fwhm_dec_deg
            center_x, center_y = self.radec2xy([self.ra], [self.dec])
            x = np.array([center_x])
            y = np.array([center_y])
        else:
            nsectors_dec = int(np.ceil(nsectors_ra / np.sin(self.mean_el_rad)))
            width_ra = self.fwhm_ra_deg / nsectors_ra
            width_dec = self.fwhm_dec_deg / nsectors_dec
            width_x = width_ra / abs(self.wcs.wcs.cdelt[0])
            width_y = width_dec / abs(self.wcs.wcs.cdelt[1])
            center_x, center_y = self.radec2xy([self.ra], [self.dec])
            min_x = center_x - width_x / 2.0 * (nsectors_ra - 1)
            max_x = center_x + width_x / 2.0 * (nsectors_ra - 1)
            min_y = center_y - width_y / 2.0 * (nsectors_dec - 1)
            max_y = center_y + width_y / 2.0 * (nsectors_dec - 1)
            x = np.linspace(min_x, max_x, nsectors_ra)
            y = np.linspace(min_y, max_y, nsectors_dec)
            x, y = np.meshgrid(x, y)
            self.log.info('Using {0} imaging sectors ({1} in RA, {2} in Dec)'.format(
                          nsectors_ra*nsectors_dec, nsectors_ra, nsectors_dec))

        # Define sectors
        self.sectors = []
        n = 0
        for i in range(nsectors_ra):
            for j in range(nsectors_dec):
                if nsectors_ra == 1 and nsectors_dec == 1:
                    name = 'field'
                else:
                    name = 'sector_{0}'.format(n)
                n += 1
                ra, dec = self.xy2radec([x[j, i]], [y[j, i]])
                self.sectors.append(Sector(name, ra[0], dec[0], width_ra, width_dec, self))

        for this_sector in self.sectors:
            # For each sector, check for intersection with other sectors
            for other_sector in self.sectors:
                if this_sector is not other_sector and this_sector.poly.intersects(other_sector.poly):
                    this_sector.poly = this_sector.poly.difference(other_sector.poly)

            # Make sector region and vertices files
            this_sector.make_vertices_file()
            this_sector.make_region_file(os.path.join(self.working_dir, 'regions',
                                                      '{}_region_ds9.reg'.format(this_sector.name)))

            # Set the imaging parameters for selfcal
            this_sector.set_imaging_parameters(self.parset['imaging_specific']['selfcal_cellsize_arcsec'],
                                               self.parset['imaging_specific']['selfcal_robust'],
                                               0.0,
                                               self.parset['imaging_specific']['selfcal_min_uv_lambda'],
                                               None,
                                               0.15,
                                               self.parset['imaging_specific']['wsclean_bl_averaging'],
                                               self.parset['imaging_specific']['selfcal_multiscale_scales_pixel'],
                                               self.parset['imaging_specific']['use_idg'],
                                               self.parset['imaging_specific']['idg_mode'])

            # Transfer flagging parameters so they are also used during imaging
            this_sector.flag_abstime = self.flag_abstime
            this_sector.flag_baseline = self.flag_baseline
            this_sector.flag_freqrange = self.flag_freqrange
            this_sector.flag_expr = self.flag_expr

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

    def check_selfcal_convergence(self):
        """
        Checks whether selfcal has converged or not on a sector-by-sector basis

        Returns
        -------
        result : bool
            True if all sectors have converged, False if not
        """
        return False

    def make_mosaic(self, iter, image_id=None):
        """
        Make mosaic of the sector images

        Parameters
        ----------
        iter : int
            Iteration index
        image_id : str, optional
            Imaging ID
        """
        if len(self.sectors) > 1:
            # Blank the sector images
            blanked_images = []
            for sector in self.sectors:
                input_image_file = sector.get_output_image_filename(image_id)
                vertices_file = sector.vertices_file
                output_image_file = input_image_file * '_blanked'
                blanked_images.append(output_image_file)
                blank_image.main(input_image_file, vertices_file, output_image_file)

            # Make the mosaic
            outfile = os.path.join(self.parset['dir_working'], 'results', 'image_{}'.format(iter),
                                   'field', 'field_mosaic.fits')
            self.output_image_filename = outfile
            mosaic_images.main(blanked_images, outfile)
        else:
            self.output_image_filename = self.sectors[0].get_output_image_filename(image_id)
