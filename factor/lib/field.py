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
from lofarpipe.support.utilities import create_directory
from shapely.geometry import Point, Polygon
import rtree

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
        self.log.debug('Scanning observations...')
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
        self.log.info('Reading sky model...')
        if type(skymodel) is str:
            skymodel = lsmtool.load(skymodel)
        source_skymodel = skymodel.copy()

        # Group by thresholding and write out the source sky model (only used when there
        # is more than one imaging sector)
        if self.parset['imaging_specific']['nsectors_ra'] > 0:
            self.log.info('Identifying sources...')
            source_skymodel.group('threshold', FWHM='60.0 arcsec')
            source_skymodel_filt = source_skymodel.copy()
            source_skymodel_filt.remove('Patch = patch_*', force=True) # Remove sources that did not threshold
            self.source_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'source_skymodel.txt')
            source_skymodel.write(self.source_skymodel_file, clobber=True)
            self.source_skymodel = source_skymodel_filt
        if regroup:
            flux = self.parset['direction_specific']['patch_target_flux_jy']
            self.log.info('Grouping sky model to form calibration patches...')
            calibration_skymodel = source_skymodel
        else:
            calibration_skymodel = skymodel

        # Now tesselate to get patches of the target flux and write out calibration sky model
        if regroup:
            self.log.info('Grouping sky model to form calibration patches of ~ {} Jy each...'.format(len(flux)))
            calibration_skymodel.group(algorithm='tessellate', targetFlux=flux, method='mid', byPatch=True)
        self.log.info('Using {} calibration patches'.format(len(calibration_skymodel.getPatchNames())))
        self.calibration_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'calibration_skymodel.txt')
        calibration_skymodel.write(self.calibration_skymodel_file, clobber=True)
        self.calibration_skymodel = calibration_skymodel

    def update_skymodels(self, iter):
        """
        Updates the source and calibration sky models from the output sector sky model(s)
        """
        # Concat all output sector sky models
        sector_skymodels = [sector.get_output_skymodel_filename() for sector in self.sectors]
        skymodel = lsmtool.load(sector_skymodels[0])
        sector_skymodels.pop(0)
        for s2 in sector_skymodels:
            skymodel.concatenate(s2)
        self.make_skymodels(skymodel)

        # Re-adjust sector boundaries and update their sky models
        self.adjust_sector_boundaries()
        for sector in self.sectors:
            sector.make_skymodels()

    def define_sectors(self):
        """
        Defines the imaging sectors
        """
        nsectors_ra = self.parset['imaging_specific']['nsectors_ra']
        nsectors_dec = int(np.ceil(nsectors_ra / np.sin(self.mean_el_rad)))
        if nsectors_ra == 1 and nsectors_dec == 1:
            nsectors_ra == 0
        if nsectors_ra == 0:
            # Use a single sector
            nsectors_dec = 1
            width_ra = self.fwhm_ra_deg
            width_dec = self.fwhm_dec_deg
            center_x, center_y = self.radec2xy([self.ra], [self.dec])
            x = np.array([center_x])
            y = np.array([center_y])
        else:
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

        # Initialize the sectors
        self.sectors = []
        n = 0
        for i in range(nsectors_ra):
            for j in range(nsectors_dec):
                if nsectors_ra == 0:
                    name = 'field'
                else:
                    name = 'sector_{0}'.format(n)
                n += 1
                ra, dec = self.xy2radec([x[j, i]], [y[j, i]])
                self.sectors.append(Sector(name, ra[0], dec[0], width_ra, width_dec, self))

        # Adjust sector boundaries to avoid known sources and update their sky models
        if nsectors_ra > 0:
            self.adjust_sector_boundaries()
            self.log.info('Making sector sky models (for predicting)...')
            for sector in self.sectors:
                sector.make_skymodels()

        for sector in self.sectors:
            # Set the imaging parameters for selfcal
            sector.set_imaging_parameters(self.parset['imaging_specific']['selfcal_cellsize_arcsec'],
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
            sector.flag_abstime = self.flag_abstime
            sector.flag_baseline = self.flag_baseline
            sector.flag_freqrange = self.flag_freqrange
            sector.flag_expr = self.flag_expr

    def find_intersecting_sources(self):
        """
        Finds sources in the source sky model that intersect with the sector polygons
        """
        idx = rtree.index.Index()
        skymodel = self.source_skymodel
        RA = skymodel.getColValues('Ra')
        Dec = skymodel.getColValues('Dec')
        x, y = self.radec2xy(RA, Dec)
        sizes = skymodel.getPatchSizes(units='degree', weight=False)

        for i, (xs, ys, ss) in enumerate(zip(x, y, sizes)):
            xmin = xs - (ss / 2.0 / self.wcs_pixel_scale)
            xmax = xs + (ss / 2.0 / self.wcs_pixel_scale)
            ymin = ys - (ss / 2.0 / self.wcs_pixel_scale)
            ymax = ys + (ss / 2.0 / self.wcs_pixel_scale)
            idx.insert(i, (xmin, ymin, xmax, ymax))

        # For each sector side, query the index to find intersections
        intersecting_ind = []
        for sector in self.sectors:
            xmin, ymin, xmax, ymax = sector.initial_poly.bounds
            side1 = (xmin-1, ymin, xmin+1, ymax)
            intersecting_ind.extend(list(idx.intersection(side1)))
            side2 = (xmax-1, ymin, xmax+1, ymax)
            intersecting_ind.extend(list(idx.intersection(side2)))
            side3 = (xmin, ymin-1, xmax, ymin+1)
            intersecting_ind.extend(list(idx.intersection(side3)))
            side4 = (xmin, ymax-1, xmax, ymax+1)
            intersecting_ind.extend(list(idx.intersection(side4)))

        # Make point polys
        xfilt = np.array(x)[(np.array(intersecting_ind),)]
        yfilt = np.array(y)[(np.array(intersecting_ind),)]
        sfilt = np.array(sizes)[(np.array(intersecting_ind),)]
        points = [Point(xp, yp).buffer(sp/self.wcs_pixel_scale) for xp, yp, sp in
                  zip(xfilt, yfilt, sfilt)]
        return points

    def adjust_sector_boundaries(self):
        """
        Adjusts the imaging sector boundaries for overlaping sources
        """
        self.log.info('Adusting sector boudaries to avoid sources...')
        intersecting_source_polys = self.find_intersecting_sources()

        for sector in self.sectors:
            print(sector.name)
            for i in range(3):
                print(i)
                # Adjust boundaries for intersection with sources
                prev_poly = Polygon(sector.poly)
                for p2 in intersecting_source_polys:
                    if sector.poly.contains(p2.centroid):
                        # If point is inside, union the sector poly with the source one
                        self.log.info('inside')
                        sector.poly = sector.poly.union(p2)
                    else:
                        # If centroid of point is outside, difference the sector poly with
                        # the source one
                        self.log.info('outside')
                        sector.poly = sector.poly.difference(p2)
                if sector.poly.equals(prev_poly):
                    print('break')
                    break

            # Make sector region and vertices files
            sector.make_vertices_file()
            sector.make_region_file(os.path.join(self.working_dir, 'regions',
                                                 '{}_region_ds9.reg'.format(sector.name)))

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

        self.wcs_pixel_scale = 0.0027777778 # degrees/pixel (= 10"/pixel)
        w = WCS(naxis=2)
        w.wcs.crpix = [1000, 1000]
        w.wcs.cdelt = np.array([-self.wcs_pixel_scale, self.wcs_pixel_scale])
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
                output_image_file = input_image_file + '_blanked'
                blanked_images.append(output_image_file)
                blank_image.main(input_image_file, vertices_file, output_image_file)

            # Make the mosaic
            outfile = os.path.join(self.parset['dir_working'], 'results', 'image_{}'.format(iter),
                                   'field_mosaic.fits')
            self.output_image_filename = outfile
            mosaic_images.main(blanked_images, outfile)
        else:
            self.output_image_filename = self.sectors[0].get_output_image_filename(image_id)

        # Create sym links to image files
        dst_dir = os.path.join(self.parset['dir_working'], 'images', 'image_{}'.format(iter))
        create_directory(dst_dir)
        dst = os.path.join(dst_dir, 'field-MFS-image.fits')
        if os.path.exists(dst):
            os.unlink(dst)
        os.symlink(self.direction.get_output_image_filename(), dst)

