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
import glob


class Field(object):
    """
    The Field object stores parameters needed for processing of the field

    Parameters
    ----------
    parset : dict
        Parset with processing parameters
    minimal : bool
        If True, only initialize the minimal required parameters
    """
    def __init__(self, parset, mininmal=False):
        # Initialize basic attributes
        self.name = 'field'
        self.log = logging.getLogger('factor:{}'.format(self.name))
        self.parset = parset.copy()
        self.working_dir = self.parset['dir_working']
        self.ms_filenames = self.parset['mss']
        self.numMS = len(self.ms_filenames)
        self.data_colname = 'DATA'
        self.flag_abstime = self.parset['flag_abstime']
        self.flag_baseline = self.parset['flag_baseline']
        self.flag_freqrange = self.parset['flag_freqrange']
        self.flag_expr = self.parset['flag_expr']
        self.input_h5parm = self.parset['input_h5parm']
        self.solve_min_uv_lambda = self.parset['calibration_specific']['solve_min_uv_lambda']
        self.solve_tecandphase = self.parset['calibration_specific']['solve_tecandphase']
        self.approximatetec = self.parset['calibration_specific']['approximatetec']
        self.propagatesolutions = self.parset['calibration_specific']['propagatesolutions']
        self.maxapproxiter = self.parset['calibration_specific']['maxapproxiter']
        self.maxiter = self.parset['calibration_specific']['maxiter']
        self.stepsize = self.parset['calibration_specific']['stepsize']
        self.tolerance = self.parset['calibration_specific']['tolerance']

        if not mininmal:
            # Scan MS files to get observation info
            self.scan_observations()

            # Make calibration and source sky models by grouping the initial sky model
            self.make_skymodels(self.parset['input_skymodel'], self.parset['regroup_input_skymodel'])

            # Set up imaging sectors
            self.makeWCS()
            self.define_sectors()

    def scan_observations(self):
        """
        Checks input MS files and initializes the associated Observation objects
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

    def get_obs_parameters(self, parameter):
        """
        Returns list of parameters for all observations

        Parameters
        ----------
        parameter : str
            Name of parameter to return

        Returns
        -------
        parameters : list
            List of parameters of each observation
        """
        return sum([obs.parameters[parameter] for obs in self.observations], [])

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

        # Check if any sky models included in Factor are within the region of the
        # input sky model. If so, concatenate them with the input sky model
        if self.parset['calibration_specific']['use_included_skymodels']:
            max_separation_deg = self.fwhm_deg * 2.0
            factor_lib_dir = os.path.dirname(os.path.abspath(__file__))
            skymodel_dir = os.path.join(os.path.split(factor_lib_dir)[0], 'skymodels')
            skymodels = glob.glob(os.path.join(skymodel_dir, '*.skymodel'))
            concat_skymodels = []
            for skymodel in skymodels:
                try:
                    s = lsmtool.load(skymodel)
                    dist_deg = s.getDistance(self.ra, self.dec)
                    if any(dist_deg < max_separation_deg):
                        concat_skymodels.append(skymodel)
                except IOError:
                    pass
            matching_radius_deg = 30.0 / 3600.0  # => 30 arcsec
            for s in concat_skymodels:
                skymodel.concatenate(s, matchBy='position', radius=matching_radius_deg,
                                     keep='from2', inheritPatches=True)

        # Group by thresholding
        self.log.info('Identifying sources...')
        source_skymodel = skymodel.copy()
        source_skymodel.group('threshold', FWHM='60.0 arcsec', threshold=0.05)
        self.source_skymodel = source_skymodel

        # Now tesselate to get patches of the target flux and write out calibration sky model
        if regroup:
            flux = self.parset['calibration_specific']['patch_target_flux_jy']
            self.log.info('Grouping sky model to form calibration patches of ~ {} Jy each...'.format(flux))
            source_skymodel.group(algorithm='tessellate', targetFlux=flux, method='mid', byPatch=True)
            calibration_skymodel = source_skymodel
        else:
            calibration_skymodel = skymodel
        self.log.info('Using {} calibration patches'.format(len(calibration_skymodel.getPatchNames())))
        self.calibration_skymodel_file = os.path.join(self.working_dir, 'skymodels', 'calibration_skymodel.txt')
        calibration_skymodel.write(self.calibration_skymodel_file, clobber=True)
        self.calibration_skymodel = calibration_skymodel

    def update_skymodels(self, iter):
        """
        Updates the source and calibration sky models from the output sector sky model(s)

        Parameters
        ----------
        iter : int
            Iteration index
        """
        # Concat all output sector sky models
        self.log.info('Updating sky model...')
        if self.parset['strategy'] == 'sectorselfcal':
            # Use new models from the imaged sectors only
            sector_skymodels = [sector.image_skymodel_file for sector in self.imaging_sectors]
            sector_names = [sector.name for sector in self.imaging_sectors]

            # Each sector is a single patch, so don't regroup
            regroup = False
        else:
            # Use models from all sectors, whether imaged or not
            sector_skymodels = []
            for sector in self.sectors:
                if sector.is_outlier:
                    sector_skymodels.append(sector.predict_skymodel_file)
                else:
                    sector_skymodels.append(sector.image_skymodel_file)
            sector_names = [sector.name for sector in self.sectors]

            # Regroup to make patches with the target flux density
            regroup = True
        for i, (sm, sn) in enumerate(zip(sector_skymodels, sector_names)):
            if i == 0:
                skymodel = lsmtool.load(sm)
                skymodel.group('single', root=sn)
            else:
                skymodel2 = lsmtool.load(sm)
                skymodel2.group('single', root=sn)
                skymodel.concatenate(skymodel2)
        self.make_skymodels(skymodel, regroup=regroup)

        # Re-adjust sector boundaries and update their sky models
        if not self.parset['strategy'] == 'sectorselfcal':
            self.adjust_sector_boundaries()
        for sector in self.sectors:
            sector.make_skymodel()

    def define_sectors(self):
        """
        Defines the imaging sectors
        """
        self.imaging_sectors = []

        # Determine whether we use a user-supplied list of sectors or a grid
        if len(self.parset['imaging_specific']['sector_center_ra_list']) > 0:
            # Use user-supplied list
            sector_center_ra_list = self.parset['imaging_specific']['sector_center_ra_list']
            sector_center_dec_list = self.parset['imaging_specific']['sector_center_dec_list']
            sector_width_ra_deg_list = self.parset['imaging_specific']['sector_width_ra_deg_list']
            sector_width_dec_deg_list = self.parset['imaging_specific']['sector_width_dec_deg_list']
            sector_do_multiscale_list = self.parset['imaging_specific']['sector_do_multiscale_list']
            n = 1
            for ra, dec, width_ra, width_dec in zip(sector_center_ra_list, sector_center_dec_list,
                                                    sector_width_ra_deg_list, sector_width_dec_deg_list):
                name = 'sector_{0}'.format(n)
                self.imaging_sectors.append(Sector(name, ra, dec, width_ra, width_dec, self))
                n += 1
            self.log.info('Using {0} user-defined imaging sector(s)'.format(len(self.imaging_sectors)))
            # TODO: check whether flux density in each sector meets minimum and warn if not?
        else:
            # Make a regular grid of sectors
            if self.parset['imaging_specific']['grid_center_ra'] is None:
                image_ra = self.ra
            else:
                image_ra = self.parset['imaging_specific']['grid_center_ra']
            if self.parset['imaging_specific']['grid_center_dec'] is None:
                image_dec = self.dec
            else:
                image_dec = self.parset['imaging_specific']['grid_center_dec']
            if self.parset['imaging_specific']['grid_width_ra_deg'] is None:
                image_width_ra = self.fwhm_ra_deg
            else:
                image_width_ra = self.parset['imaging_specific']['grid_width_ra_deg']
            if self.parset['imaging_specific']['grid_width_dec_deg'] is None:
                image_width_dec = self.fwhm_dec_deg
            else:
                image_width_dec = self.parset['imaging_specific']['grid_width_dec_deg']

            nsectors_ra = self.parset['imaging_specific']['grid_nsectors_ra']
            if nsectors_ra == 0:
                # Force a single sector
                nsectors_ra = 1
                nsectors_dec = 1
            else:
                nsectors_dec = int(np.ceil(nsectors_ra / np.sin(self.mean_el_rad)))
            if nsectors_ra == 1 and nsectors_dec == 1:
                # Make a single sector
                nsectors_dec = 1
                width_ra = image_width_ra
                width_dec = image_width_dec
                center_x, center_y = self.radec2xy([image_ra], [image_dec])
                x = np.array([center_x])
                y = np.array([center_y])
            else:
                # Make the grid
                width_ra = image_width_ra / nsectors_ra
                width_dec = image_width_dec / nsectors_dec
                width_x = width_ra / abs(self.wcs.wcs.cdelt[0])
                width_y = width_dec / abs(self.wcs.wcs.cdelt[1])
                center_x, center_y = self.radec2xy([image_ra], [image_dec])
                min_x = center_x - width_x / 2.0 * (nsectors_ra - 1)
                max_x = center_x + width_x / 2.0 * (nsectors_ra - 1)
                min_y = center_y - width_y / 2.0 * (nsectors_dec - 1)
                max_y = center_y + width_y / 2.0 * (nsectors_dec - 1)
                x = np.linspace(min_x, max_x, nsectors_ra)
                y = np.linspace(min_y, max_y, nsectors_dec)
                x, y = np.meshgrid(x, y)
                self.log.info('Using {0} imaging sectors ({1} in RA, {2} in Dec)'.format(
                              nsectors_ra*nsectors_dec, nsectors_ra, nsectors_dec))

            # Initialize the sectors in the grid
            n = 1
            for i in range(nsectors_ra):
                for j in range(nsectors_dec):
                    name = 'sector_{0}'.format(n)
                    ra, dec = self.xy2radec([x[j, i]], [y[j, i]])
                    self.imaging_sectors.append(Sector(name, ra[0], dec[0], width_ra, width_dec, self))
                    n += 1

        # Adjust sector boundaries to avoid known sources and update their sky models
        self.adjust_sector_boundaries()
        self.log.info('Making sector sky models (for predicting)...')
        for sector in self.imaging_sectors:
            sector.make_skymodel()

        # Set the imaging parameters for each imaging sector
        sector_do_multiscale_list = self.parset['imaging_specific']['sector_do_multiscale_list']
        for i, sector in enumerate(self.imaging_sectors):
            if len(sector_do_multiscale_list) > 0:
                do_multiscale = sector_do_multiscale_list[i]
            else:
                do_multiscale = None
            sector.set_imaging_parameters(do_multiscale)

            # Transfer any field flagging/calibration parameters so they are also used
            # during imaging
            sector.flag_abstime = self.flag_abstime
            sector.flag_baseline = self.flag_baseline
            sector.flag_freqrange = self.flag_freqrange
            sector.flag_expr = self.flag_expr
            sector.solve_tecandphase = self.solve_tecandphase

        # Make outlier sectors containing any remaining calibration sources (not
        # included in any sector sky model). These sectors are not imaged; they are only
        # used in prediction and subtraction
        self.outlier_sectors = []
        outlier_skymodel = self.make_outlier_skymodel()
        nsources = len(outlier_skymodel)
        if nsources > 0:
            nnodes = 10  # TODO: tune to number of available nodes and/or memory?
            for i in range(nnodes):
                outlier_sector = Sector('outlier_{0}'.format(i), self.ra, self.dec, 1.0, 1.0, self)
                outlier_sector.is_outlier = True
                outlier_sector.predict_skymodel = outlier_skymodel.copy()
                startind = i * int(nsources/nnodes)
                if i == nnodes-1:
                    endind = nsources
                else:
                    endind = startind + int(nsources/nnodes)
                outlier_sector.predict_skymodel.select(np.array(range(startind, endind)))
                outlier_sector.make_skymodel()
                self.outlier_sectors.append(outlier_sector)

        # Finally, make a list of all sectors
        self.sectors = self.imaging_sectors[:] + self.outlier_sectors

    def find_intersecting_sources(self):
        """
        Finds sources that intersect with the intial sector boundaries

        Returns
        -------
        intersecting_source_polys: list of Polygons
            List of source polygons that intersect one or more sector boundaries
        """
        idx = rtree.index.Index()
        skymodel = self.source_skymodel
        RA, Dec = skymodel.getPatchPositions(asArray=True)
        x, y = self.radec2xy(RA, Dec)
        sizes = skymodel.getPatchSizes(units='degree')
        minsize = 1  #  minimum allowed source size in pixels
        sizes = [max(minsize, s/2.0/self.wcs_pixel_scale) for s in sizes]  # radii in pixels

        for i, (xs, ys, ss) in enumerate(zip(x, y, sizes)):
            xmin = xs - ss
            xmax = xs + ss
            ymin = ys - ss
            ymax = ys + ss
            idx.insert(i, (xmin, ymin, xmax, ymax))

        # For each sector side, query the index to find any intersections
        intersecting_ind = []
        buffer = 2  # how many pixels away from each side to check
        for sector in self.imaging_sectors:
            xmin, ymin, xmax, ymax = sector.initial_poly.bounds
            side1 = (xmin-buffer, ymin, xmin+buffer, ymax)
            intersecting_ind.extend(list(idx.intersection(side1)))
            side2 = (xmax-buffer, ymin, xmax+buffer, ymax)
            intersecting_ind.extend(list(idx.intersection(side2)))
            side3 = (xmin, ymin-buffer, xmax, ymin+buffer)
            intersecting_ind.extend(list(idx.intersection(side3)))
            side4 = (xmin, ymax-buffer, xmax, ymax+buffer)
            intersecting_ind.extend(list(idx.intersection(side4)))

        # Make polygons for intersecting sources, with a size = 1.5 * radius of source
        if len(intersecting_ind) > 0:
            xfilt = np.array(x)[(np.array(intersecting_ind),)]
            yfilt = np.array(y)[(np.array(intersecting_ind),)]
            sfilt = np.array(sizes)[(np.array(intersecting_ind),)]
            intersecting_source_polys = [Point(xp, yp).buffer(sp*1.5) for
                                         xp, yp, sp in zip(xfilt, yfilt, sfilt)]
        else:
            intersecting_source_polys = []
        return intersecting_source_polys

    def adjust_sector_boundaries(self):
        """
        Adjusts the imaging sector boundaries for overlaping sources
        """
        self.log.info('Adusting sector boudaries to avoid sources...')
        intersecting_source_polys = self.find_intersecting_sources()
        for sector in self.imaging_sectors:
            # Make sure all sectors start from their initial polygons
            sector.poly = sector.initial_poly

        for sector in self.imaging_sectors:
            # Adjust boundaries for intersection with sources
            for p2 in intersecting_source_polys:
                if sector.poly.contains(p2.centroid):
                    # If point is inside, union the sector poly with the source one
                    sector.poly = sector.poly.union(p2)
                else:
                    # If centroid of point is outside, difference the sector poly with
                    # the source one
                    sector.poly = sector.poly.difference(p2)

            # Check whether the sector has been broken into two or more separate polygons
            if type(sector.poly) is not Polygon:
                # If so, keep largest one only
                sector.poly = sector.poly[np.argmax([p.area for p in sector.poly])]

            # Make sector region and vertices files
            sector.make_vertices_file()
            sector.make_region_file(os.path.join(self.working_dir, 'regions',
                                                 '{}_region_ds9.reg'.format(sector.name)))

    def make_outlier_skymodel(self):
        """
        Make a sky model of any outlier calibration sources, not included in any
        imaging sector
        """
        all_source_names = self.calibration_skymodel.getColValues('Name').tolist()
        sector_source_names = []
        for sector in self.imaging_sectors:
            skymodel = lsmtool.load(sector.predict_skymodel_file)
            sector_source_names.extend(skymodel.getColValues('Name').tolist())
        outlier_ind = np.array([all_source_names.index(sn) for sn in all_source_names
                                if sn not in sector_source_names])
        outlier_skymodel = self.calibration_skymodel.copy()
        outlier_skymodel.select(outlier_ind, force=True)
        return outlier_skymodel

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

        self.wcs_pixel_scale = 10.0 / 3600.0  # degrees/pixel (= 10"/pixel)
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

    def make_mosaic(self, iter):
        """
        Make mosaic of the sector images

        Parameters
        ----------
        iter : int
            Iteration index
        """
        dst_dir = os.path.join(self.parset['dir_working'], 'images', 'image_{}'.format(iter))
        create_directory(dst_dir)
        field_image_filename = os.path.join(dst_dir, 'field-MFS-image.fits')
        if os.path.exists(field_image_filename):
            os.remove(field_image_filename)
        if len(self.imaging_sectors) > 1:
            # Blank the sector images before making mosaic
            self.log.info('Making mosiac of sector images...')
            blanked_images = []
            for sector in self.imaging_sectors:
                input_image_file = sector.image_file
                vertices_file = sector.vertices_file
                output_image_file = input_image_file + '_blanked'
                blanked_images.append(output_image_file)
                blank_image.main(input_image_file, output_image_file, vertices_file=vertices_file)

            # Make the mosaic
            mosaic_images.main(blanked_images, field_image_filename)
        else:
            # No need to mosaic a single image; just copy it
            output_image_filename = self.imaging_sectors[0].image_file
            if os.path.exists(output_image_filename):
                os.remove(output_image_filename)
            os.system('cp {0} {1}'.format(output_image_filename, field_image_filename))
