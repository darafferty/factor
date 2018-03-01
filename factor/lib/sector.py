"""
Definition of the Sector class that holds parameters for each iamge sector
set
"""
import logging
import numpy as np
from astropy.coordinates import Angle
from shapely.geometry import Point, Polygon, box
from shapely.prepared import prep
import pickle
import lsmtool
import os
import copy


class Sector(object):
    """
    The Sector object contains various parameters for a sector of the field

    Parameters
    ----------
    name : str
        Name of sector
    ra : float
        RA in degrees of sector center
    dec : float
        Dec in degrees of sector center
    width_ra : float
        Width of sector in RA degrees
    width_dec : float
        Width of sector in Dec in degrees
    field : Field object
        Field object
    """
    def __init__(self, name, ra, dec, width_ra, width_dec, field):
        self.name = name
        self.log = logging.getLogger('factor:{0}'.format(self.name))
        if type(ra) is str:
            ra = Angle(ra).to('deg').value
        if type(dec) is str:
            dec = Angle(dec).to('deg').value
        self.ra = ra
        self.dec = dec
        self.width_ra = width_ra
        self.width_dec = width_dec
        self.field = field
        self.vertices_file = os.path.join(field.working_dir, 'regions', '{}_vertices.pkl'.format(self.name))
        self.region_file = '[]'
        self.output_image_filename = {}
        self.output_skymodel_filename = {}

        # Make copies of the observation objects, as each sector may have its own
        # observation-specific settings
        self.observations = []
        for obs in field.observations:
            obs.log = None
            cobs = copy.deepcopy(obs)
            cobs.log = logging.getLogger('factor:{}'.format(cobs.name))
            self.observations.append(cobs)

        # Define the sector polygon vertices and sky model
        self.define_vertices()
        self.make_skymodel()
        self.adjust_boundary()

    def set_imaging_parameters(self, cellsize_arcsec, robust, taper_arcsec,
                               min_uv_lambda, max_uv_lambda, max_peak_smearing,
                               wsclean_bl_averaging=False, multiscale_scales_pixel=None,
                               use_idg=True, idg_mode='hybrid'):
        """
        Sets the imaging parameters for given values

        Parameters
        ----------
        cellsize_arcsec : float
            Pixel size in arcsec for imaging
        robust : float
            Briggs robust parameter for imaging
        taper_arcsec : float
            Taper in arcsec for imaging
        min_uv_lambda : float
            Minimum uv cut in lamdba
        max_uv_lambda : float
            Maximum uv cut in lamdba
        max_peak_smearing : float
            Maximum allowed peak flux density reduction
        wsclean_bl_averaging : bool, optional
            Use baseline-dependent averaging in WSClean
        multiscale_scales_pixel : list of float, optional
            List of scales for WSClean when multi-scale is used
        use_idg : bool, optional
            Use image domain gridder in WSClean
        idg_mode : str, optional
            IDG mode to use
        """
        self.cellsize_arcsec = cellsize_arcsec
        self.cellsize_deg = cellsize_arcsec / 3600.0
        self.robust = robust
        self.taper_arcsec = taper_arcsec
        self.min_uv_lambda = min_uv_lambda
        self.max_uv_lambda = max_uv_lambda
        self.use_idg = use_idg
        self.idg_mode = idg_mode

        # Set ID for current imaging parameters
        self.current_image_id = self.get_image_id(cellsize_arcsec, robust, taper_arcsec,
                                                  min_uv_lambda, max_uv_lambda)

        # Set image size
        self.imsize = [int(self.width_ra / self.cellsize_deg * 1.1),
                       int(self.width_dec / self.cellsize_deg * 1.1)]
        if self.use_idg:
            # IDG does not yet support rectangular images
            self.imsize = [max(self.imsize), max(self.imsize)]
        self.wsclean_imsize = '{0} {1}'.format(self.imsize[0], self.imsize[1])
        self.log.debug('Image size is {0} x {1} pixels'.format(
                       self.imsize[0], self.imsize[1]))

        # Set number of output channels to get 5 channels, but of no less than ~ 2 MHz each
        target_nchannels = 4
        tot_bandwidth = 0.0
        for obs in self.observations:
            # Find observation with largest bandwidth
            obs_bandwidth = obs.numchannels * obs.channelwidth
            if obs_bandwidth > tot_bandwidth:
                tot_bandwidth = obs_bandwidth
        self.wsclean_nchannels = min(target_nchannels, int(np.ceil(tot_bandwidth / 2e6)))

        # Set number of iterations and threshold
        scaling_factor = np.sqrt(np.float(tot_bandwidth / 2e6))
        self.wsclean_niter = int(12000 * scaling_factor)

        # Set multiscale: get source sizes and check for large sources
        self.multiscale = None
        large_size_arcmin = 4.0
        if self.multiscale is None:
            sizes_arcmin, _, _ = self.get_source_sizes_arcmin('inside')
            if sizes_arcmin is not None and any([s > large_size_arcmin for s in sizes_arcmin]):
                self.multiscale = True
            else:
                self.multiscale = False
        if self.multiscale:
            self.multiscale_scales_pixel = multiscale_scales_pixel
            self.wsclean_niter /= 2 # fewer iterations are needed
            self.log.debug("Will do multiscale cleaning.")
        else:
            self.multiscale_scales_pixel = 0

        # Set the observation-specific parameters
        for obs in self.observations:
            # Set filename for model-subtracted data that matches the one made by the
            # calibrate pipeline
            ms_subtracted_filename = '{0}.sector_{1}_sub'.format(obs.ms_filename,
                                                            self.name.split('_')[1])
            # Set imaging parameters
            obs.set_imaging_parameters(cellsize_arcsec, max_peak_smearing,
                                       self.width_ra, self.width_dec, ms_subtracted_filename)

        # Set BL-dependent averaging
        if wsclean_bl_averaging:
            timestep_sec = (self.observations[0].timepersample *
                            self.observations[0].imaging_parameters['image_timestep'])
            self.wsclean_nwavelengths = self.get_nwavelengths(self.cellsize_deg,
                                                              timestep_sec)
        else:
            self.wsclean_nwavelengths = 0

    def get_image_id(self, cellsize_arcsec, robust, taper_arcsec, min_uv_lambda, max_uv_lambda):
        """
        Returns the imaging ID for given values

        Parameters
        ----------
        cellsize_arcsec : float
            Pixel size in arcsec for imaging
        robust : float
            Briggs robust parameter for imaging
        taper_arcsec : float
            Taper in arcsec for imaging
        min_uv_lambda : float
            Minimum uv cut in lamdba
        max_uv_lambda : float
            Maximum uv cut in lamdba
        """
        return '{0}_{1}_{2}_{3}_{4}'.format(cellsize_arcsec, robust, taper_arcsec,
                                            min_uv_lambda, max_uv_lambda)

    def store_output_image_filename(self, filename):
        """
        Stores path to output image filename

        Parameters
        ----------
        filename : str
            Path to file
        """
        self.output_image_filename[self.current_image_id] = filename

    def store_output_skymodel_filename(self, filename):
        """
        Stores path to output skymodel filename

        Parameters
        ----------
        filename : str
            Path to file
        """
        self.output_skymodel_filename[self.current_image_id] = filename

    def get_output_image_filename(self, image_id=None):
        """
        Returns path to output image for given ID

        Parameters
        ----------
        image_id : str, optional
            Imaging ID
        """
        if image_id is None:
            image_id = self.current_image_id
        try:
            return self.output_image_filename[image_id]
        except KeyError:
            return None

    def get_output_skymodel_filename(self, image_id=None):
        """
        Returns path to output sky model for given ID

        Parameters
        ----------
        image_id : str, optional
            Imaging ID
        """
        if image_id is None:
            image_id = self.current_image_id
        try:
            return self.output_skymodel_filename[image_id]
        except KeyError:
            return None

    def get_nwavelengths(self, cellsize_deg, timestep_sec):
        """
        Returns nwavelengths for WSClean BL-based averaging

        The value depends on the integration time given the specified maximum
        allowed smearing. We scale it from the imaging cell size assuming normal
        sampling as:

        max baseline in nwavelengths = 1 / theta_rad ~= 1 / (cellsize_deg * 3 * pi / 180)
        nwavelengths = max baseline in nwavelengths * 2 * pi * integration time in seconds / (24 * 60 * 60) / 4

        Parameters
        ----------
        cellsize_deg : float
            Pixel size of image in degrees
        timestep_sec : float
            Length of one timestep in seconds

        """
        max_baseline = 1 / (3 * cellsize_deg * np.pi / 180)
        wsclean_nwavelengths_time = int(max_baseline * 2*np.pi * timestep_sec /
            (24 * 60 * 60) / 4)
        return wsclean_nwavelengths_time

    def make_skymodel(self):
        """
        Makes a sky model for the sector from the parent field calibration sky model
        """
        skymodel = lsmtool.load(self.field.calibration_skymodel_file)

        # Make list of sources in full sky model
        RA = skymodel.getColValues('Ra')
        Dec = skymodel.getColValues('Dec')
        x, y = self.field.radec2xy(RA, Dec)
        x = np.array(x)
        y = np.array(y)

        # Identify sources near poly bounding box
        minx, miny, maxx, maxy = self.poly.bounds
        delx = (maxx - minx) * 0.2
        minx -= delx
        maxx += delx
        dely = (maxy - miny) * 0.2
        miny -= dely
        maxy += dely
        inner_poly = self.get_inner_poly()
        if inner_poly is not None:
            iminx, iminy, imaxx, imaxy = inner_poly.bounds
            inside_ind = np.where( (x > iminx) & (x < imaxx) &
                                   (y > iminy) & (y < imaxy) )
            near_ind = np.where( ((x > minx) & (x < iminx) & (y > miny) & (y < maxy)) |
                                 ((x > imaxx) & (x < maxx) & (y > miny) & (y < maxy)) |
                                 ((y > miny) & (y < iminy) & (x > minx) & (x < maxx)) |
                                 ((y > imaxy) & (y < maxy) & (x > minx) & (x < maxx)) )
        else:
            inside_ind = [[]]
            near_ind = np.where((x > minx) & (y > miny) &
                                (x < maxx) & (y < maxy))

        points = []
        near_boundary = np.zeros(len(skymodel), dtype=bool)
        inside = np.zeros(len(skymodel), dtype=bool)
        for i, (xp, yp) in enumerate(zip(x, y)):
            if i in near_ind[0]:
                p = Point((xp, yp))
                p.index = i
                points.append(p)
                near_boundary[i] = True
            if i in inside_ind[0]:
                inside[i] = True

        # Make sky model with sources near boundary only (for source avoidance)
        boundary_skymodel = skymodel.copy()
        boundary_skymodel.select(near_boundary)
        self.num_sources_near_boundary = len(boundary_skymodel)
        if self.num_sources_near_boundary > 0:
            self.boundary_skymodel_file = os.path.join(self.field.working_dir, 'skymodels',
                                                   '{}_boundary_skymodel.txt'.format(self.name))
            boundary_skymodel.write(self.boundary_skymodel_file, clobber=True)

        # Find all sources that are inside the sector
        prepared_polygon = prep(self.poly)
        intersecting_points = filter(prepared_polygon.contains, points)
        for p in intersecting_points:
            inside[p.index] = True
        skymodel.select(inside)

        # Write sky model to file
        skymodel.setPatchPositions(method='wmean')
        self.calibration_skymodel_file = os.path.join(self.field.working_dir, 'skymodels',
                                                      '{}_skymodel.txt'.format(self.name))
        skymodel.write(self.calibration_skymodel_file, clobber=True)

        # Save list of patches (directions) in the format written by DDECal in the h5parm
        self.patches = '[{}]'.format(','.join(['[{}]'.format(p) for p in skymodel.getPatchNames()]))

        # Find nearest patch to sector center
        patch_dist = skymodel.getDistance(self.ra, self.dec, byPatch=True).tolist()
        patch_names = skymodel.getPatchNames()
        self.central_patch = patch_names[patch_dist.index(min(patch_dist))]

    def get_source_sizes_arcmin(self, filter='inside'):
        """
        Returns list of source sizes in arcmin

        Parameters
        ----------
        filter : str, optional
            One of "inside" or "near boundary"

        Returns
        -------
        sizes : list
            List of source sizes in arcmin
        """
        lsmtool._logging.setLevel('debug')
        if filter == 'inside':
            skymodel = lsmtool.load(self.calibration_skymodel_file)
        elif filter == 'near boundary':
            skymodel = lsmtool.load(self.boundary_skymodel_file)

        skymodel.group('threshold', FWHM='60.0 arcsec')
        sizes = skymodel.getPatchSizes(units='arcmin', weight=False)
        RA, Dec = skymodel.getPatchPositions(asArray=True)
        return sizes, RA, Dec

    def get_obs_parameters(self, parameter):
        """
        Returns list of imaging parameters for all observations

        Parameters
        ----------
        parameter : str
            Name of imaging parameter to return

        Returns
        -------
        parameters : list
            List of parameters, with one entry for each observation
        """
        return [obs.imaging_parameters[parameter] for obs in self.observations]

    def define_vertices(self):
        """
        Determines the vertices of the sector polygon
        """
        # Define initial polygon as a rectangle
        sx, sy = self.field.radec2xy([self.ra], [self.dec])
        ra_width_pix = self.width_ra / abs(self.field.wcs.wcs.cdelt[0])
        dec_width_pix = self.width_dec / abs(self.field.wcs.wcs.cdelt[1])
        x0 = sx[0] - ra_width_pix / 2.0
        y0 = sy[0] - dec_width_pix / 2.0
        poly_verts = [(x0, y0), (x0, y0+dec_width_pix),
                      (x0+ra_width_pix, y0+dec_width_pix),
                      (x0+ra_width_pix, y0), (x0, y0)]
        poly = Polygon(poly_verts)
        self.poly = poly

    def adjust_boundary(self):
        """
        Adjusts the boundary of the sector for known sources
        """
        if self.num_sources_near_boundary == 0:
            return

        # Find nearby sources in input sky model and adjust sector boundaries
        # if necessary
        poly = self.poly
        sizes, RA, Dec = self.get_source_sizes_arcmin('near boundary')

        # Make buffered points for all sources
        points = []
        sx, sy = self.field.radec2xy(RA, Dec)
        for xp, yp, sp in zip(sx, sy, sizes):
            radius = sp * 1.2 / 2.0 / self.field.wcs.wcs.cdelt[0] # size of source in pixels
            points.append(Point((xp, yp)).buffer(radius))

        # Alter sector polygon to avoid sources in the input sky model
        niter = 0
        while niter < 3:
            niter += 1
            prepared_polygon = prep(poly)
            intersecting_points = filter(prepared_polygon.intersects, points)

            # Adjust sector polygon for each source that intersects it
            for p2 in intersecting_points:
                if poly.contains(p2.centroid):
                    # If point is inside, union the polys
                    poly = poly.union(p2)
                else:
                    # If centroid of point is outside, difference the polys
                    poly = poly.difference(p2)
        self.poly = poly

    def get_inner_poly(self):
        """
        Return the poly that fits inside the sector without touching it
        """
        # Find poly that fits inside
        minx, miny, maxx, maxy = self.poly.bounds
        stepx = (maxx - minx) / 10
        stepy = (maxy - miny) / 10
        success = False
        for i in range(10):
            # Shrink inner_poly until it's fully within poly
            minx += stepx
            maxx -= stepx
            miny += stepy
            maxy -= stepy
            inner_poly = box(minx, miny, maxx, maxy)
            if self.poly.contains(inner_poly):
                success = True
                break
        if success:
            return inner_poly
        else:
            return None

    def get_vertices_radec(self):
        """
        Return the vertices as RA, Dec for the sector boundary
        """
        ra, dec = self.field.xy2radec(self.poly.exterior.coords.xy[0].tolist(),
                                 self.poly.exterior.coords.xy[1].tolist())
        vertices = [np.array(ra), np.array(dec)]

        return vertices

    def make_vertices_file(self):
        """
        Make a vertices file for the sector boundary
        """
        vertices = self.get_vertices_radec()

        with open(self.vertices_file, 'wb') as f:
            pickle.dump(vertices, f)

    def make_region_file(self, outputfile, region_format='ds9'):
        """
        Make a ds9 or CASA region file for the sector boundary

        Parameters
        ----------
        outputfile : str
            Name of output region file
        region_format : str, optional
            Format of region file: 'ds9' or 'casa'
        """
        vertices = self.get_vertices_radec()

        if region_format == 'casa':
            lines = ['#CRTFv0\n\n']
            xylist = []
            RAs = vertices[0][0:-1] # trim last point, as it is a repeat of the first
            Decs = vertices[1][0:-1]
            for x, y in zip(RAs, Decs):
                xylist.append('[{0}deg, {1}deg]'.format(x, y))
            lines.append('poly[{0}]\n'.format(', '.join(xylist)))

            with open(outputfile, 'wb') as f:
                f.writelines(lines)
        elif region_format == 'ds9':
            lines = []
            lines.append('# Region file format: DS9 version 4.0\nglobal color=green '
                         'font="helvetica 10 normal" select=1 highlite=1 edit=1 '
                         'move=1 delete=1 include=1 fixed=0 source=1\nfk5\n')
            xylist = []
            RAs = vertices[0]
            Decs = vertices[1]
            for x, y in zip(RAs, Decs):
                xylist.append('{0}, {1}'.format(x, y))
            lines.append('polygon({0})\n'.format(', '.join(xylist)))
            lines.append('point({0}, {1}) # point=cross width=2 text={{{2}}}\n'.
                format(self.ra, self.dec, self.name))

            with open(outputfile, 'wb') as f:
                f.writelines(lines)
        else:
            self.log.error('Region format not understood.')

