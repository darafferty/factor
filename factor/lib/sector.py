"""
Definition of the Sector class that holds parameters for each iamge sector
set
"""
import logging
import numpy as np
from astropy.coordinates import Angle
from shapely.geometry import Point, Polygon
from shapely.prepared import prep
from PIL import Image, ImageDraw
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

        # Define the initial sector polygon vertices and sky model
        self.define_vertices()

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
            sizes_arcmin = self.get_source_sizes_arcmin('inside')
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

    def make_skymodels(self):
        """
        Makes predict and source sky models from the parent field sky models
        """
        # Load and filter the predict sky model
        skymodel = self.field.calibration_skymodel.copy()
        skymodel = self.filter_skymodel(skymodel)

        # Write filtered sky model to file
        skymodel.setPatchPositions(method='wmean')
        self.predict_skymodel_file = os.path.join(self.field.working_dir, 'skymodels',
                                                      '{}_predict_skymodel.txt'.format(self.name))
        skymodel.write(self.predict_skymodel_file, clobber=True)

        # Save list of patches (directions) in the format written by DDECal in the h5parm
        self.patches = '[{}]'.format(','.join(['[{}]'.format(p) for p in skymodel.getPatchNames()]))

        # Find nearest patch to sector center
        patch_dist = skymodel.getDistance(self.ra, self.dec, byPatch=True).tolist()
        patch_names = skymodel.getPatchNames()
        self.central_patch = patch_names[patch_dist.index(min(patch_dist))]

        # Load and filter the source sky model
        skymodel = self.field.source_skymodel.copy()
        skymodel = self.filter_skymodel(skymodel)

        # Write filtered sky model to file
        skymodel.setPatchPositions(method='wmean')
        self.source_skymodel_file = os.path.join(self.field.working_dir, 'skymodels',
                                                      '{}_source_skymodel.txt'.format(self.name))
        skymodel.write(self.source_skymodel_file, clobber=True)

    def filter_skymodel(self, skymodel):
        """
        Filters input skymodel to select only sources that lie inside the sector

        Parameters
        ----------
        skymodel : LSMTool skymodel object
            Input sky model

        Returns
        -------
        filtered_skymodel : LSMTool skymodel object
            Filtered sky model
        """
        # Make list of sources
        RA = skymodel.getColValues('Ra')
        Dec = skymodel.getColValues('Dec')
        x, y = self.field.radec2xy(RA, Dec)
        x = np.array(x)
        y = np.array(y)
        xsize = int(1.1 * (max(x) - min(x)))
        ysize = int(1.1 * (max(y) - min(y)))
        inside = np.zeros(len(skymodel), dtype=bool)
        prepared_polygon = prep(self.poly)

        # Unmask everything outside of the polygon + its border (outline)
        mask = Image.new('L', (xsize, ysize), 0)
        verts = [(yv, xv) for xv, yv in zip(self.poly.exterior.coords.xy[0],
                                            self.poly.exterior.coords.xy[1])]
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=1)
        inside_ind = np.where(np.array(mask)[(x.astype(int), y.astype(int))])
        inside[inside_ind] = True

        # Now check sources in the border precisely
        mask = Image.new('L', (xsize, ysize), 0)
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=0)
        border_ind = np.where(np.array(mask)[(x.astype(int), y.astype(int))])
        points = [Point(xs, ys) for xs, ys in zip(x[border_ind], y[border_ind])]
        for i, p in enumerate(points):
            p.index = border_ind[0][i]
        outside_points = filter(lambda v: not prepared_polygon.contains(v), points)
        0/0
        for outside_point in outside_points:
            inside[outside_point.index] = False
        skymodel.select(inside)
        return skymodel

    def get_source_sizes_arcmin(self):
        """
        Returns list of source sizes in arcmin

        Returns
        -------
        sizes : list
            List of source sizes in arcmin
        """
        lsmtool._logging.setLevel('debug')
        skymodel = lsmtool.load(self.source_skymodel_file)
        sizes = skymodel.getPatchSizes(units='arcmin', weight=False)
        return sizes

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
        self.initial_poly = poly
        self.poly = Polygon(poly)

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

