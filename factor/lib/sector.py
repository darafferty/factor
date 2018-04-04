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
        self.is_outlier = False

        # Make copies of the observation objects, as each sector may have its own
        # observation-specific settings
        self.observations = []
        for obs in field.observations:
            obs.log = None # deepcopy cannot copy the log object
            cobs = copy.deepcopy(obs)
            obs.log = logging.getLogger('factor:{}'.format(obs.name))
            cobs.log = logging.getLogger('factor:{}'.format(cobs.name))
            self.observations.append(cobs)

        # Define the initial sector polygon vertices and sky model
        self.intialize_vertices()

    def set_predict_parameters(self):
        """
        Sets the predict parameters
        """
        for obs in self.observations:
            obs.set_predict_parameters(self.name, self.patches)

    def set_imaging_parameters(self, do_multiscale=None):
        """
        Sets the imaging parameters

        Parameters
        ----------
        do_multiscale : str
            If True, multiscale clean is done. If None, multiscale clean is done when a
            large source is detected
        """
        self.cellsize_arcsec = self.field.parset['imaging_specific']['cellsize_arcsec']
        self.cellsize_deg = self.cellsize_arcsec / 3600.0
        self.robust = self.field.parset['imaging_specific']['robust']
        self.taper_arcsec = self.field.parset['imaging_specific']['taper_arcsec']
        self.min_uv_lambda = self.field.parset['imaging_specific']['min_uv_lambda']
        self.max_uv_lambda = self.field.parset['imaging_specific']['max_uv_lambda']
        self.use_idg = self.field.parset['imaging_specific']['use_idg']
        self.idg_mode = self.field.parset['imaging_specific']['idg_mode']

        # Set image size based on current sector polygon
        xmin, ymin, xmax, ymax = self.poly.bounds
        self.width_ra = (xmax - xmin) * self.field.wcs_pixel_scale # deg
        self.width_dec = (ymax - ymin) * self.field.wcs_pixel_scale # deg
        self.imsize = [int(self.width_ra / self.cellsize_deg * 1.1),
                       int(self.width_dec / self.cellsize_deg * 1.1)]
        if self.use_idg:
            # IDG does not yet support rectangular images
            self.imsize = [max(self.imsize), max(self.imsize)]

            # IDG has problems with small images
            self.imsize = [max(1500, self.imsize[0]), max(1500, self.imsize[0])]

        self.wsclean_imsize = '{0} {1}'.format(self.imsize[0], self.imsize[1])
        self.log.debug('Image size is {0} x {1} pixels'.format(
                       self.imsize[0], self.imsize[1]))

        # Set number of output channels to get 5 channels, but of no less than ~ 2 MHz each
        target_nchannels = 5
        tot_bandwidth = 0.0
        for obs in self.observations:
            # Find observation with largest bandwidth
            obs_bandwidth = obs.numchannels * obs.channelwidth
            if obs_bandwidth > tot_bandwidth:
                tot_bandwidth = obs_bandwidth
        self.wsclean_nchannels = min(target_nchannels, int(np.ceil(tot_bandwidth / 2e6)))

        # Set number of iterations and threshold
        total_time_hr = 0.0
        for obs in self.observations:
            # Find total observation time in hours
            total_time_hr += (obs.endtime - obs.starttime) / 3600.0
        scaling_factor = np.sqrt(np.float(tot_bandwidth / 2e6) * total_time_hr / 8.0)
        self.wsclean_niter = int(12000 * scaling_factor)

        # Set multiscale: get source sizes and check for large sources
        self.multiscale = do_multiscale
        if self.multiscale is None:
            large_size_arcmin = 4.0 # threshold source size for multiscale to be activated
            sizes_arcmin = self.source_sizes * 60.0
            if sizes_arcmin is not None and any([s > large_size_arcmin for s in sizes_arcmin]):
                self.multiscale = True
            else:
                self.multiscale = False
        if self.multiscale:
            self.multiscale_scales_pixel = self.field.parset['imaging_specific']['multiscale_scales_pixel']
            self.wsclean_niter /= 2 # fewer iterations are needed
            self.log.debug("Will do multiscale cleaning.")
        else:
            self.multiscale_scales_pixel = 0

        # Set the observation-specific parameters
        max_peak_smearing = self.field.parset['imaging_specific']['max_peak_smearing']
        for obs in self.observations:
            # Set imaging parameters
            obs.set_imaging_parameters(self.cellsize_arcsec, max_peak_smearing,
                                       self.width_ra, self.width_dec)

        # Set BL-dependent averaging
        do_bl_averaging = False # does not yet work with IDG
        if do_bl_averaging:
            timestep_sec = (self.observations[0].timepersample *
                            self.observations[0].parameters['image_timestep'])
            self.wsclean_nwavelengths = self.get_nwavelengths(self.cellsize_deg,
                                                              timestep_sec)
        else:
            self.wsclean_nwavelengths = 0

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
        Makes predict sky model from the field calibration sky model
        """
        # Filter the predict sky model
        if self.is_outlier:
            # For outlier sector, we use the sky model made earlier
            skymodel = self.predict_skymodel
        else:
            skymodel = self.field.calibration_skymodel.copy()
            skymodel = self.filter_skymodel(skymodel)

        # Write filtered sky model to file
        self.predict_skymodel_file = os.path.join(self.field.working_dir, 'skymodels',
                                                      '{}_predict_skymodel.txt'.format(self.name))
        skymodel.write(self.predict_skymodel_file, clobber=True)

        # Save list of patches (directions) in the format written by DDECal in the h5parm
        self.patches = '[{}]'.format(','.join(['[{}]'.format(p) for p in skymodel.getPatchNames()]))

        # Find nearest patch to flux-weighted center of the sector sky model
        if not self.is_outlier:
            tmp_skymodel = skymodel.copy()
            tmp_skymodel.group('single')
            ra, dec = tmp_skymodel.getPatchPositions(method='wmean', asArray=True)
            patch_dist = skymodel.getDistance(ra[0], dec[0], byPatch=True).tolist()
            patch_names = skymodel.getPatchNames()
            self.central_patch = patch_names[patch_dist.index(min(patch_dist))]

            # Filter the field source sky model and store source sizes
            all_source_names = self.field.source_skymodel.getColValues('Name').tolist()
            source_names = skymodel.getColValues('Name')
            in_sector = np.array([all_source_names.index(sn) for sn in source_names])
            source_skymodel = self.field.source_skymodel.copy()
            source_skymodel.select(in_sector)
            self.source_sizes = source_skymodel.getPatchSizes(units='degree')

        # Set the parameters for predict
        self.set_predict_parameters()

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

        # Keep only those sources inside the sector bounding box
        inside = np.zeros(len(skymodel), dtype=bool)
        xmin, ymin, xmax, ymax = self.poly.bounds
        inside_ind = np.where( (x >= xmin) & (x <= xmax) & (y >= ymin ) & (y <= ymax ))
        inside[inside_ind] = True
        skymodel.select(inside)
        RA = skymodel.getColValues('Ra')
        Dec = skymodel.getColValues('Dec')
        x, y = self.field.radec2xy(RA, Dec)
        x = np.array(x)
        y = np.array(y)

        # Now check the actual sector boundary against filtered sky model
        xpadding = int(0.1 * (max(x) - min(x)))
        ypadding = int(0.1 * (max(y) - min(y)))
        xshift = int(min(x)) - xpadding
        yshift = int(min(y)) - ypadding
        xsize = int(np.ceil(max(x) - min(x))) + 2*xpadding
        ysize = int(np.ceil(max(y) - min(y))) + 2*ypadding
        x -= xshift
        y -= yshift
        prepared_polygon = prep(self.poly)

        # Unmask everything outside of the polygon + its border (outline)
        inside = np.zeros(len(skymodel), dtype=bool)
        mask = Image.new('L', (xsize, ysize), 0)
        verts = [(xv-xshift, yv-yshift) for xv, yv in zip(self.poly.exterior.coords.xy[0],
                                            self.poly.exterior.coords.xy[1])]
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=1)
        inside_ind = np.where(np.array(mask).transpose()[(x.astype(int), y.astype(int))])
        inside[inside_ind] = True

        # Now check sources in the border precisely
        mask = Image.new('L', (xsize, ysize), 0)
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=0)
        border_ind = np.where(np.array(mask).transpose()[(x.astype(int), y.astype(int))])
        points = [Point(xs, ys) for xs, ys in zip(x[border_ind], y[border_ind])]
        for i, p in enumerate(points):
            p.index = border_ind[0][i]
        outside_points = filter(lambda v: not prepared_polygon.contains(v), points)
        for outside_point in outside_points:
            inside[outside_point.index] = False
        skymodel.select(inside)
        return skymodel

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
            List of parameters, with one entry for each observation
        """
        return [obs.parameters[parameter] for obs in self.observations]

    def intialize_vertices(self):
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

