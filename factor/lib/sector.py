"""
Definition of the Sector class that holds parameters for each iamge sector
set
"""
import logging
import numpy as np
from astropy.coordinates import Angle
from shapely.geometry import Point, Polygon
from shapely.prepared import prep
from shapely.wkt import dumps
import lsmtool
import os


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
        self.observations = field.observations[:]
        self.vertices_file = os.path.join(field.working_dir, 'regions', '{}_vertices.txt'.format(self.name))
        self.region_file = ''

        # Define the sector polygon vertices and sky model
        self.define_vertices()
        self.make_skymodel()

    def set_imaging_parameters(self, cellsize_arcsec, robust, taper_arcsec,
                               min_uv_lambda, max_uv_lambda, max_peak_smearing):
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
        """
        self.cellsize_arcsec = cellsize_arcsec
        self.cellsize_deg = cellsize_arcsec / 3600.0
        self.robust = robust
        self.taper_arcsec = taper_arcsec
        self.min_uv_lambda = min_uv_lambda
        self.max_uv_lambda = max_uv_lambda

        # Set image size
        self.imsize = [self.width_ra / self.cellsize_deg,
                       self.width_dec / self.cellsize_deg]
        self.log.debug('Image size is {0} x {1} pixels'.format(
                       self.imsize[0], self.imsize[1]))

        # Set number of iterations and threshold
        scaling_factor = 2.0 * np.sqrt(iter+1)
        self.wsclean_niter = int(12000 * scaling_factor)

        # Set number of output channels to get ~ 4 MHz per channel
        tot_bandwidth = 0.0
        for obs in self.observations:
            # Find largest bandwidth
            obs_bandwidth = obs.numchannels * obs.channelwidth
            if obs_bandwidth > tot_bandwidth:
                tot_bandwidth = obs_bandwidth
        self.wsclean_nchannels = max(1, int(np.ceil(tot_bandwidth / 4e6)))

        # Set multiscale: get source sizes and check for large sources
        self.multiscale = None
        large_size_arcmin = 4.0
        if self.multiscale is None:
            sizes_arcmin = self.get_source_sizes()
            if sizes_arcmin is not None and any([s > large_size_arcmin for s in sizes_arcmin]):
                self.multiscale = True
            else:
                self.multiscale = False
        if self.multiscale:
            self.wsclean_multiscale = '-multiscale,'
            self.wsclean_niter /= 2 # fewer iterations are needed
            self.log.debug("Will do multiscale cleaning.")
        else:
            self.wsclean_multiscale = ''

        # Set the observation-specific parameters
        for obs in self.observations:
            # Set filename for model-subtracted data that matches the one made by the
            # calibrate pipeline
            ms_subtracted_filename = '{0}.sector_{1}'.format(obs.ms_filename,
                                                            self.name.split('_')[1])
            # Set imaging parameters
            obs.set_imaging_parameters(cellsize_arcsec, max_peak_smearing,
                                       self.width_ra, self.width_dec, ms_subtracted_filename)

    def make_skymodel(self):
        """
        Makes a sky model for the sector from the parent field sky model
        """
        skymodel = lsmtool.load(self.field.skymodel_file)

        # Make list of sources in full sky model
        RA = skymodel.getColValues('Ra')
        Dec = skymodel.getColValues('Dec')
        x, y = self.field.radec2xy(RA, Dec)
        points = []
        for i, (xp, yp) in enumerate(zip(x, y)):
            p = Point((xp, yp))
            p.index = i
            points.append(p)

        # Find sources that are inside the sector
        vertices = self.get_vertices_radec()
        xv, yv = self.field.radec2xy(vertices[0], vertices[1])
        poly = Polygon([(xp, yp) for xp, yp in zip(xv, yv)])
        prepared_polygon = prep(poly)
        intersecting_points = filter(prepared_polygon.contains, points)
        inside = np.zeros(len(skymodel), dtype=bool)
        for p in intersecting_points:
            inside[p.index] = True
        skymodel.select(inside, force=True)

        # Write sky model to file
        skymodel.setPatchPositions(method='wmean')
        self.skymodel_file = os.path.join(self.field.working_dir, 'skymodels', '{}_skymodel.txt'.format(self.name))
        skymodel.write(self.skymodel_file, clobber=True)

        # Save list of patches (directions) in the format written by DDECal in the h5parm
        self.patches = '[{}]'.format(','.join(['[{}]'.format(p) for p in skymodel.getPatchNames()]))

        # Find nearest patch to sector center
        patch_dist = skymodel.getDistance(self.ra, self.dec, byPatch=True).tolist()
        patch_names = skymodel.getPatchNames()
        self.central_patch = patch_names[patch_dist.index(min(patch_dist))]

    def get_source_sizes(self):
        """
        Returns list of source sizes in arcmin
        """
        skymodel = lsmtool.load(self.skymodel_file)
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

        # Find nearby sources in input sky model and adjust sector boundaries
        # if necessary
        lsmtool._logging.setLevel('debug')
        s = lsmtool.load(self.field.skymodel_file)
        dists = s.getDistance(self.ra, self.dec, byPatch=False)
        radius = np.hypot(self.width_ra, self.width_dec)
        s.select(dists < radius)
        if len(s) > 0:
            # Group components in input sky model by thresholding after
            # convolving model with 1-arcmin beam
            s.group('threshold', FWHM='60.0 arcsec', root='facet', threshold=0.01)
            s.remove('Patch = patch_*', force=True) # Remove sources that did not threshold
            RA, Dec = s.getPatchPositions(asArray=True)
            sx, sy = self.field.radec2xy(RA, Dec)
            sizes = s.getPatchSizes(units='degree').tolist()

            # Make buffered points for all sources in the input sky model
            points = []
            for xp, yp, sp in zip(sx, sy, sizes):
                radius = sp * 1.2 / 2.0 / self.field.wcs.wcs.cdelt[0] # size of source in pixels
                points.append(Point((xp, yp)).buffer(sp))

            # Alter sector polygon to avoid sources in the input sky model
            niter = 0
            while niter < 3:
                niter += 1

                prepared_polygon = prep(poly)
                intersecting_points = filter(prepared_polygon.intersects, points)

                # Adjust sector polygon for each source that intersects it
                for p2 in intersecting_points:
                    if poly.contains(p2.centroid):
                        # If centroid of point is outside, difference the polys
                        poly = poly.difference(p2)
                    else:
                        # If point is inside, union the polys
                        poly = poly.union(p2)
        self.poly = poly

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
        with open(self.vertices_file, 'wb') as f:
            dumps(self.poly)

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

