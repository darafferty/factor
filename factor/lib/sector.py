"""
Definition of the Sector class that holds parameters for each iamge sector
set
"""
import logging
import numpy as np
from astropy.coordinates import Angle


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
    observations : list of Observation objects
        List of observations for the field
    """
    def __init__(self, name, ra, dec, width_ra, width_dec, observations):
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
        self.observations = observations[:]

    def set_imaging_parameters(self, cellsize_arcsec, robust, taper_arcsec,
                               min_uv_lambda, max_uv_lambda, max_peak_smearing):
        """
        Sets the imaging parameters for given values

        Parameters
        ----------
        field : Field object
            Field object
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
        self.robust = robust
        self.taper_arcsec = taper_arcsec
        self.min_uv_lambda = min_uv_lambda
        self.max_uv_lambda = max_uv_lambda

        # Set image size
        self.imsize = [self.width_ra*3600.0/self.cellsize_arcsec,
                       self.width_dec*3600.0/self.cellsize_arcsec]
        self.log.debug('Image size is {0} x {1} pixels'.format(
                       self.imsize[0], self.imsize[1]))

        # Set number of iterations and threshold
        scaling_factor = 2.0 * np.sqrt(iter+1)
        self.wsclean_niter = int(12000 * scaling_factor)

        # Set multiscale: get source sizes and check for large sources
        if self.contains_target:
            self.multiscale = True

        large_size_arcmin = 4.0
        if self.multiscale is None:
            sizes_arcmin = self.get_source_sizes()
            if sizes_arcmin is not None and any([s > large_size_arcmin for s in sizes_arcmin]):
                self.multiscale = True
            else:
                self.multiscale = False
        if self.multiscale:
            self.wsclean_multiscale = '-multiscale,'
            self.wsclean_full1_image_niter /= 2 # fewer iterations are needed
            self.wsclean_full2_image_niter /= 2 # fewer iterations are needed
            self.log.debug("Will do multiscale cleaning.")
        else:
            self.wsclean_multiscale = ''

        # Set the observation-specific parameters
        for obs in self.observations:
            obs.set_imaging_parameters(cellsize_arcsec, max_peak_smearing,
                                       self.width_ra, self.width_dec)

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


