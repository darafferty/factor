#! /usr/bin/env python
"""
Script to regrid a FITS image
"""
from factor.lib import miscellaneous as misc
from factor.lib.image import FITSImage
from reproject import reproject_interp
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
import numpy as np


def main(input_image, template_image, vertices_file, output_image, skip=False):
    """
    Regrid a FITS image

    Parameters
    ----------
    input_image : str
        Filename of input FITS image to blank
    template_image : str
        Filename of mosaic template FITS image
    vertices_file : str
        Filename of file with vertices
    output_image : str
        Filename of output FITS image
    skip : bool
        If True, skip all processing
    """
    skip = misc.string2bool(skip)
    if skip:
        return

    # Read template header and data
    regrid_hdr = pyfits.open(template_image)[0].header
    isum = pyfits.open(template_image)[0].data
    isum[:] = np.nan
    shape_out = isum.shape
    wcs_out = pywcs(regrid_hdr)

    # Read input image and blank outside its polygon
    d = FITSImage(input_image)
    d.vertices_file = vertices_file
    d.blank()
    wcs_in = d.get_wcs()

    # Define the subarray of the output image that fully encloses the reprojected input
    # image
    ny, nx = d.img_data.shape
    xc = np.array([-0.5, nx - 0.5, nx - 0.5, -0.5])
    yc = np.array([-0.5, -0.5, ny - 0.5, ny - 0.5])
    xc_out, yc_out = wcs_out.world_to_pixel(wcs_in.pixel_to_world(xc, yc))
    imin = max(0, int(np.floor(xc_out.min() + 0.5)))
    imax = min(shape_out[1], int(np.ceil(xc_out.max() + 0.5)))
    jmin = max(0, int(np.floor(yc_out.min() + 0.5)))
    jmax = min(shape_out[0], int(np.ceil(yc_out.max() + 0.5)))

    # Set up output projection
    wcs_out_indiv = wcs_out.deepcopy()
    wcs_out_indiv.wcs.crpix[0] -= imin
    wcs_out_indiv.wcs.crpix[1] -= jmin
    shape_out_indiv = (jmax - jmin, imax - imin)

    # Reproject, place into output image, and write out final FITS file
    ind = slice(jmin, jmax), slice(imin, imax)
    isum[ind] = reproject_interp((d.img_data, wcs_in), output_projection=wcs_out_indiv,
                                 shape_out=shape_out_indiv, return_footprint=False)
    d.img_data = isum
    d.img_hdr = regrid_hdr
    d.write(output_image)
