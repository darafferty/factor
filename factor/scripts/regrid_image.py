#! /usr/bin/env python
"""
Script to regrid a FITS image
"""
from factor.lib import miscellaneous as misc
from factor.lib.image import FITSImage
from reproject import reproject_interp
from astropy.io import fits as pyfits


def main(input_image, template_image, vertices_file, output_image, do_weights=False,
         skip=False):
    """
    Regrid a FITS image

    Parameters
    ----------
    input_image : str
        Filename of input image to blank
    template_image : str
        Filename of output image
    vertices_file : str
        Filename of file with vertices
    output_image : str
        Filename of template mosaic image
    do_weights : bool, optional
        If True, save weights
    skip : bool
        If True, skip all processing
    """
    do_weights = misc.string2bool(do_weights)
    skip = misc.string2bool(skip)
    if skip:
        return

    # Read template header
    regrid_hdr = pyfits.open(template_image)[0].header
    d = FITSImage(input_image)
    d.vertices_file = vertices_file
    d.blank()
    d.calc_weight()
    r, footprint = reproject_interp((d.img_data, d.img_hdr), regrid_hdr)
    d.img_data = r
    d.img_hdr = regrid_hdr
    d.write(output_image)
    if do_weights:
        w, footprint = reproject_interp((d.weight_data, d.img_hdr), regrid_hdr)
        d.img_data = w
        d.img_hdr = regrid_hdr
        d.write(output_image+'weights')
