#! /usr/bin/env python
"""
Script to make a mosiac
"""
from factor.lib import miscellaneous as misc
from astropy.io import fits as pyfits
import numpy as np
import os


def main(input_image_list, template_image, output_image, skip=False):
    """
    Make a mosaic image

    Parameters
    ----------
    input_image_list : list
        List of filenames of input images to mosaic
    template_image : str
        Filename of mosaic template image
    output_image : str
        Filename of output image
    skip : bool
        If True, skip all processing
    """
    input_image_list = misc.string2list(input_image_list)
    skip = misc.string2bool(skip)
    if skip:
        os.system('cp {0} {1}'.format(input_image_list[0], output_image))
        return

    # Load template and sector images and add them to mosaic
    regrid_hdr = pyfits.open(template_image)[0].header
    isum = pyfits.open(template_image)[0].data
    wsum = np.zeros_like(isum)
    for sector_image in input_image_list:
        r = pyfits.open(sector_image)[0].data
        w = np.ones_like(r)
        w[np.isnan(r)] = 0
        r[np.isnan(r)] = 0
        isum += r
        wsum += w
    isum /= wsum
    del wsum, r, w
    isum[np.isnan(isum)] = np.nan
    hdu = pyfits.PrimaryHDU(header=regrid_hdr, data=isum)
    hdu.writeto(output_image, overwrite=True)
