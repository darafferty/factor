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
    xslices = [slice(0, int(isum.shape[0] / 2.0)),
               slice(int(isum.shape[0] / 2.0), isum.shape[0])]
    yslices = [slice(0, int(isum.shape[1] / 2.0)),
               slice(int(isum.shape[1] / 2.0), isum.shape[1])]
    for xs in xslices:
        for ys in yslices:
            wsum = np.zeros_like(isum[xs, ys])
            for sector_image in input_image_list:
                r = pyfits.open(sector_image)[0].section[xs, ys]
                w = np.ones_like(r)
                w[np.isnan(r)] = 0
                r[np.isnan(r)] = 0
                isum[xs, ys] += r
                wsum += w
            isum[xs, ys] /= wsum
    del wsum, r, w
    isum[np.isnan(isum)] = np.nan
    hdu = pyfits.PrimaryHDU(header=regrid_hdr, data=isum)
    hdu.writeto(output_image, overwrite=True)
