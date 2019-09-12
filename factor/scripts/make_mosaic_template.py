#! /usr/bin/env python
"""
Script to make a template image for mosaicking
"""
from factor.lib.image import FITSImage
from factor.lib import miscellaneous as misc
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
import numpy as np


def main(input_image_list, vertices_file_list, output_image, skip=False):
    """
    Make a mosaic template image

    Parameters
    ----------
    input_image_list : list
        List of filenames of input images to mosaic
    vertices_file_list : list
        List of filenames of input vertices files
    output_image : str
        Filename of output image
    skip : bool
        If True, skip all processing
    """
    input_image_list = misc.string2list(input_image_list)
    vertices_file_list = misc.string2list(vertices_file_list)
    skip = misc.string2bool(skip)
    if skip:
        return

    # Set up images used in mosaic
    directions = []
    for image_file, vertices_file in zip(input_image_list, vertices_file_list):
        d = FITSImage(image_file)
        d.vertices_file = vertices_file
        d.blank()
        directions.append(d)

    # Prepare header for final gridding
    mra = np.mean(np.array([d.get_wcs().wcs.crval[0] for d in directions]))
    mdec = np.mean(np.array([d.get_wcs().wcs.crval[1] for d in directions]))

    # Make a reference WCS and use it to find the extent in pixels
    # needed for the combined image
    rwcs = pywcs(naxis=2)
    rwcs.wcs.ctype = directions[0].get_wcs().wcs.ctype
    rwcs.wcs.cdelt = directions[0].get_wcs().wcs.cdelt
    rwcs.wcs.crval = [mra, mdec]
    rwcs.wcs.crpix = [1, 1]
    xmin, xmax, ymin, ymax = 0, 0, 0, 0
    for d in directions:
        w = d.get_wcs()
        ys, xs = np.where(d.img_data)
        axmin, aymin, axmax, aymax = xs.min(), ys.min(), xs.max(), ys.max()
        del xs, ys
        for x, y in ((axmin, aymin), (axmax, aymin), (axmin, aymax), (axmax, aymax)):
            ra, dec = [float(f) for f in w.wcs_pix2world(x, y, 0)]
            nx, ny = [float(f) for f in rwcs.wcs_world2pix(ra, dec, 0)]
            xmin, xmax, ymin, ymax = min(nx, xmin), max(nx, xmax), min(ny, ymin), max(ny, ymax)
    xsize = int(xmax-xmin)
    ysize = int(ymax-ymin)
    rwcs.wcs.crpix = [-int(xmin)+1, -int(ymin)+1]
    regrid_hdr = rwcs.to_header()
    regrid_hdr['NAXIS'] = 2
    regrid_hdr['NAXIS1'] = xsize
    regrid_hdr['NAXIS2'] = ysize
    for ch in ('BMAJ', 'BMIN', 'BPA'):
        regrid_hdr[ch] = pyfits.open(directions[0].imagefile)[0].header[ch]
    regrid_hdr['ORIGIN'] = 'Rapthor'
    regrid_hdr['UNITS'] = 'Jy/beam'
    isum = np.zeros([ysize, xsize])
    hdu = pyfits.PrimaryHDU(header=regrid_hdr, data=isum)
    hdu.writeto(output_image, overwrite=True)
