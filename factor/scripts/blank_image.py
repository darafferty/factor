#! /usr/bin/env python
"""
Script to blank regions (with zeros or NaNs) in a fits image. Can also be used to make
a clean mask
"""
import argparse
from argparse import RawTextHelpFormatter
from .lib import miscellaneous as misc
import numpy as np
import sys
from astropy.io import fits as pyfits
from astropy import wcs


def main(input_image, output_image, vertices_file=None, reference_ra_deg=None,
         reference_dec_deg=None, cellsize_deg=None, imsize=None, make_blank_image=False,
         region_file=None):
    """
    Blank a region in an image

    Parameters
    ----------
    input_image : str
        Filename of input image to blank
    output_image : str
        Filename of output image
    vertices_file : str, optional
        Filename of file with vertices
    reference_ra_deg : float, optional
        RA for center of output mask image
    reference_dec_deg : float, optional
        Dec for center of output mask image
    cellsize_deg : float, optional
        Size of a pixel in degrees
    imsize : int, optional
        Size of image as "xsize ysize"
    region_file : str, optional
        Filename of region file in CASA format to use as a mask
    make_blank_image : bool, optional
        If True, a blank template image is made. In this case, reference_ra_deg
        and reference_dec_deg must be specified
    region_file : str, optional
        Filename of region file in CASA format to use as the mask
    """
    if type(make_blank_image) is str:
        if make_blank_image.lower() == 'true':
            make_blank_image = True
        else:
            make_blank_image = False

    if make_blank_image:
        print('Making empty image...')
        if reference_ra_deg is not None and reference_dec_deg is not None:
            reference_ra_deg = float(reference_ra_deg)
            reference_dec_deg = float(reference_dec_deg)
            temp_image = output_image + '.tmp'
            ximsize = int(imsize.split(' ')[0])
            yimsize = int(imsize.split(' ')[1])
            misc.make_template_image(temp_image, reference_ra_deg, reference_dec_deg,
                                     ximsize=ximsize, yimsize=yimsize,
                                     cellsize_deg=float(cellsize_deg))
        else:
            print('ERROR: a reference position must be given to make an empty template image')
            sys.exit(1)

    if vertices_file is not None:
        # Construct polygon
        if make_blank_image:
            header = pyfits.getheader(temp_image, 0)
        else:
            header = pyfits.getheader(input_image, 0)
        w = wcs.WCS(header)
        RAind = w.axis_type_names.index('RA')
        Decind = w.axis_type_names.index('DEC')
        vertices = misc.read_vertices(vertices_file)
        RAverts = vertices[0]
        Decverts = vertices[1]
        verts = []
        for RAvert, Decvert in zip(RAverts, Decverts):
            ra_dec = np.array([[0.0, 0.0, 0.0, 0.0]])
            ra_dec[0][RAind] = RAvert
            ra_dec[0][Decind] = Decvert
            verts.append((w.wcs_world2pix(ra_dec, 0)[0][RAind], w.wcs_world2pix(ra_dec, 0)[0][Decind]))

        if make_blank_image:
            hdu = pyfits.open(temp_image, memmap=False)
        else:
            hdu = pyfits.open(input_image, memmap=False)
        data = hdu[0].data

        # Transpose the array (since FITS standard is [0, 0, y, x]), rasterize the poly,
        # the transpose back
        data_rasertize = data[0, 0, :, :].T
        data_rasertize = misc.rasterize(verts, data_rasertize)
        data[0, 0, :, :] = data_rasertize.T

        hdu[0].data = data
        hdu.writeto(output_image, overwrite=True)


if __name__ == '__main__':
    descriptiontext = "Blank regions of an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_image_file', help='Filename of input image')
    parser.add_argument('vertices_file', help='Filename of vertices file')
    parser.add_argument('output_image_file', help='Filename of output image')
    parser.add_argument('image_is_wsclean_model', help='True if input is WSClean model root', type=bool, default=False)
    parser.add_argument('-b', '--blank_value', help='value for blank pixesl', type=str, default='zero')
    args = parser.parse_args()
    main(args.input_image_file, args.vertices_file, args.output_image_file,
         blank_value=args.blank_value, image_is_wsclean_model=args.image_is_wsclean_model)
