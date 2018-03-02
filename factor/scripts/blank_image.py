#! /usr/bin/env python
"""
Script to blank regions (with zeros or NaNs) in a fits image
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os
import pickle
import glob
from shapely.geometry import Point, Polygon, box
from shapely.prepared import prep
from astropy.io import fits as pyfits
from astropy import wcs
from PIL import Image, ImageDraw


def read_vertices(filename):
    """
    Returns facet vertices stored in input file
    """
    with open(filename, 'r') as f:
        vertices = pickle.load(f)
    return vertices


def main(input_image_file, vertices_file, output_image_file, blank_value='zero',
    image_is_wsclean_model=False):
    """
    Blank a region in an image

    Parameters
    ----------
    input_image_file : str
        Filename of input image to blank
    vertices_file : str, optional
        Filename of file with vertices (must be a pickle file containing
        a dictionary with the vertices in the 'vertices' entry)
    output_image_file : str
        Filename of output image
    blank_value : str, optional
        Value for blanks (one of 'zero' or 'nan')
    image_is_wsclean_model : bool, optional
        If True, the input and output image files are treated as the root name
        of a WSClean model image (or images)

    """
    if type(image_is_wsclean_model) is str:
        if image_is_wsclean_model.lower() == 'true':
            image_is_wsclean_model = True
        else:
            image_is_wsclean_model = False

    if image_is_wsclean_model:
        input_image_files = glob.glob(input_image_file+'*-model.fits')
        output_image_files = [f.replace(input_image_file, output_image_file) for f in input_image_files]
    else:
        input_image_files = [input_image_file]
        output_image_files = [output_image_file]

    if blank_value == 'zero':
        blank_val = 0.0
    elif blank_value == 'nan':
        blank_val = np.nan
    else:
        print('Blank value type "{}" not understood.'.format(blank_with))
        sys.exit(1)

    # Construct polygon of facet region
    header = pyfits.getheader(input_image_files[0], 0)
    w = wcs.WCS(header)
    RAind = w.axis_type_names.index('RA')
    Decind = w.axis_type_names.index('DEC')
    vertices = read_vertices(vertices_file)
    RAverts = vertices[0]
    Decverts = vertices[1]
    verts = []
    for RAvert, Decvert in zip(RAverts, Decverts):
        ra_dec = np.array([[0.0, 0.0, 0.0, 0.0]])
        ra_dec[0][RAind] = RAvert
        ra_dec[0][Decind] = Decvert
        verts.append((w.wcs_world2pix(ra_dec, 0)[0][RAind], w.wcs_world2pix(ra_dec, 0)[0][Decind]))
    poly = Polygon(verts)
    prepared_polygon = prep(poly)

    for input_image, output_image in zip(input_image_files, output_image_files):
        hdu = pyfits.open(input_image, memmap=False)
        data = hdu[0].data

        # Mask everything outside of the polygon + plus its border (outline)
        mask = Image.new('L', (data.shape[2], data.shape[3]), 0)
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=1)
        data[0, 0, :, :] *= mask

        # Now check the border precisely
        mask = Image.new('L', (data.shape[2], data.shape[3]), 0)
        ImageDraw.Draw(mask).polygon(verts, outline=1, fill=0)
        masked_ind = np.where(np.array(mask).transpose())
        points = [Point(xm, ym) for xm, ym in zip(masked_ind[0], masked_ind[1])]
        outside_points = filter(lambda v: not prepared_polygon.contains(v), points)
        for outside_point in outside_points:
            data[0, 0, int(outside_point.y), int(outside_point.x)] = 0

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
