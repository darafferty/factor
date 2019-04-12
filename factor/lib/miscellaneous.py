"""
Miscellaneous functions
"""
import numpy as np
import sys
import pickle
from shapely.geometry import Point, Polygon
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


def make_template_image(image_name, reference_ra_deg, reference_dec_deg,
    ximsize=512, yimsize=512, cellsize_deg=0.000417, freqs=None, times=None,
    antennas=None, aterm_type='tec', fill_val=0):
    """
    Make a blank image and save it to disk

    Parameters
    ----------
    image_name : str
        Filename of output image
    reference_ra_deg : float, optional
        RA for center of output mask image
    reference_dec_deg : float, optional
        Dec for center of output mask image
    imsize : int, optional
        Size of output image
    cellsize_deg : float, optional
        Size of a pixel in degrees
    freqs : list
        Frequencies to use to construct extra axes (for IDG a-term images)
    times : list
        Times to use to construct extra axes (for IDG a-term images)
    antennas : list
        Antennas to use to construct extra axes (for IDG a-term images)
    aterm_type : str
        One of 'tec' or 'gain'
    fill_val : int
        Value with which to fill the data
    """
    if freqs is not None and times is not None and antennas is not None:
        if aterm_type == 'tec':
            # TEC solutions
            # data is [RA, DEC, ANTENNA, FREQ, TIME].T
            nants = len(antennas)
            ntimes = len(times)
            nfreqs = len(freqs)
            shape_out = [ntimes, nfreqs, nants, yimsize, ximsize]
            ref_freq = freqs[0]
        else:
            # Gain solutions
            # data is [RA, DEC, MATRIX, ANTENNA, FREQ, TIME].T
            nants = len(antennas)
            ntimes = len(times)
            nfreqs = 1  # IDG does not support freq yet for gains
            shape_out = [ntimes, nfreqs, nants, 4, yimsize, ximsize]
            ref_freq = np.mean(freqs)
    else:
        # Normal FITS image
        # data is [STOKES, FREQ, DEC, RA]
        shape_out = [1, 1, yimsize, ximsize]
        ref_freq = 150e6

    hdu = pyfits.PrimaryHDU(np.ones(shape_out, dtype=np.float32)*fill_val)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header

    # Add RA, Dec info
    header['CRVAL1'] = reference_ra_deg
    header['CDELT1'] = -cellsize_deg
    header['CRPIX1'] = ximsize / 2.0
    header['CUNIT1'] = 'deg'
    header['CTYPE1'] = 'RA---SIN'
    header['CRVAL2'] = reference_dec_deg
    header['CDELT2'] = cellsize_deg
    header['CRPIX2'] = yimsize / 2.0
    header['CUNIT2'] = 'deg'
    header['CTYPE2'] = 'DEC--SIN'

    # Add STOKES info or ANTENNA info
    if antennas is None:
        header['CRVAL3'] = 1.0
        header['CDELT3'] = 1.0
        header['CRPIX3'] = 1.0
        header['CUNIT3'] = ''
        header['CTYPE3'] = 'STOKES'
    else:
        header['CRVAL3'] = 0.0
        header['CDELT3'] = 1.0
        header['CRPIX3'] = 1.0
        header['CUNIT3'] = ''
        header['CTYPE3'] = 'ANTENNA'

    # Add frequency info
    header['RESTFRQ'] = ref_freq
    header['CRVAL4'] = ref_freq
    header['CDELT4'] = 3e8
    header['CRPIX4'] = 1.0
    header['CUNIT4'] = 'Hz'
    header['CTYPE4'] = 'FREQ'

    # Add time info
    if times is not None:
        ref_time = times[0]
        if ntimes > 1:
            deltas = times[1:] - times[:-1]
            del_time = np.min(deltas)
        else:
            del_time = 1.0
        header['CRVAL5'] = ref_time
        header['CDELT5'] = del_time
        header['CRPIX5'] = 1.0
        header['CUNIT5'] = 's'
        header['CTYPE5'] = 'TIME'

    # Add matrix info
    if aterm_type == 'gain':
        header['CRVAL6'] = 1.0
        header['CDELT6'] = 1.0
        header['CRPIX6'] = 1.0
        header['CUNIT6'] = ''
        header['CTYPE6'] = 'MATRIX'

    # Add equinox
    header['EQUINOX'] = 2000.0

    # Add telescope
    header['TELESCOP'] = 'LOFAR'

    hdulist[0].header = header
    hdulist.writeto(image_name, overwrite=True)
    hdulist.close()


def rasterize(verts, data):
    """
    Rasterize a polygon into a data array

    Parameters
    ----------
    verts : list of (x, y) tuples
        List of input vertices of polygon to rasterize
    data : 2-D array
        Array into which rasterize polygon
    """
    poly = Polygon(verts)
    prepared_polygon = prep(poly)

    # Mask everything outside of the polygon plus its border (outline) with zeros
    # (inside polygon plus border are ones)
    mask = Image.new('L', (data.shape[0], data.shape[1]), 0)
    ImageDraw.Draw(mask).polygon(verts, outline=1, fill=1)
    data *= mask

    # Now check the border precisely
    mask = Image.new('L', (data.shape[0], data.shape[1]), 0)
    ImageDraw.Draw(mask).polygon(verts, outline=1, fill=0)
    masked_ind = np.where(np.array(mask).transpose())
    points = [Point(xm, ym) for xm, ym in zip(masked_ind[0], masked_ind[1])]
    outside_points = filter(lambda v: prepared_polygon.disjoint(v), points)
    for outside_point in outside_points:
        data[int(outside_point.y), int(outside_point.x)] = 0

    return data
