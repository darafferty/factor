#! /usr/bin/env python
"""
Script to create a mosaic from facet images
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables
import casacore.images as pim
from casacore import quanta
import numpy as np
from astropy.io import fits as pyfits
import os
import itertools
import multiprocessing
from factor.cluster import get_total_memory


def regrid_star(inputs):
    """
    Simple helper function for regrid
    """
    return regrid(*inputs)


def regrid(im, ind, ma, outshape):
    """
    Creates mosaic

    Parameters
    ----------
    images : str or list of str
        List of filenames of facet images. May be given as a list or as a string
        (e.g., '[image1, image2]'. Each image must be blanked with zeros outside
        of the facet region
    outfits : str
        Filename of output FITS mosaic
    maxwidth : int, optional
        Maximum number of pixels to consider for the width of the mosaic
        [default 0 = unlimited] This can be helpful at high declination.
    """
    image = pim.image(im)

    # Get Stokes axis. Ensure we are working with the Stokes parameter requested.
    sptcoords = image.coordinates().get_coordinate('spectral')
    nc = sptcoords.get_axis_size()
    stkcoords = image.coordinates().get_coordinate('stokes')
    if stkcoords.get_axis_size() == 1:
        assert(stkcoords.get_stokes()[0] == 'I')
    else:
        stks = stkcoords.get_stokes().index('I')
        image = image.subimage(blc=(0, stks), trc=(nc-1, stks), dropdegenerate=False)

    return image.regrid(ind, ma, outshape=outshape).getdata()


def main(images, outfits, maxwidth=0):
    """
    Creates mosaic

    Parameters
    ----------
    images : str or list of str
        List of filenames of facet images. May be given as a list or as a string
        (e.g., '[image1, image2]'. Each image must be blanked with zeros outside
        of the facet region
    outfits : str
        Filename of output FITS mosaic
    maxwidth : int, optional
        Maximum number of pixels to consider for the width of the mosaic
        [default 0 = unlimited] This can be helpful at high declination.
    """
    if type(images) is str:
        images = images.strip('[]').split(',')
        images = [im.strip() for im in images]

    psf_fwhm = [] # resolution
    frequency = [] # frequency of images (should be equal?)
    for i in range(len(images)):
        this_pim = pim.image(images[i])
        info_dict = this_pim.info()['imageinfo']['restoringbeam']
        bpar_ma = quanta.quantity(info_dict['major']).get_value('deg')
        bpar_mi = quanta.quantity(info_dict['minor']).get_value('deg')
        bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
        psf_fwhm.append([bpar_ma, bpar_mi, bpar_pa])
        frequency.append(this_pim.info()['coordinates']['spectral2']['restfreq'])

    psf_fwhm = np.array(psf_fwhm)
    frequency = np.array(frequency)
    mean_psf_fwhm = np.mean(psf_fwhm, axis=0)
    mean_frequency = np.mean(frequency)

    # Initialize some vectors
    declims = [] # store the limits of the declination axes
    raleft = []
    raright = []
    rainc = [] # store the r.a. increments in case they differ
    decinc = [] # store the dec increments in case they differ
    pims = [] # stores the casacore images of the data

    # Get image frames for input images
    for im in images:
        image = pim.image(im)
        sptcoords = image.coordinates().get_coordinate('spectral')
        nc = sptcoords.get_axis_size()

        # Get Stokes axis. Ensure we are working with the Stokes parameter requested.
        stkcoords = image.coordinates().get_coordinate('stokes')
        if stkcoords.get_axis_size() == 1:
            assert(stkcoords.get_stokes()[0] == 'I')
        else:
            stks = stkcoords.get_stokes().index('I')
            image = image.subimage(blc=(0, stks), trc=(nc-1, stks), dropdegenerate=False)
        ns = 1

        dircoords = image.coordinates().get_coordinate('direction')
        nx = dircoords.get_axis_size(axis=1)
        ny = dircoords.get_axis_size(axis=0)
        inc = dircoords.get_increment()
        ref = dircoords.get_referencepixel()
        val = dircoords.get_referencevalue()
        # wsclean image header is weird
        if val[1]<0:
            val[1]+=2*np.pi
        ra_axis = (range(nx)-ref[1])*inc[1]+val[1]
        dec_axis = (range(ny)-ref[0])*inc[0]+val[0]
        rainc.append(inc[1])
        decinc.append(inc[0])
        declims.append(min(dec_axis))
        declims.append(max(dec_axis))
        mean_ra = np.mean(ra_axis)
        raleft.append((ra_axis[0] - mean_ra) * np.cos(val[0]) + mean_ra)
        raright.append(mean_ra - (mean_ra - ra_axis[-1]) * np.cos(val[0]))
        pims.append(image)

    # Generate the mosaic coordinate frame
    master_dec = np.arange(min(declims), max(declims), min(decinc))
    if max(raleft)-min(raright) > 5.*np.pi/3.: # crossed RA=0
        for i in range(len(raright)):
            raright[i] = raright[i]-2.*np.pi
    master_ra = np.arange(max(raleft), min(raright), max(rainc))
    lmra = len(master_ra)
    if maxwidth != 0:
        if lmra > maxwidth:
            xboundary = (lmra-maxwidth)/2
            master_ra = master_ra[xboundary:-xboundary]
    ma = pims[-1].coordinates()
    ma['direction'].set_referencepixel([len(master_dec)/2, len(master_ra)/2])
    ma['direction'].set_increment([decinc[np.argmin(np.abs(decinc))], rainc[np.argmin(np.abs(rainc))]])
    ma['direction'].set_referencevalue([master_dec[len(master_dec)/2], master_ra[len(master_ra)/2]])

    # Initialize the output image
    master_im = np.zeros((len(master_dec), len(master_ra)))
    mem_req_gb = master_im.nbytes / 1024**3 * 6  # each regrid takes ~ 6 x master_im size
    mem_avil_gb = get_total_memory()
    nimg_simul = int(mem_avil_gb / mem_req_gb) - 6

    # Reproject the images onto the master grid
    ind = [2,3]
    outshape = (int(nc), int(ns), len(master_dec), len(master_ra))
    pool = multiprocessing.Pool(nimg_simul)
    results = pool.map(regrid_star, itertools.izip(images, itertools.repeat(ind),
                                                   itertools.repeat(ma),
                                                   itertools.repeat(outshape)))
    pool.close()
    pool.join()
    for imdata in results:
        master_im += np.squeeze(imdata)
    del(results)
    blank = np.ones_like(master_im) * np.nan
    master_im = np.where(master_im, master_im, blank)

    # Write output
    arrax = np.zeros( (1,1, len(master_im[:,0]), len(master_im[0,:])) )
    arrax[0,0,:,:] = master_im
    new_pim = pim.image('', shape=(1, 1, len(master_dec), len(master_ra)), coordsys=ma)
    new_pim.putdata(arrax)
    new_pim.tofits(outfits, overwrite=True)

    # need to add new beam info (not sure if this is possible with casacore)
    hdu = pyfits.open(outfits, mode='update', memmap=False)
    header = hdu[0].header
    header['BMAJ'] = mean_psf_fwhm[0]
    header['BMIN'] = mean_psf_fwhm[1]
    header['BPA'] = mean_psf_fwhm[2]
    header['BUNIT'] = pims[-1].info()['unit']
    header['RESTFRQ'] = mean_frequency
    header['RESTFREQ'] = mean_frequency
    newhdu = pyfits.PrimaryHDU(data=hdu[0].data, header=header)
    newhdu.writeto(outfits, overwrite=True)


if __name__ == '__main__':
    descriptiontext = "Create a mosaic from facet images.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('images', help='List of filenames of facet images')
    parser.add_argument('outfits', help='Output name of mosaic fits file')
    parser.add_argument('-m','--maxwidth', help='Maximum number of pixels to '
        'consider for the width of the mosaic [default 0 = unlimited] This can '
        'be helpful at high declination.', default=0, type=int)

    args = parser.parse_args()
    main(args.images, args.outfits, maxwidth=args.maxwidth)
