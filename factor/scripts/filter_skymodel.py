#! /usr/bin/env python
"""
Script to filter a sky model with an image
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import numpy as np
# workaround for a bug / ugly behavior in matplotlib / pybdsf
# (some installtaions of matplotlib.pyplot fail if there is no X-Server available
# which is not caught in pybdsf.)
try:
    import matplotlib
    matplotlib.use('Agg')
except (RuntimeError, ImportError):
    pass
try:
    import bdsf
except ImportError:
    from lofar import bdsm as bdsf


def main(input_image, input_skymodel, output_skymodel, threshisl=3.0, threshpix=5.0,
         rmsbox=(150, 50), rmsbox_bright=(35, 7), adaptive_rmsbox=True,
         use_adaptive_threshold=False, adaptive_thresh=150.0):
    """
    Filter the input sky model so that they lie in islands in the image

    Parameters
    ----------
    input_image : str
        Filename of input image to blank
    input_skymodel : str
        Filename of input makesourcedb sky model
    output_skymodel : str
        Filename of output makesourcedb sky model
    threshisl : float, optional
        Value of thresh_isl PyBDSF parameter
    threshpix : float, optional
        Value of thresh_pix PyBDSF parameter
    rmsbox : tuple of floats, optional
        Value of rms_box PyBDSF parameter
    rmsbox_bright : tuple of floats, optional
        Value of rms_box_bright PyBDSF parameter
    adaptive_rmsbox : tuple of floats, optional
        Value of adaptive_rms_box PyBDSF parameter
    use_adaptive_threshold : bool, optional
        If True, use an adaptive threshold estimated from the negative values in
        the image
    adaptive_thresh : float, optional
        If adaptive_rmsbox is True, this value sets the threshold above
        which a source will use the small rms box
    """
    if rmsbox is not None and isinstance(rmsbox, basestring):
        rmsbox = eval(rmsbox)

    if isinstance(rmsbox_bright, basestring):
        rmsbox_bright = eval(rmsbox_bright)

    if isinstance(adaptive_rmsbox, basestring):
        if adaptive_rmsbox.lower() == 'true':
            adaptive_rmsbox = True
        else:
            adaptive_rmsbox = False

    if isinstance(use_adaptive_threshold, basestring):
        if use_adaptive_threshold.lower() == 'true':
            use_adaptive_threshold = True
        else:
            use_adaptive_threshold = False

    if use_adaptive_threshold:
        # Get an estimate of the rms
        img = bdsf.process_image(input_image, mean_map='zero', rms_box=rmsbox,
                                 thresh_pix=threshpix, thresh_isl=threshisl,
                                 thresh='hard', adaptive_rms_box=adaptive_rmsbox,
                                 adaptive_thresh=adaptive_thresh, rms_box_bright=rmsbox_bright,
                                 rms_map=True, quiet=True, stop_at='isl')

        # Find min and max pixels
        max_neg_val = abs(np.min(img.ch0_arr))
        max_neg_pos = np.where(img.ch0_arr == np.min(img.ch0_arr))
        max_pos_val = abs(np.max(img.ch0_arr))
        max_pos_pos = np.where(img.ch0_arr == np.max(img.ch0_arr))

        # Estimate new thresh_isl from min pixel value's sigma, but don't let
        # it get higher than 1/2 of the peak's sigma
        threshisl_neg = 2.0 * max_neg_val / img.rms_arr[max_neg_pos][0]
        max_sigma = max_pos_val / img.rms_arr[max_pos_pos][0]
        if threshisl_neg > max_sigma / 2.0:
            threshisl_neg = max_sigma / 2.0

        # Use the new threshold only if it is larger than the user-specified one
        if threshisl_neg > threshisl:
            threshisl = threshisl_neg

    img = bdsf.process_image(input_image, mean_map='zero', rms_box=rmsbox,
                             thresh_pix=threshpix, thresh_isl=threshisl,
                             thresh='hard', adaptive_rms_box=adaptive_rmsbox,
                             adaptive_thresh=adaptive_thresh, rms_box_bright=rmsbox_bright,
                             rms_map=True, quiet=True, stop_at='isl')

    if img.nisl > 0:
        maskfile = input_image + '.mask'
        img.export_image(outfile=maskfile, clobber=True, img_type='island_mask')

        s = lsmtool.load(input_skymodel)
        s.select('{} == True'.format(maskfile))  # keep only those in PyBDSF masked regions
        s.group(maskfile)  # group the sky model by mask islands
        s.write(output_skymodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Filter a sky model with an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_image', help='Filename of input image')
    parser.add_argument('input_skymodel', help='Filename of input sky model')
    parser.add_argument('output_skymodel', help='Filename of output sky model')
    parser.add_argument('-i', '--threshisl', help='', type=float, default=3.0)
    parser.add_argument('-p', '--threshpix', help='', type=float, default=5.0)
    parser.add_argument('-r', '--rmsbox', help='rms box width and step (e.g., "(60, 20)")',
        type=str, default='(60, 20)')
    parser.add_argument('--rmsbox_bright', help='rms box for bright sources(?) width and step (e.g., "(60, 20)")',
        type=str, default='(60, 20)')
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)

    args = parser.parse_args()
    main(args.input_image, args.input_skymodel, args.output_skymodel,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
         rmsbox_bright=args.rmsbox_bright, adaptive_rmsbox=args.adaptive_rmsbox)
