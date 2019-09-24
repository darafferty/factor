#! /usr/bin/env python
"""
Script to filter a sky model with an image
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
import numpy as np
import os
import bdsf
from factor.lib import miscellaneous as misc
import casacore.tables as pt


def main(input_image, input_skymodel_nonpb, input_skymodel_pb, output_root,
         threshisl=3.0, threshpix=5.0, rmsbox=(150, 50), rmsbox_bright=(35, 7),
         adaptive_rmsbox=True, use_adaptive_threshold=False, adaptive_thresh=150.0,
         beamMS=None):
    """
    Filter the input sky model so that they lie in islands in the image

    Parameters
    ----------
    input_image : str
        Filename of input image to blank
    input_skymodel_nonpb : str
        Filename of input makesourcedb sky model, without primary-beam correction
    input_skymodel_pb : str, optional
        Filename of input makesourcedb sky model, with primary-beam correction
    output_root : str
        Root of filename of output makesourcedb sky models. Output filenames will be
        output_root+'.apparent_sky' and output_root+'.true_sky'
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
    if rmsbox is not None and isinstance(rmsbox, str):
        rmsbox = eval(rmsbox)
    if isinstance(rmsbox_bright, str):
        rmsbox_bright = eval(rmsbox_bright)
    adaptive_rmsbox = misc.string2bool(adaptive_rmsbox)
    use_adaptive_threshold = misc.string2bool(use_adaptive_threshold)
    if isinstance(beamMS, str):
        beamMS = misc.string2list(beamMS)

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

        s = lsmtool.load(input_skymodel_nonpb)
        s.select('{} == True'.format(maskfile))  # keep only those in PyBDSF masked regions
        s.group(maskfile)  # group the sky model by mask islands
        s.write(output_root+'.apparent_sky', clobber=True)

        if not os.path.exists(input_skymodel_pb):
            # No true-sky model available, so use apparent-sky one and correct for beam
            if len(beamMS) > 1:
                ms_times = []
                for ms in beamMS:
                    tab = pt.table(ms, ack=False)
                    ms_times.append(np.mean(tab.getcol('TIME')))
                    tab.close()
                ms_times_sorted = sorted(ms_times)
                mid_time = ms_times_sorted[int(len(ms_times)/2)]
                beam_ind = ms_times.index(mid_time)
            else:
                beam_ind = 0
            s = lsmtool.load(input_skymodel_nonpb, beamMS=beamMS[beam_ind])
            applyBeam = True
            invertBeam = True
            adjustSI = True
        else:
            s = lsmtool.load(input_skymodel_pb)
            applyBeam = False
            invertBeam = False
            adjustSI = False
        s.select('{} == True'.format(maskfile))  # keep only those in PyBDSF masked regions
        s.group(maskfile)  # group the sky model by mask islands
        s.write(output_root+'.true_sky', clobber=True, applyBeam=applyBeam,
                invertBeam=invertBeam, adjustSI=adjustSI)


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
    parser.add_argument('--rmsbox_bright', help='rms box for bright sources, width and step (e.g., "(60, 20)")',
                        type=str, default='(60, 20)')
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)

    args = parser.parse_args()
    main(args.input_image, args.input_skymodel, args.output_skymodel,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
         rmsbox_bright=args.rmsbox_bright, adaptive_rmsbox=args.adaptive_rmsbox)
