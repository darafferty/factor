#! /usr/bin/env python
"""
Script to convert selfcal solutions to single gain table
"""
import argparse
from argparse import RawTextHelpFormatter
from losoto.h5parm import h5parm
from scipy.interpolate import interp1d
import shutil
import numpy as np
import sys
import os


def calc_total(tec_phase, amp_val, phase_val, fast_times, slow_times, slow_freqs, final_freqs):
    """
    Calculates total phases and amplitudes

    Parameters
    ----------
    tec_phase : array
        Fast phase values from TEC
    amp_val : array
        Slow amplitude values
    phase_val : array
        Slow phase values
    fast_times : array
        Fast phase times
    slow_times : array
        Slow amp/phase times
    final_freqs : array
        Final frequencies
    """
    # Interpolate in time and frequency
    f_amp = interp1d(slow_times, amp_val, kind='nearest',
                     axis=0, bounds_error=False, fill_value='extrapolate')
    vals = np.squeeze(f_amp(fast_times))
    f_amp = interp1d(slow_freqs, vals, kind='nearest',
                     axis=1, bounds_error=False, fill_value='extrapolate')
    total_amp = np.squeeze(f_amp(final_freqs))

    f_phase = interp1d(slow_times, phase_val, kind='nearest',
                       axis=0, bounds_error=False, fill_value='extrapolate')
    vals = np.squeeze(f_phase(fast_times))
    f_phase = interp1d(slow_freqs, vals, kind='nearest',
                       axis=1, bounds_error=False, fill_value='extrapolate')
    total_phase = np.mod((np.squeeze(f_phase(final_freqs)) + tec_phase + np.pi), 2*np.pi) - np.pi

    return total_amp, total_phase


def main(fast_h5parm, slow_h5parm, output_file, freqstep=1, scratch_dir=None):
    """
    Converts multiple selfcal tables to single gain table

    Parameters
    ----------
    fast_h5parm : str
        File with fast phase (TEC and CommonScalarPhase) solutions
    fast_5parm : str
        File with slow gain solutions
    output_file : str
        Output filename
    freqstep : int
        Frequency step to divide up frequency width of solutions
    scratch_dir : str, optional
        Scratch directory for temp storage

    """
    freqstep = int(freqstep)

    # Copy to scratch directory if specified
    if scratch_dir is not None:
        fast_h5parm_orig = fast_h5parm
        fast_h5parm = os.path.join(scratch_dir, os.path.basename(fast_h5parm_orig))
        slow_h5parm_orig = slow_h5parm
        slow_h5parm = os.path.join(scratch_dir, os.path.basename(slow_h5parm_orig))
        output_file_orig = output_file
        output_file = os.path.join(scratch_dir, os.path.basename(output_file_orig))
        shutil.copytree(fast_h5parm_orig, fast_h5parm)
        shutil.copytree(slow_h5parm_orig, slow_h5parm)

    fast_h5 = h5parm(fast_h5parm)
    fast_solset = fast_h5.getSolset('sol000')
    tec_soltab = fast_solset.getSoltab('tec000', useCache=True)
    slow_h5 = h5parm(slow_h5parm)
    slow_solset = slow_h5.getSolset('sol000')
    amp_soltab = slow_solset.getSoltab('amplitude000', useCache=True)
    phase_soltab = slow_solset.getSoltab('phase000', useCache=True)

    if os.path.exists(output_file):
        os.remove(output_file)
    output_h5 = h5parm(output_file, readonly=False)
    output_solset = output_h5.makeSolset(solsetName = 'sol000', addTables=False)
    fast_solset.obj._f_copy_children(output_solset.obj, recursive=True, overwrite=True)
    slow_solset.obj._f_copy_children(output_solset.obj, recursive=True, overwrite=True)

    # Get various quantities over which we must iterate
    station_names = tec_soltab.ant[:]
    fast_times = tec_soltab.time[:]
    slow_times = amp_soltab.time[:]
    slow_freqs = amp_soltab.freq[:]

    # Determine final frequency grid
    freqwidth = slow_freqs[1] - slow_freqs[0]
    final_freqwidth = freqwidth / freqstep
    if freqstep > 1:
        final_freqs = []
        for freq in slow_freqs:
            low_freq = freq - freqwidth / 2
            for i in range(freqstep):
                final_freqs.append(low_freq + final_freqwidth * (i + 0.5))
        final_freqs = np.array(final_freqs)
    else:
        final_freqs = slow_freqs

    # Make output soltabs
    ntim = len(fast_times)
    nfre = len(final_freqs)
    ants = amp_soltab.ant[:]
    nant = len(ants)
    dirs = amp_soltab.dir[:]
    ndir = len(dirs)
    pols = amp_soltab.pol[:]
    npol = len(pols)
    output_amp_soltab = output_solset.makeSoltab('amplitude', 'totalamplitude',
                            axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                            axesVals=[fast_times, final_freqs, amp_soltab.ant[:],
                                      amp_soltab.dir[:], amp_soltab.pol[:]],
                            vals=np.zeros((ntim, nfre, nant, ndir, npol)),
                            weights=np.ones((ntim, nfre, nant, ndir, npol)))
    output_amp_soltab.setCache(np.zeros((ntim, nfre, nant, ndir, npol)),
                               np.ones((ntim, nfre, nant, ndir, npol)))
    output_amp_soltab.useCache = True
    output_phase_soltab = output_solset.makeSoltab('phase', 'totalphase',
                            axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                            axesVals=[fast_times, final_freqs, amp_soltab.ant[:],
                                      amp_soltab.dir[:], amp_soltab.pol[:]],
                            vals=np.zeros((ntim, nfre, nant, ndir, npol)),
                            weights=np.ones((ntim, nfre, nant, ndir, npol)))
    output_phase_soltab.setCache(np.zeros((ntim, nfre, nant, ndir, npol)),
                                 np.ones((ntim, nfre, nant, ndir, npol)))
    output_phase_soltab.useCache = True

    # Add the phase and amp corrections together
    for dir in dirs:
        for pol in pols:
            tec_phase_list = []
            amp_list = []
            phase_list = []
            for station in ants:
                tec_soltab.setSelection(dir=[dir], ant=[station])
                tec_phase = (-8.44797245e9 * np.column_stack([np.squeeze(tec_soltab.val[:])]*nfre) /
                             np.resize(final_freqs, (ntim, nfre)))
                amp_soltab.setSelection(dir=[dir], ant=[station], pol=[pol])
                phase_soltab.setSelection(dir=[dir], ant=[station], pol=[pol])
                time_ind = amp_soltab.axesNames.index('time')

                total_amp, total_phase = calc_total(tec_phase, amp_soltab.val,
                                                    phase_soltab.val, fast_times,
                                                    slow_times, slow_freqs, final_freqs)
                output_amp_soltab.setSelection(dir=[dir], ant=[station], pol=[pol])
                output_amp_soltab.setValues(total_amp)
                output_phase_soltab.setSelection(dir=[dir], ant=[station], pol=[pol])
                output_phase_soltab.setValues(total_phase)
    output_amp_soltab.flush()
    output_phase_soltab.flush()
    output_h5.close()
    fast_h5.close()
    slow_h5.close()

    # Copy output to original path and delete copies if scratch directory is specified
    if scratch_dir is not None:
        if os.path.exists(output_file_orig):
            shutil.rmtree(output_file_orig)
        shutil.copytree(output_file, output_file_orig)
        shutil.rmtree(output_file)
        shutil.rmtree(fast_h5parm)
        shutil.rmtree(slow_h5parm)


if __name__ == '__main__':
    descriptiontext = "Converts multiple selfcal tables to single gain table.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fast_selfcal_h5parm', help='name of the h5parm with fast solutions')
    parser.add_argument('slow_selfcal_h5parm', help='name of the h5parm with slow solutions')
    parser.add_argument('output_file', help='name of the output file')
    args = parser.parse_args()

    main(args.fast_selfcal_h5parm, args.slow_selfcal_h5parm, args.output_file)
