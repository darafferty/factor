"""
Module that holds all screen-related functions
"""
import os
import numpy as np
import logging
import sys
import losoto.operations as operations
from losoto.h5parm import h5parm
import scipy.interpolate


log = logging.getLogger('factor:screens')


def generate_screens(fastphase_h5parm, slowgain_h5parm, bands):
    """
    Generate screens from input solutions
    """
    H_fast = h5parm(fastphase_h5parm, readonly=False)
    solset_fast = H_fast.getSolset('sol000')
    H_slow = h5parm(slowgain_h5parm, readonly=False)
    solset_slow = H_slow.getSolset('sol000')

    # calculate total phases at two frequencies
    tecsoltab = solset_fast.getSoltab('tec000')
    scphsoltab = solset_fast.getSoltab('scalarphase000')
    tec_vals = np.array(tecsoltab.val)
    scphase_vals = np.array(scphsoltab.val)
    source_names = tecsoltab.dir[:]
    times = tecsoltab.time[:]
    station_names = tecsoltab.ant[:]
    ref_anntennaid = 0
    freq = 140e6
    r = (-8.4479745e9 * tec_vals/freq) + scphase_vals # order is [time, ant, dir, freq]
    for idx in range(len(station_names)):
        r[:, idx, :, :] = r[:, idx, :, :] - r[:, ref_anntennaid, :, :]
    w = np.ones(r.shape)
    solset_fast.makeSoltab('phase', 'totalphase140MHz',
            axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times,
            station_names, source_names, np.array([freq])], vals=r, weights=w)
    freq = 160e6
    r = (-8.4479745e9 * tec_vals/freq) + scphase_vals # order is [time, ant, dir]
    for idx in range(len(station_names)):
        r[:, idx, :, :] = r[:, idx, :, :] - r[:, ref_anntennaid, :, :]
    w = np.ones(r.shape)
    solset_fast.makeSoltab('phase', 'totalphase160MHz',
            axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times,
            station_names, source_names, np.array([freq])], vals=r, weights=w)

    # estimate weights, first from past phases to get relative weights direction-to-direction,
    # then from slow phases to get time evolution of weights
    soltab = solset_slow.getSoltab('phase000')
    operations.reweight.run(soltab, 'window', nmedian=7, nstddev=11, flagBad=True)
    slow_times = soltab.time[:]
    slow_weights_phase = np.average(soltab.weight[:], axis=1) # average over frequency
    slow_weights_phase = np.average(slow_weights_phase, axis=-1) # average over pol

    soltab_freq1 = solset_fast.getSoltab('totalphase140MHz')
    operations.reweight.run(soltab_freq1, 'window', nmedian=3, nstddev=251, flagBad=True)
    fast_times = soltab_freq1.time[:]
    w = scipy.interpolate.interp1d(slow_times, slow_weights_phase, kind='linear',
        axis=0, fill_value='extrapolate')(fast_times)
    wr = np.reshape(w, soltab_freq1.weight[:].shape)
    wr = wr * np.median(soltab_freq1.weight[:], axis=0)/np.median(wr, axis=0)
    soltab_freq1.setValues(wr, weight=True)

    soltab_freq2 = solset_fast.getSoltab('totalphase160MHz')
    operations.reweight.run(soltab_freq2, 'window', nmedian=3, nstddev=251, flagBad=True)
    fast_times = soltab_freq2.time[:]
    wr = np.reshape(w, soltab_freq2.weight[:].shape)
    wr = wr * np.median(soltab_freq2.weight[:], axis=0)/np.median(wr, axis=0)
    soltab_freq2.setValues(wr, weight=True)

    # fit screens
    operations.phasescreen.run(soltab_freq1, 'phasescreen140MHz', niter=3, nsigma=3)
    operations.phasescreen.run(soltab_freq2, 'phasescreen160MHz', niter=3, nsigma=3)

    # plot screens
    if plot_screens:
        operations.plotscreen.run(soltab_freq1, prefix='140MHzphase')
        operations.plotscreen.run(soltab_freq2, prefix='160MHzphase')

    # Calculate phases at one or more frequencies
    soltab1 = solset.getSoltab('phasescreen140MHz')
    soltab2 = solset.getSoltab('phasescreen160MHz')
    source_dict = solset_fast.getSou()
    frequencies = []
    for band in bands:
        frequencies.extend(band.chan_freqs_hz)
    operations.screenphases.run(soltab_freq1, soltab_freq2, source_dict,
        frequencies, 'outputphases')
