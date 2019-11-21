#! /usr/bin/env python
"""
Script to fit gain solutions with screens
"""
import losoto.operations as operations
from losoto.h5parm import h5parm
import numpy as np
import scipy.interpolate as si
from factor.lib import miscellaneous as misc
from loess import loess_2d, loess_1d


def normalize_values(soltab):
    """
    Normalize values so that mean is equal to unity for amplitudes and zero for phases

    Parameters
    ----------
    soltab : solution table
        Input table with solutions

    Returns
    -------
    parms, weights : arrays
        The normalized parameters and weights
    """
    soltype = soltab.getType()
    parms = soltab.val[:]  # ['time', 'freq', 'ant', 'dir', 'pol']
    weights = soltab.weight[:]
    initial_flagged_indx = np.logical_or(np.isnan(parms), weights == 0.0)
    initial_unflagged_indx = np.logical_and(~np.isnan(parms), weights != 0.0)
    parms[initial_flagged_indx] = np.nan

    # Normalize each station, time, and frequency separately to unity, but
    # maintain the relative offset between polarizations. Note that we work in log space
    # for these operations if the solutions are amplitudes
    if soltype == 'amplitude':
        parms = np.log(parms)
    for dir in range(len(soltab.dir[:])):
        for s in range(len(soltab.ant[:])):
            for t in range(len(soltab.time[:])):
                for f in range(len(soltab.freq[:])):
                    norm_factor = np.nanmean(parms[t, f, s, dir, :][initial_unflagged_indx[t, f, s, dir, :]])
                    parms[t, f, s, dir, :] -= norm_factor

    # Make sure flagged solutions are still flagged
    if soltype == 'amplitude':
        parms = 10**parms
    parms[initial_flagged_indx] = np.nan
    weights[initial_flagged_indx] = 0.0

    return parms, weights


def find_bandpass_correction(soltab, parms_normalized):
    """
    Find the correction to the bandpass for each station

    Parameters
    ----------
    soltab : solution table
        Input table with solutions
    parms_normalized : array
        The normalized parameters

    Returns
    -------
    parms, weights : arrays
        The parameters with bandpass corrections and weights
    """
    parms = soltab.val[:]  # ['time', 'freq', 'ant', 'dir', 'pol']
    parms /= parms_normalized  # divide out XX-YY offsets
    weights = soltab.weight[:]
    initial_flagged_indx = np.logical_or(np.isnan(parms), weights == 0.0)
    initial_unflagged_indx = np.logical_and(~np.isnan(parms), weights != 0.0)
    parms[initial_flagged_indx] = np.nan

    # Normalize each direction, time, and station to have a mean of 1 over all frequencies
    # and polarizations. Note that we work in log space for these operations
    parms = np.log(parms)
    for dir in range(len(soltab.dir[:])):
        for s in range(len(soltab.ant[:])):
            for t in range(len(soltab.time[:])):
                norm_factor = np.nanmean(parms[t, :, s, dir, :][initial_unflagged_indx[t, :, s, dir, :]])
                parms[t, :, s, dir, :] -= norm_factor

    # Take the median over all directions, times, and pols for each station at each freq
    # to reduce the noise on the bandpass correction, then smooth in frequency with LOESS.
    for s in range(len(soltab.ant[:])):
        for f in range(len(soltab.freq[:])):
            parms[:, f, s, :, :] = np.nanmedian(parms[:, f, s, :, :])
        nanind = np.where(~np.isnan(parms[0, :, s, 0, 0]))
        fs, ps, ws = loess_1d.loess_1d(soltab.freq[nanind]/100e8, parms[0, :, s, 0, 0][nanind], frac=0.5, degree=2)
        parms[0, :, s, 0, 0][nanind] = ps
    parms = 10**parms
    parms *= parms_normalized  # put XX-YY offsets back

    # Make sure flagged solutions are still flagged
    parms[initial_flagged_indx] = np.nan
    weights[initial_flagged_indx] = 0.0

    return parms, weights


def smooth(soltab, stddev_threshold=0.25, freq_sampling=5, time_sampling=2):
    """
    Smooth scalarphases. The smoothing is done in real and imaginary space

    Parameters
    ----------
    soltab : solution table
        Input table with solutions
    stddev_threshold : float, optional
        The threshold stddev below which no smoothing is done
    freq_sampling : int, optoinal
        Sampling stride to use for frequency when doing LOESS smooth
    time_sampling : int, optoinal
        Sampling stride to use for time when doing LOESS smooth

    Returns
    -------
    parms, weights : arrays
        The parameters with bandpass corrections and weights
    """
    parms = soltab.val[:]  # ['time', 'freq', 'ant', 'dir']
    weights = soltab.weight[:]

    # Find gaps in time and treat each block separately
    times = soltab.times[:]
    delta_times = times[1:] - times[:-1]  # time at center of solution interval
    timewidth = np.min(delta_times)
    gaps = np.where(delta_times > timewidth*1.2)
    gaps_ind = gaps[0] + 1
    gaps_ind = np.append(gaps_ind, np.array([len(times)]))

    for dir in range(len(soltab.dir[:])):
        # Find standard deviation of the real part of a core station and determine
        # whether we need to smooth this direction or not
        csindx = 2  # should always be a core station
        sdev = np.stddev(np.cos(parms[:, :, csindx, dir]))
        if sdev >= stddev_threshold:
            # Set frac depending on sdev
            frac = min(0.8, 0.3 * sdev / stddev_threshold)
            for s in range(len(soltab.ant[:])):
                g_start = 0
                for gnum, g_stop in enumerate(gaps_ind):
                    # Define slices for frequency and time sampling
                    freq_slice = slice(0, -1, freq_sampling)
                    time_slice = slice(g_start, g_stop, time_sampling)

                    # Do the smoothing with LOESS
                    yv, xv = np.meshgrid(soltab.freq[freq_slice], times[time_slice])
                    yv /= np.min(yv)
                    xv -= np.min(xv)
                    zreal = np.cos(parms[time_slice, freq_slice, s, dir])
                    zsreal, wreal = loess_2d.loess_2d(xv.flatten(), yv.flatten(), zreal.flatten(),
                                                      rescale=True, frac=frac, degree=2)
                    zimag = np.sin(parms[time_slice, freq_slice, s, dir])
                    zsimag, wimag = loess_2d.loess_2d(xv.flatten(), yv.flatten(), zimag.flatten(),
                                                      rescale=True, frac=frac, degree=2)

                    # Interpolate back to original grid if needed
                    if freq_sampling > 1 or time_sampling > 1:
                        zr = zsreal.reshape((len(times[time_slice]), len(soltab.freq[freq_slice])))
                        f = si.interp1d(times[time_slice], zr, axis=0, kind='linear', fill_value='extrapolate')
                        zr1 = f(times)
                        f = si.interp1d(soltab.freq[freq_slice], zr1, axis=1, kind='linear', fill_value='extrapolate')
                        zr = f(soltab.freq)
                        zi = zsimag.reshape((len(times[time_slice]), len(soltab.freq[freq_slice])))
                        f = si.interp1d(times[time_slice], zi, axis=0, kind='linear', fill_value='extrapolate')
                        zi1 = f(times)
                        f = si.interp1d(soltab.freq[freq_slice], zi1, axis=1, kind='linear', fill_value='extrapolate')
                        zi = f(soltab.freq)
                    else:
                        zr = zsreal.reshape((len(times), len(soltab.freq)))
                        zi = zsimag.reshape((len(times), len(soltab.freq)))
                    parms[:, :, s, dir] = np.arctan2(zr, zi)
                    g_start = g_stop

    return parms, weights


def remove_soltabs(solset, soltabnames):
    """
    Remove soltab
    """
    for soltabname in soltabnames:
        try:
            soltab = solset.getSoltab(soltabname)
            soltab.delete()
        except:
            pass


def main(h5parmfile, solsetname='sol000', ampsoltabname=None,
         phsoltabname=None, ref_id=0, fit_screens=False, calculate_weights=False,
         smooth_phases=False, normalize=False, find_bandpass=False):
    """
    Fit screens to gain solutions

    Parameters
    ----------
    h5parmfile : str
        Filename of h5parm
    solsetname : str, optional
        Name of solset
    ampsoltabname : str, optional
        Name of TEC soltab
    phsoltabname : str, optional
        Name of error soltab
    ref_id : int, optional
        Index of reference station
    """
    ref_id = int(ref_id)
    normalize = misc.string2bool(normalize)
    fit_screens = misc.string2bool(fit_screens)
    calculate_weights = misc.string2bool(calculate_weights)
    smooth_phases = misc.string2bool(smooth_phases)
    find_bandpass = misc.string2bool(find_bandpass)

    # Read in solutions
    H = h5parm(h5parmfile, readonly=False)
    solset = H.getSolset(solsetname)
    if ampsoltabname is not None:
        ampsoltab = solset.getSoltab(ampsoltabname)
        amp = np.array(ampsoltab.val)
        damp = np.ones(amp.shape)

        if ampsoltabname != 'origamplitude000':
            ampsoltab.rename('origamplitude000', overwrite=True)
        if normalize:
            amp, damp = normalize_values(ampsoltab)
            if find_bandpass:
                amp, damp = find_bandpass_correction(ampsoltab, amp)
        solset.makeSoltab('amplitude', 'amplitude000', axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                          axesVals=[ampsoltab.time[:], ampsoltab.freq[:], ampsoltab.ant[:],
                          ampsoltab.dir[:], ampsoltab.pol[:]], vals=amp, weights=damp)

    if phsoltabname is not None:
        phsoltab = solset.getSoltab(phsoltabname)
        ph = np.array(phsoltab.val)
        dph = np.ones(ph.shape)
        axis_names = phsoltab.getAxesNames()

        # Identify any duplicate times and remove
        times = phsoltab.time[:]
        delta_times = times[1:] - times[:-1]
        nodupind = np.where(delta_times > 0.1)
        if len(nodupind[0]) < len(times):
            times = times[nodupind]
            if 'pol' in axis_names:
                ph = np.squeeze(ph[nodupind, :, :, :, :])
                dph = np.squeeze(dph[nodupind, :, :, :, :])
            else:
                ph = np.squeeze(ph[nodupind, :, :, :])
                dph = np.squeeze(dph[nodupind, :, :, :])
            if phsoltabname != 'origphase000':
                phsoltab.rename('origphase000', overwrite=True)
            if 'pol' in axis_names:
                solset.makeSoltab('phase', 'phase000', axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                                  axesVals=[times, phsoltab.freq[:], phsoltab.ant[:],
                                  phsoltab.dir[:], phsoltab.pol[:]], vals=ph, weights=dph)
            else:
                solset.makeSoltab('phase', 'phase000', axesNames=['time', 'freq', 'ant', 'dir'],
                                  axesVals=[times, phsoltab.freq[:], phsoltab.ant[:],
                                  phsoltab.dir[:]], vals=ph, weights=dph)
            phsoltab = solset.getSoltab('phase000')
            ph = np.array(phsoltab.val)
            dph = np.ones(ph.shape)

        if phsoltabname != 'origphase000':
            phsoltab.rename('origphase000', overwrite=True)
        if smooth_phases:
            ph, dph = smooth(phsoltab)
        if normalize:
            ph, dph = normalize_values(phsoltab)
        if 'pol' in axis_names:
            solset.makeSoltab('phase', 'phase000', axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                              axesVals=[times, phsoltab.freq[:], phsoltab.ant[:],
                              phsoltab.dir[:], phsoltab.pol[:]], vals=ph, weights=dph)
        else:
            solset.makeSoltab('phase', 'phase000', axesNames=['time', 'freq', 'ant', 'dir'],
                              axesVals=[times, phsoltab.freq[:], phsoltab.ant[:],
                              phsoltab.dir[:]], vals=ph, weights=dph)

    if fit_screens:
        # Find weights
        if calculate_weights:
            operations.reweight.run(ampsoltab, 'window', nmedian=3, nstddev=501)
            operations.reweight.run(phsoltab, 'window', nmedian=3, nstddev=501)

        # Fit screens
        remove_soltabs(solset, ['amplitudescreen000', 'amplitudescreen000resid'])
        operations.stationscreen.run(ampsoltab, 'amplitudescreen000', niter=1, nsigma=5,
            refAnt=ref_id, order=20, scale_order=False)
        remove_soltabs(solset, ['phasescreen000', 'phasescreen000resid'])
        operations.stationscreen.run(phsoltab, 'phasescreen000', niter=1, nsigma=5,
            refAnt=ref_id, order=20, scale_order=False)

    H.close()
