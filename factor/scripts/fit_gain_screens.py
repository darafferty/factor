#! /usr/bin/env python
"""
Script to fit gain solutions with screens
"""
import argparse
from argparse import RawTextHelpFormatter
import losoto.operations as operations
from losoto.h5parm import h5parm
import numpy as np
import itertools
import multiprocessing
import math
import shutil
import multiprocessing
from scipy.interpolate import LSQUnivariateSpline, interp1d, interp2d
import sys
import scipy.ndimage
import astropy.convolution
from factor.lib import miscellaneous as misc


def std(inputData, Zero=False, axis=None, dtype=None):
    """
    Robust estimator of the standard deviation of a data set.

    Based on the robust_sigma function from the AstroIDL User's Library.

    .. versionchanged:: 1.0.3
        Added the 'axis' and 'dtype' keywords to make this function more
        compatible with np.std()
    """
    epsilon = 1.0e-20
    if axis is not None:
        fnc = lambda x: std(x, dtype=dtype)
        sigma = np.apply_along_axis(fnc, axis, inputData)
    else:
        data = inputData.ravel()
        if type(data).__name__ == "MaskedArray":
            data = data.compressed()
        if dtype is not None:
            data = data.astype(dtype)

        if Zero:
            data0 = 0.0
        else:
            data0 = np.median(data)
        maxAbsDev = np.median(np.abs(data-data0)) / 0.6745
        if maxAbsDev < epsilon:
            maxAbsDev = (np.abs(data-data0)).mean() / 0.8000
        if maxAbsDev < epsilon:
            sigma = 0.0
            return sigma

        u = (data-data0) / 6.0 / maxAbsDev
        u2 = u**2.0
        good = np.where( u2 <= 1.0 )
        good = good[0]
        if len(good) < 3:
            print("WARNING:  Distribution is too strange to compute standard deviation")
            sigma = -1.0
            return sigma

        numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
        nElements = (data.ravel()).shape[0]
        denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
        sigma = nElements*numerator / (denominator*(denominator-1.0))
        if sigma > 0:
            sigma = math.sqrt(sigma)
        else:
            sigma = 0.0

    return sigma

def findscatter(datavector):
    shifted_vec = np.roll(datavector, 1)
    scatter = np.nanmedian(abs(shifted_vec - datavector))
    return scatter


def findscatter_time(dataarray):
    scattervec = []
    for freq in range(0,len(dataarray[:,0])):
      scatter = findscatter(dataarray[freq,:])
      scattervec.append(scatter)
    return np.nanmedian(scattervec)


def findscatter_freq(dataarray):
    scattervec = []
    for time in range(0,len(dataarray[0,:])):
      scatter = findscatter(dataarray[:,time])
      scattervec.append(scatter)
    return np.nanmedian(scattervec)


def findnoisevec(datavector):
    shifted_vec = np.roll(datavector, 1)
    scatter_vec = (abs(shifted_vec - datavector))
    scatter_vec = scipy.ndimage.filters.median_filter(scatter_vec,9, mode='mirror')

    # now smooth
    gauss = astropy.convolution.Gaussian1DKernel(stddev=4.0)
    scatter_vec = astropy.convolution.convolve(scatter_vec,gauss , boundary='extend')

    # normalize scatter_vec
    scatter_vec = scatter_vec/np.mean(scatter_vec)

    return scatter_vec


def spline1D(amp_orig):
    # to compute knot points
    f = lambda m, n: [i*n//m + n//(2*m) for i in range(m)]

    if amp_orig is None:
        return None, None, None, None, None, None, None

    # expand array and mirror full array around edges
    ndata = len(amp_orig)
    amp = np.zeros(ndata+2*ndata)
    amp[ndata:ndata+ndata] = amp_orig

    for i in range(0, ndata):
        # Mirror at left edge.
        idx = min(ndata-1, ndata-i)
        amp[i] = amp_orig[idx]
        # Mirror at right edge
        idx = max(0, ndata-2-i)
        amp[ndata+ndata+i] = amp_orig[idx]

    # Find flagged values
    flagged = np.where(amp == 1.0)

    # work in log-sapce
    amp_orig_ext = np.copy(amp)
    amp = np.log10(amp)
    weights = (0.*np.copy(amp)) + 1 # initialize weights to 1

    # filter bad data and determine average scatter of amplitudes
    scatter = findscatter(amp)
    # remove some really bad stuff, by putting weights to zero.
    idxbadi1 = np.where(amp > (np.median(amp) + (35.*std(amp))))
    weights[idxbadi1] = 1e-10 # small value, zero generates NaN in spline
    idxbadi2 = np.where(amp < (np.median(amp) - (35.*std(amp))))
    weights[idxbadi2] = 1e-10  # small value, zero generates NaN in spline

    # Set weights for flagged values
    weights[flagged] = 1e-10  # small value, zero generates NaN in spline

    # make the noisevec
    if len(amp) > 30:  # so at least 30/3 = 10 good data points
        # create noise vector
        noisevec = findnoisevec(amp)
    else:
        noisevec = (np.copy(amp) * 0.) + 1.0 # just make constant noise, if we have too little datapoints

    if scatter < 0.005:
        #Interior knots t must satisfy Schoenberg-Whitney conditions
        scatter = 0.005 # otherwise we fit more parameters than we have data points
    knotfactor = 0.5e3*scatter  # normalize based on trial and error
    timevec = np.arange(0,len(amp))
    knotvec = f(np.int(len(amp)/knotfactor),len(amp))

    # simple optimization knot selection for vectors that have at least 30 data points
    # based on the noisevector
    # removes even numbered knots if the noise is high
    knotvec_copy = np.copy(knotvec) # otherwise tcopy is updated as well
    if len(timevec) > 30 and len(knotvec) > 2:
        for counter, knot in enumerate(knotvec_copy):
            if (counter % 2 == 0) and noisevec[knot] > 1.5: # even index and large noise
                knotvec.remove(knot)

    # asign midpoint if not enough data points/20
    if len (knotvec) < 3: # because we are working with a 3x larger mirrored array
        knotvec = [np.int(len(timevec)*0.25),np.int(len(timevec)/2),np.int(len(timevec)*0.75)]

    splineorder =  5 #  default
    if len(knotvec) == 3 and scatter > 0.1:
        splineorder = 3 # reduce order, data is  bad
        if scatter > 0.2:
            splineorder = 1 # very bad data
    spl2 = LSQUnivariateSpline(timevec, amp, knotvec, w=weights, k=splineorder)

    # now find bad data devatiating from the fit 15 x scatter
    residual = np.abs(spl2(timevec)-amp)
    idx      = np.where(residual > 15.*scatter)

    # second iteration
    if np.any(idx):
        ampcopy = np.copy(amp)
        ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
        spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

    residual = np.abs(spl2(timevec)-amp)
    idx      = np.where(residual > 8.*scatter)

    # third iteration
    if np.any(idx):
        ampcopy = np.copy(amp)
        ampcopy[idx] = spl2(timevec[idx]) # replace bad amplitudes by model
        spl2 = LSQUnivariateSpline(timevec, ampcopy, knotvec,  w=weights, k=splineorder)

    # again look at residual, go back to original amp again, find deviating data > 3x scatter
    residual = np.abs(spl2(timevec)-amp)
    idx      = np.where(residual > 3.*scatter)
    # replace the bad data with model
    model    =spl2(timevec)
    amp[idx] = model[idx]

    # go out of log-space
    idxnodata = np.where(np.logical_or(amp > 1.0, amp < -10.0))
    amp[idxnodata] = 0.0
    amp[flagged] = 0.0
    model[flagged] = 0.0
    amp = 10**amp

    amp_clean = amp[ndata:ndata + ndata]

    idxbad = np.where(amp_clean != amp_orig)
    n_knots = np.int(np.ceil(np.float(len(knotvec))/3.)) # approxmiate, just for plot

    # return cleaned amplitudes, model, scatter, number of knots, indices of replaced outliers
    return amp_clean, 10**(model[ndata:ndata + ndata]), noisevec[ndata:ndata + ndata], scatter, n_knots, idxbad, weights[ndata:ndata + ndata]


def pad_2Darray(a, width, mode):
    pad_shape = (a.shape[0]*3, a.shape[1]*3)
    pad_a = np.zeros(pad_shape)

    # center
    pad_a[a.shape[0]:2*a.shape[0], a.shape[1]:2*a.shape[1]] = a

    # four corners
    pad_a[0:a.shape[0], 0:a.shape[1]] = a[::-1, ::-1]
    pad_a[0:a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, ::-1]
    pad_a[2*a.shape[0]:3*a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, ::-1]
    pad_a[2*a.shape[0]:3*a.shape[0], 0:a.shape[1]] = a[::-1, ::-1]

    # middle edges
    pad_a[0:a.shape[0], a.shape[1]:2*a.shape[1]] = a[:, ::-1]
    pad_a[a.shape[0]:2*a.shape[0], 2*a.shape[1]:3*a.shape[1]] = a[::-1, :]
    pad_a[2*a.shape[0]:3*a.shape[0], a.shape[1]:2*a.shape[1]] = a[:, ::-1]
    pad_a[a.shape[0]:2*a.shape[0], 0:a.shape[1]] = a[::-1, :]

    return pad_a


def median2Dampfilter(amp_orig):
    try:
        from np import pad
    except ImportError:
        pad = pad_2Darray

    # pad array by reflection around axis
    orinal_size = np.shape(amp_orig)
    amp = pad(amp_orig, ((np.shape(amp_orig)[0],np.shape(amp_orig)[0]),
        (np.shape(amp_orig)[1],np.shape(amp_orig)[1])), mode='reflect')
    flagged = np.where(np.logical_or(amp == 1.0, amp == 0.0))

    # take the log
    amp = np.log10(amp)

    # Set flagged values to NaN
    amp[flagged] = np.nan

    # create median filtered array, ignoring NaNs
    amp_median = scipy.ndimage.filters.generic_filter(amp, np.nanmedian, (3,5)) # so a bit more smoothing along the time-axis

    # find scatter
    scatter_freq = findscatter_freq(amp)
    scatter_time = findscatter_time(amp)
    scatter = 0.5*(scatter_freq+scatter_time) # average x-y scatter

    # find bad data
    idxbad = np.where((np.abs(amp - amp_median)) > scatter*3.)
    baddata = np.copy(amp)*0.0
    baddata[idxbad] = 1.0

    # replace the bad data points
    amp_cleaned = np.copy(amp)
    amp_cleaned[idxbad] = amp_median[idxbad]

    # raise to the power
    amp = 10**amp
    amp_median = 10**amp_median
    amp_cleaned = 10**amp_cleaned

    #back to original size
    amp_median = amp_median[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]
    baddata   = baddata[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]
    amp_cleaned = amp_cleaned[orinal_size[0]:2*orinal_size[0],orinal_size[1]:2*orinal_size[1]]

    return amp_cleaned, amp_median, baddata


def smooth(soltab, smooth_amplitudes=True, normalize=True):

    if soltab.getType() == 'amplitude':
        # Check for flagged solutions. If found, set to 1
        parms = soltab.val[:]  # ['time', 'freq', 'ant', 'dir', 'pol']
        weights = soltab.weight[:]
        ntimes = len(soltab.time[:])
        nfreqs = len(soltab.freq[:])
        initial_flagged_indx = np.logical_or(np.isnan(parms), weights == 0.0)
        initial_unflagged_indx = np.logical_and(~np.isnan(parms), weights != 0.0)
        parms[initial_flagged_indx] = 1.0

        # Do smoothing
        if smooth_amplitudes:
            for dir in range(len(soltab.dir[:])):
                for pol in range(len(soltab.pol[:])):
                    for ant in range(len(soltab.ant[:])):
                        channel_amp_orig = parms[:, :, ant, dir, pol]

                        # Set up array for pool
                        channel_amp_interp = []
                        for chan in range(nfreqs):
                            unflagged_times = np.where(channel_amp_orig[:, chan] != 1.0)
                            flagged_times = np.where(channel_amp_orig[:, chan] == 1.0)
                            if np.any(flagged_times):
                                channel_amp_orig[flagged_times, chan] = 1.0
                            if np.any(unflagged_times):
                                channel_amp_interp.append(channel_amp_orig[:, chan])
                            else:
                                channel_amp_interp.append(None)

                        # Do 1-D smooth
                        if ntimes > 5:
                            pool = multiprocessing.Pool()
                            results = pool.map(spline1D, channel_amp_interp)
                            pool.close()
                            pool.join()

                            for chan, (amp_cleaned, model, noisevec, scatter, n_knots, idxbad, w) in enumerate(results):
                                # put back the results
                                if amp_cleaned is None:
                                    amp_cleaned = channel_amp_orig[:, chan]
                                parms[:, chan, ant, dir, pol] = np.copy(amp_cleaned)

                        # Do 2-D smooth
                        if nfreqs > 5:
                            channel_amp_orig = parms[:, :, ant, dir, pol]

                            # Smooth
                            channel_amp_interp = []
                            unflagged_sols = np.where(channel_amp_orig != 1.0)
                            if np.any(unflagged_sols):
                                # Set flagged solutions to 1.0
                                flagged_sols = np.where(np.logical_or(channel_amp_orig == 1.0, channel_amp_orig <= 0.0))
                                channel_amp_orig[flagged_sols] = 1.0

                                # Filter
                                amp_cleaned, amp_median, baddata = median2Dampfilter(channel_amp_orig.transpose([1, 0]))
                                amp_cleaned = amp_cleaned.transpose([1, 0])
                                amp_cleaned[flagged_sols] = 1.0
                                parms[:, :, ant, dir, pol] = np.copy(amp_cleaned)

        # Normalize the amplitude solutions to a mean of one across all channels
        if normalize:
            for dir in range(len(soltab.dir[:])):
                # First find the normalization factor from unflagged solutions
                norm_factor = 1.0/(np.nanmean(parms[:, :, :, dir, :][initial_unflagged_indx[:, :, :, dir, :]]))
                print("smooth_amps_spline.py: Normalization factor for direction {0} is {1}".format(dir, norm_factor))
                parms[:, :, :, dir, :] *= norm_factor

                # Clip extremely low amplitude solutions to prevent very high
                # amplitudes in the corrected data
                unflagged = np.where(~np.isnan(parms[:, :, :, dir, :]))
                low_ind = np.where(parms[:, :, :, dir, :][unflagged] < 0.2)
                parms[:, :, :, dir, :][unflagged][low_ind] = 0.2

        # Make sure flagged solutions are still flagged
        parms[initial_flagged_indx] = np.nan
        weights[initial_flagged_indx] = 0.0

    elif soltab.getType() == 'phase':
        # Flag solutions with large deviations from the mean
        nsigma = [2.0, 2.0, 3.0, 3.0, 4.0]
        parms = soltab.val[:]
        weights = soltab.weight[:]
        for dir in range(len(soltab.dir[:])):
            for pol in range(len(soltab.pol[:])):
                for ant in range(len(soltab.ant[:])):
                    phase_orig = parms[:, :, ant, dir, pol]
                    outliers = []
                    for nsig in nsigma:
                        ind = np.where(np.abs(phase_orig-np.mean(phase_orig, axis=0)) >
                                       nsig*np.std(phase_orig, axis=0))
                        if ind[0].size == 0:
                            break
                        else:
                            outliers.append(ind)
                            phase_orig[ind] = np.mean(phase_orig)
                    for ind in outliers:
                        parms[:, :, ant, dir, pol][ind] = np.NaN
                        weights[:, :, ant, dir, pol][ind] = 0.0

    else:
        print('Solution type must be phase or amplitude')
        sys.exit(1)

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


def main(h5parmfile, solsetname='sol000', ampsoltabname='amplitude000',
         phsoltabname='phase000', outsoltabroot='_screensols', ref_id=0,
         fit_screens=False, calculate_weights=False, smooth_amplitudes=False,
         smooth_phases=False, normalize=False):
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
    outsoltabroot : str, optional
        Root name for output soltabs
    ref_id : int, optional
        Index of reference station
    """
    ref_id = int(ref_id)
    normalize = misc.string2bool(normalize)
    fit_screens = misc.string2bool(fit_screens)
    calculate_weights = misc.string2bool(calculate_weights)
    smooth_amplitudes = misc.string2bool(smooth_amplitudes)
    smooth_phases = misc.string2bool(smooth_phases)

    # Read in solutions
    H = h5parm(h5parmfile, readonly=False)
    solset = H.getSolset(solsetname)
    ampsoltab = solset.getSoltab(ampsoltabname)
    amp = np.array(ampsoltab.val)
    damp = np.ones(amp.shape)
    phsoltab = solset.getSoltab(phsoltabname)
    ph = np.array(phsoltab.val)
    dph = np.ones(ph.shape)

    ampsoltab.rename('origamplitude000', overwrite=True)
    if smooth_amplitudes or normalize:
        amp, damp = smooth(ampsoltab, smooth_amplitudes=smooth_amplitudes, normalize=normalize)
    solset.makeSoltab('amplitude', 'amplitude000', axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                      axesVals=[ampsoltab.time[:], ampsoltab.freq[:], ampsoltab.ant[:],
                      ampsoltab.dir[:], ampsoltab.pol[:]], vals=amp, weights=damp)

    phsoltab.rename('origphase000', overwrite=True)
    if smooth_phases:
        ph, dph = smooth(phsoltab)
    solset.makeSoltab('phase', 'phase000', axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                      axesVals=[phsoltab.time[:], phsoltab.freq[:], phsoltab.ant[:],
                      phsoltab.dir[:], phsoltab.pol[:]], vals=ph, weights=dph)

    if fit_screens:
        if smooth_amplitudes:
            ampsoltab = solset.getSoltab('amplitude000')
        if ampsoltab:
            phsoltab = solset.getSoltab('phase000')

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
