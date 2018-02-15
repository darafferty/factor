#! /usr/bin/env python
"""
Script to fit TEC solutions with screens
"""
import argparse
from argparse import RawTextHelpFormatter
import losoto.operations as operations
from losoto.h5parm import h5parm
import numpy as np
import matplotlib.pyplot as plt
from losoto.operations_lib import normalize_phase
from scipy import optimize
import pickle
import itertools
import multiprocessing


def convert_mvt(mvt):
    """
    Converts casacore MVTime to MJD

    Parameters
    ----------
    mvt : str
        MVTime

    Returns
    -------
    mjd : float
        MJD in seconds
    """
    t = Time(mjd_sec / 3600 / 24, format='mjd', scale='utc')
    date, hour = t.iso.split(' ')
    year, month, day = date.split('-')
    d = t.datetime
    month = d.ctime().split(' ')[1]

    return '{0}{1}{2}/{3}'.format(day, month, year, hour)


def _rolling_window_lastaxis(a, window):
    """Directly taken from Erik Rigtorp's post to numpy-discussion.
    <http://www.mail-archive.com/numpy-discussion@scipy.org/msg29450.html>"""
    import numpy as np

    if window < 1:
       raise ValueError, "`window` must be at least 1."
    if window > a.shape[-1]:
       raise ValueError, "`window` is too long."
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def most_common(lst):
    return max(set(lst), key=lst.count)


def remove_jumps_pool(inputs):
    """
    Simple helper function for pool.map
    """
    return remove_jumps(*inputs)


def remove_jumps(tec_vals, jump_val, nsamples):

    # Get random start positions
    startpoints = []
    nstddev = 21
    pad_width = [(0, 0)] * len(tec_vals.shape)
    pad_width[-1] = ((nstddev-1)/2, (nstddev-1)/2)
    pad_vals = np.pad(tec_vals, pad_width, 'constant', constant_values=(np.nan,))
    stddev = np.nanstd(_rolling_window_lastaxis(pad_vals, nstddev), axis=-1)
    noisy_ind = np.where(stddev > 3.0*np.mean(stddev))
    while len(startpoints) < nsamples:
        # Select only those start positions that lie in regions of low scatter
        # with values near zero
        startpoints = np.random.random_integers(0, len(tec_vals)-1, nsamples*10)
        startpoints = [s for s in startpoints if (s not in noisy_ind[0] and tec_vals[s] < jump_val/2.0)]
        if len(startpoints) > nsamples:
            startpoints = startpoints[:nsamples]

    # For each start position, move backward and forward in time, removing jumps
    min_err = jump_val / 5.0
    max_err = jump_val / 2.0
    tec_samples = []
    for s, startindx in enumerate(startpoints):
        tec = tec_vals.copy()
        if startindx > 0:
            # Search backwards
            skip = False
            for i in range(startindx, 0, -1):
                # Remove jumps
                if not skip:
                    diff = np.fmod(tec[i]-tec[i-1], jump_val)
                    tec[i-1] = tec[i] - diff
                    err = 0.025 #max(0.02, np.mean(stddev))
                    if approx_equal(abs(diff), jump_val, tol=err):
                        tec[i-1] += np.sign(diff) * jump_val
                    if i > 1:
                        # Check diff to next solution as well
                        diff = np.fmod(tec[i]-tec[i-2], jump_val)
                        tec[i-2] = tec[i] - diff
                        if approx_equal(abs(diff), jump_val, tol=err):
                            tec[i-2] += np.sign(diff) * jump_val
                            skip = False
                else:
                    skip = False

        if startindx < len(tec_vals)-1:
            # Search forwards
            skip = False
            for i in range(startindx, len(tec_vals)-1, 1):
                # Remove jumps
                if not skip:
                    diff = np.fmod(tec[i]-tec[i+1], jump_val)
                    tec[i+1] = tec[i] - diff
                    err = 0.025 #max(0.02, np.mean(stddev))
                    if approx_equal(abs(diff), jump_val, tol=err):
                        tec[i+1] += np.sign(diff) * jump_val
                    if i < len(tec_vals)-2:
                        # Check diff to next solution as well
                        diff = np.fmod(tec[i]-tec[i+2], jump_val)
                        tec[i+2] = tec[i] - diff
                        if approx_equal(abs(diff), jump_val, tol=err):
                            tec[i+2] += np.sign(diff) * jump_val
                            skip = False
                else:
                    skip = False
        tec_samples.append(tec)

    # Take the median as the best guess
    tec_samples = np.array(tec_samples)
    tec = np.median(np.array(tec_samples), axis=0)
#     tecind = np.argmin(np.abs(np.array(tec_samples)), axis=0)
#     for i in range(len(tec)):
#         tec[i] = tec_samples[tecind[i], i]
#         tec[i] = most_common(tec_samples[:, i].tolist())
    nmed = 5
    pad_width = [(0, 0)] * len(tec_vals.shape)
    pad_width[-1] = ((nmed-1)/2, (nmed-1)/2)
    pad_vals = np.pad(tec, pad_width, 'constant', constant_values=(np.nan,))
    medtec = np.nanmedian(_rolling_window_lastaxis(pad_vals, nmed), axis=-1)

    nstddev = 15
    pad_width = [(0, 0)] * len(tec_vals.shape)
    pad_width[-1] = ((nstddev-1)/2, (nstddev-1)/2)
    pad_vals = np.pad(tec-medtec, pad_width, 'constant', constant_values=(np.nan,))
    stddev = np.nanstd(_rolling_window_lastaxis(pad_vals, nstddev), axis=-1) * 3.0
    tec_samples = []
#     tec_vals = np.array(tec).copy()
    stop_err = 0.025
    for s, startindx in enumerate(startpoints):
        tec = tec_vals.copy()
        if startindx > 0:
            # Search backwards
            skip = False
            for i in range(startindx, 0, -1):
                # Remove jumps
                if not skip:
                    diff = np.fmod(tec[i]-tec[i-1], jump_val)
                    tec[i-1] = tec[i] - diff
                    err = min(max_err, max(0.025, stddev[i-1]))
                    if err > stop_err:
                        tec[0:i-1] = np.nan
                        break
                    if approx_equal(abs(diff), jump_val, tol=err):
                        tec[i-1] += np.sign(diff) * jump_val
                    if i > 1:
                        # Check diff to next solution as well
                        old_val = tec[i-2]
                        diff = np.fmod(tec[i]-tec[i-2], jump_val)
                        tec[i-2] = tec[i] - diff
                        if approx_equal(abs(diff), jump_val, tol=err):
                            tec[i-2] += np.sign(diff) * jump_val
                        if abs(old_val - tec[i-2]) > jump_val*0.9:
                            # If we removed a jump, set middle value to average
                            tec[i-1] = (tec[i] + tec[i-2]) / 2.0
                            skip = True
                else:
                    skip = False
        if startindx < len(tec_vals)-1:
            # Search forwards
            skip = False
            for i in range(startindx, len(tec_vals)-1, 1):
                # Remove jumps
                if not skip:
                    diff = np.fmod(tec[i]-tec[i+1], jump_val)
                    tec[i+1] = tec[i] - diff
                    err = min(max_err, max(0.015, stddev[i+1]))
                    if err > stop_err:
                        tec[i+1:] = np.nan
                        break
                    if approx_equal(abs(diff), jump_val, tol=err):
                        tec[i+1] += np.sign(diff) * jump_val
                    if i < len(tec_vals)-2:
                        # Check diff to next solution as well
                        old_val = tec[i+2]
                        diff = np.fmod(tec[i]-tec[i+2], jump_val)
                        tec[i+2] = tec[i] - diff
                        if approx_equal(abs(diff), jump_val, tol=err):
                            tec[i+2] += np.sign(diff) * jump_val
                            skip = False
                        if abs(old_val - tec[i+2]) > jump_val*0.9:
                            # If we removed a jump, set middle value to average
                            tec[i+1] = (tec[i] + tec[i+2]) / 2.0
                            skip = True
                else:
                    skip = False
        tec_samples.append(tec)

    # Take the median as the best guess
#     tec_samples = np.array(tec_samples)
#     0/0
    tec = np.nanmedian(np.array(tec_samples), axis=0)
#     tecind = np.argmin(np.abs(np.array(tec_samples)), axis=0)
#     for i in range(len(tec)):
#         tec[i] = tec_samples[tecind[i], i]

    return tec, stddev


def _float_approx_equal(x, y, tol=1e-18, rel=None):
    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors are None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)


def approx_equal(x, y, *args, **kwargs):
    """approx_equal(float1, float2[, tol=1e-18, rel=1e-7]) -> True|False
    approx_equal(obj1, obj2[, *args, **kwargs]) -> True|False

    Return True if x and y are approximately equal, otherwise False.

    If x and y are floats, return True if y is within either absolute error
    tol or relative error rel of x. You can disable either the absolute or
    relative check by passing None as tol or rel (but not both).

    For any other objects, x and y are checked in that order for a method
    __approx_equal__, and the result of that is returned as a bool. Any
    optional arguments are passed to the __approx_equal__ method.

    __approx_equal__ can return NotImplemented to signal that it doesn't know
    how to perform that specific comparison, in which case the other object is
    checked instead. If neither object have the method, or both defer by
    returning NotImplemented, approx_equal falls back on the same numeric
    comparison used for floats.

    >>> almost_equal(1.2345678, 1.2345677)
    True
    >>> almost_equal(1.234, 1.235)
    False

    """
    if not (type(x) is type(y) is float):
        # Skip checking for __approx_equal__ in the common case of two floats.
        methodname = '__approx_equal__'
        # Allow the objects to specify what they consider "approximately equal",
        # giving precedence to x. If either object has the appropriate method, we
        # pass on any optional arguments untouched.
        for a,b in ((x, y), (y, x)):
            try:
                method = getattr(a, methodname)
            except AttributeError:
                continue
            else:
                result = method(b, *args, **kwargs)
                if result is NotImplemented:
                    continue
                return bool(result)
    # If we get here without returning, then neither x nor y knows how to do an
    # approximate equal comparison (or are both floats). Fall back to a numeric
    # comparison.
    return _float_approx_equal(x, y, *args, **kwargs)



def calculate_total_tec(tec_vals, scphase_vals, freq, ant_axis_ind, ref_id=0,
    bootstrap=False, station_names=None, baseline_file=None):

    # Calculate total tec
    if scphase_vals is not None:
        r = tec_vals + normalize_phase(scphase_vals) / (-8.4479745e9) * freq
    else:
        r = tec_vals

    # Adjust all to reference station
    r_ref = r.copy()
    nstations = tec_vals.shape[ant_axis_ind]
    naxes = len(tec_vals.shape)
    ref_indlist = [ref_id if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
    for idx in range(nstations):
        indlist = [idx if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
        r_ref[indlist] = r[indlist] - r[ref_indlist]

    if bootstrap:
        # For distant stations, use the nearest as the reference
        distant_stations, nearest_stations = get_bootstrap_stations()
        for distant_station, nearest_station in zip(distant_stations, nearest_stations):
            idx = station_names.tolist().index(distant_station)
            indlist = [idx if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
            nearest_ind = station_names.tolist().index(nearest_station)
            ref_indlist = [nearest_ind if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
            r_ref[indlist] = r[indlist] - r[ref_indlist]

    return r_ref


def get_bootstrap_stations():
    """
    Returns lists of distant stations and their nearest neighbor
    """
    distant_stations = ['RS509HBA', 'RS508HBA', 'RS409HBA', 'RS407HBA', 'RS310HBA', 'RS210HBA']
    nearest_stations = ['RS508HBA', 'RS406HBA', 'RS306HBA', 'RS406HBA', 'RS307HBA', 'RS208HBA']

    return (distant_stations, nearest_stations)


def bootstrap_phases_to_reference(phase, station_names, ant_axis_ind, ref_id=0):
    """
    Fit screens to phase solutions

    Parameters
    ----------
    phase : array
        Array of phases
    ref_id : int, optional
        Index of reference station

    """
    nstations = phase.shape[ant_axis_ind]
    naxes = len(phase.shape)
    distant_stations, nearest_stations = get_bootstrap_stations()

    # Reverse the order and adjust
    distant_stations.reverse()
    nearest_stations.reverse()
    for distant_station, nearest_station in zip(distant_stations, nearest_stations):
        idx = station_names.tolist().index(distant_station)
        indlist = [idx if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
        nearest_ind = station_names.tolist().index(nearest_station)
        ref_indlist = [nearest_ind if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
        phase[indlist] = phase[indlist] + phase[ref_indlist]

    return phase


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


def main(h5parmfile, starttime=None, ntimes=None, solsetname='sol000',
    tecsoltabname='tec000', errsoltabname='error000',
    scphsoltabname='phase000', outsoltabroot='_screensols', ref_id=0,
    bootstrap=False):
    """
    Fit screens to TEC solutions

    Parameters
    ----------
    h5parmfile : str
        Filename of h5parm
    starttime : str
        Start time in casacore MVTime format
    ntimes : int
        Number of times to fit
    solsetname : str, optional
        Name of solset
    tecsoltabname : str, optional
        Name of TEC soltab
    errsoltabname : str, optional
        Name of error soltab
    scphsoltabname : str, optional
        Name of scalarphase soltab
    outsoltabroot : str, optional
        Root name for output soltabs
    ref_id : int, optional
        Index of reference station
    bootstrap : bool, optional
        If True, use bootstraping of the reference station for the most distant
        stations

    """
    ref_id = int(ref_id)

    # Read in solutions
    H = h5parm(h5parmfile, readonly=False)
    solset = H.getSolset(solsetname)
    tecsoltab = solset.getSoltab(tecsoltabname)
#     tecsoltab.setSelection(time={'min':tecsoltab.time[0],'max':tecsoltab.time[100]})
#     soltab.setSelection(ant=['RS409HBA'], dir='[Patch_29]')
    tec_vals = np.array(tecsoltab.val)
    errsoltab = solset.getSoltab(errsoltabname)
    err_vals = np.array(errsoltab.val)
    source_names = tecsoltab.dir[:]
    times = tecsoltab.time[:]
    station_names = tecsoltab.ant[:]
    freq = tecsoltab.freq[0]
    if scphsoltabname in solset.getSoltabNames():
        scphsoltab = solset.getSoltab(scphsoltabname)
        scphase_vals = np.array(scphsoltab.val)
        print('fit_tec_screens.py: Using TEC and scalar phase values')
    else:
        scphase_vals = None
        print('fit_tec_screens.py: Using TEC values only')
    ant_ind = tecsoltab.getAxesNames().index('ant')

    # Calculate TEC
    tec = calculate_total_tec(tec_vals, scphase_vals, freq, ant_ind,
        ref_id=ref_id, bootstrap=bootstrap, station_names=station_names)
    dtec = np.ones(tec_vals.shape)

    # Remove jumps
    nstat = tec.shape[1]
    ndir = tec.shape[2]
    jump_val = 2*np.pi/8.4479745e9*freq
    for s in range(nstat):
        if s == 58:
            pool = multiprocessing.Pool()
            tec_pool = [tec[:, s, d, 0] for d in range(ndir)]
            results = pool.map(remove_jumps_pool, itertools.izip(tec_pool, itertools.repeat(jump_val), itertools.repeat(31)))
            pool.close()
            pool.join()

            for d, (tec_clean, dtec_clean) in enumerate(results):
                # put back the results
                tec[:, s, d, 0] = tec_clean
                dtec[:, s, d, 0] = dtec_clean

    remove_soltabs(solset, ['newtec000'])
    solset.makeSoltab('tec', 'newtec000',
            axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times,
            station_names, source_names, np.array([freq])], vals=tec, weights=dtec)
