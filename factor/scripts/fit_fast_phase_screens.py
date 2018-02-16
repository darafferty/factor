#! /usr/bin/env python
"""
Script to fit fast-phase solutions with screens
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


def calculate_total_phase(tec_vals, scphase_vals, freq, ant_axis_ind, ref_id=0,
    bootstrap=False, station_names=None, baseline_file=None):

    # Calculate total phase
    r = (-8.4479745e9 * tec_vals/freq) + scphase_vals

    # Adjust all to reference station
    nstations = tec_vals.shape[ant_axis_ind]
    naxes = len(tec_vals.shape)
    for idx in range(nstations):
        indlist = [idx if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
        ref_indlist = [ref_id if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
        r[indlist] = r[indlist] - r[ref_indlist]
    if bootstrap:
        # For distant stations, use the nearest as the reference
        distant_stations, nearest_stations = get_bootstrap_stations()
        for distant_station, nearest_station in zip(distant_stations, nearest_stations):
            idx = station_names.tolist().index(distant_station)
            indlist = [idx if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
            nearest_ind = station_names.tolist().index(nearest_station)
            ref_indlist = [nearest_ind if i == ant_axis_ind else slice(0, None) for i in range(naxes)]
            r[indlist] = r[indlist] - r[ref_indlist]

    return r


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


def main(h5parmfile, freq1=120e6, freq2=170e6, tstart=None,
    solsetname='sol000', tecsoltabname='tec000',
    scphsoltabname='scalarphase000', outsoltabroot='_screensols', ref_id=0,
    bootstrap=False, calculate_weights=False):
    """
    Fit screens to phase solutions

    Parameters
    ----------
    h5parmfile : str
        Filename of h5parm
    freq1 : float
        Low frequency in Hz at which screens will be fit
    freq1 : float
        High frequency in Hz at which screens will be fit
    starttime : str
        Start time in casacore MVTime format
    ntimes : int
        Number of times to fit
    solsetname : str, optional
        Name of solset
    tecsoltabname : str, optional
        Name of TEC soltab
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
    scphsoltab = solset.getSoltab(scphsoltabname)
    tec_vals = np.array(tecsoltab.val)
    scphase_vals = np.array(scphsoltab.val)
    source_names = tecsoltab.dir[:]
    times = tecsoltab.time[:]
    station_names = tecsoltab.ant[:]

    # Find total phases at two frequencies
    freqs = [freq1, freq2]
    ant_ind = tecsoltab.getAxesNames().index('ant')
    for i, freq in enumerate(freqs):
        remove_soltabs(solset, ['totalphase{}'.format(i+1)])
        r = calculate_total_phase(tec_vals, scphase_vals, freq, ant_ind,
            ref_id=ref_id, bootstrap=bootstrap, station_names=station_names)
        w = np.ones(r.shape)
        solset.makeSoltab('phase', 'totalphase{}'.format(i+1),
                axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times,
                station_names, source_names, np.array([freq])], vals=r, weights=w)

        # Find weights
        if calculate_weights:
            soltab = solset.getSoltab('totalphase{}'.format(i+1))
            operations.reweight.run(soltab, 'window', nmedian=3, nstddev=501)

        # Fit screens
        remove_soltabs(solset, ['phasescreen{}'.format(i+1), 'phasescreen{}resid'.format(i+1)])
        soltab = solset.getSoltab('totalphase{}'.format(i+1))
        allbut1 = [d for d in soltab.dir[:] if d != '[Patch_1]']
        soltab.setSelection(dir=allbut1, time={'min': 0.0, 'max': 4987956330.1176252})

        operations.stationscreen.run(soltab, 'phasescreen{}'.format(i+1), niter=1, nsigma=5,
            refAnt=ref_id, order=20, scale_order=False)

    # Plot
	distant_stations, nearest_stations = get_bootstrap_stations()
	soltab = solset.getSoltab('phasescreen1')
	soltab.setSelection(ant=distant_stations[5])
	ressoltab = solset.getSoltab('phasescreen1resid')
	ressoltab.setSelection(ant=distant_stations[5])
	if bootstrap:
		prefix = 'bootstrap'
	else:
		prefix = 'rs210_120'
	operations.plotscreen.run(soltab, ressoltab=ressoltab, prefix=prefix)
	soltab = solset.getSoltab('phasescreen2')
	soltab.setSelection(ant=distant_stations[5])
	ressoltab = solset.getSoltab('phasescreen2resid')
	ressoltab.setSelection(ant=distant_stations[5])
	if bootstrap:
		prefix = 'bootstrap'
	else:
		prefix = 'rs210_170'
	operations.plotscreen.run(soltab, ressoltab=ressoltab, prefix=prefix)

    # Bootstrap the screens back to the global reference
    if bootstrap:
        soltab = solset.getSoltab('phasescreen1')
        ant_ind = soltab.getAxesNames().index('ant')
        phase = bootstrap_phases_to_reference(soltab.val[:], station_names,
            ant_ind, ref_id)
        soltab.setValues(phase)
        res_soltab = solset.getSoltab('phasescreen1resid')
        resid = bootstrap_phases_to_reference(res_soltab.val[:], station_names,
            ant_ind, ref_id)
        res_soltab.setValues(resid)
        soltab = solset.getSoltab('phasescreen2')
        phase = bootstrap_phases_to_reference(soltab.val[:], station_names,
            ant_ind, ref_id)
        soltab.setValues(phase)
        res_soltab = solset.getSoltab('phasescreen2resid')
        resid = bootstrap_phases_to_reference(res_soltab.val[:], station_names,
            ant_ind, ref_id)
        res_soltab.setValues(resid)

    # Calculate values from screens
    remove_soltabs(solset, ['tec_screensols000', 'scalarphase_screensols000'])
    soltab_names = ['phasescreen{}'.format(i+1) for i in range(len(freqs))]
    soltabs = [solset.getSoltab(n) for n in soltab_names]
    tecsoltab.setSelection(dir=allbut1, time={'min': 0.0, 'max': 4987956330.1176252})
    soltabs.append(tecsoltab)
    source_dict = solset.getSou()
    operations.screenvalues.run(soltabs, source_dict, '_screensols000')
