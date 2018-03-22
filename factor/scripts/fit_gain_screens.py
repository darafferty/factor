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


def main(h5parmfile, starttime=None, ntimes=None, solsetname='sol000',
    tecsoltabname='tec000', errsoltabname='error000', outsoltabroot='_screensols',
    ref_id=0, fit_screens=False, calculate_weights=True):
    """
    Fit screens to gain solutions

    Parameters
    ----------
    h5parmfile : str
        Filename of h5parm
    starttime : str, optional
        Start time in casacore MVTime format (NYI)
    ntimes : int, optional
        Number of times to fit (NYI)
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
    """
    ref_id = int(ref_id)

    # Read in solutions
    H = h5parm(h5parmfile, readonly=False)
    solset = H.getSolset(solsetname)
    tecsoltab = solset.getSoltab(tecsoltabname)
    tec = np.array(tecsoltab.val)
    errsoltab = solset.getSoltab(errsoltabname)
    err_vals = np.array(errsoltab.val)
    source_names = tecsoltab.dir[:]
    times = tecsoltab.time[:]
    station_names = tecsoltab.ant[:]
    freq = tecsoltab.freq[0]
    ant_ind = tecsoltab.getAxesNames().index('ant')
    dtec = np.ones(tec.shape)

    # Remove jumps
    nstat = tec.shape[1]
    ndir = tec.shape[2]
    jump_val = 2*np.pi/8.4479745e9*freq
    for s in range(nstat):
        pool = multiprocessing.Pool()
        tec_pool = [tec[:, s, d, 0] for d in range(ndir)]
        results = pool.map(remove_jumps_pool, itertools.izip(tec_pool,
                           itertools.repeat(jump_val), itertools.repeat(31)))
        pool.close()
        pool.join()

        for d, (tec_clean, dtec_clean) in enumerate(results):
            # put back the results
            tec[:, s, d, 0] = tec_clean
            dtec[:, s, d, 0] = dtec_clean

    remove_soltabs(solset, ['screentec000'])
    solset.makeSoltab('tec', 'screentec000',
            axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times,
            station_names, source_names, np.array([freq])], vals=tec, weights=dtec)

    if fit_screens:
        # Rename fixed TEC soltab (otherwise it will be overwritten later)
        soltab = solset.getSoltab('screentec000')
        soltab.rename('fixedtec000')

        # Find weights
        if calculate_weights:
            operations.reweight.run(soltab, 'window', nmedian=3, nstddev=501)

        # Fit screens
        remove_soltabs(solset, ['tecscreen000', 'tecscreen000resid'])
        operations.stationscreen.run(soltab, 'tecscreen000', niter=1, nsigma=5,
            refAnt=ref_id, order=20, scale_order=False)

        # Calculate values from screens
        remove_soltabs(solset, ['tec_screensols000'])
        soltab = solset.getSoltab('tecscreen000')
        source_dict = solset.getSou()
        operations.screenvalues.run(soltab, source_dict, outsoltabroot)
        soltab = solset.getSoltab('tec{}'.format(outsoltabroot))
        soltab.rename('screentec000')
