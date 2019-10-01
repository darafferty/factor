#! /usr/bin/env python
"""
Script to subtract sector model data
"""
import casacore.tables as pt
import numpy as np
import sys
import os
import subprocess
from lofarpipe.support.data_map import DataMap
from astropy.time import Time
import dateutil.parser
import losoto
from factor.lib import miscellaneous as misc


def get_nchunks(msin, nsectors, fraction=1.0):
    """
    Determines number of chunks for available memory of node

    Parameters
    ----------
    msin : str
        Input MS file name
    nsectors : int
        Number of imaging sectors
    fraction : float
        Fraction of total memory to use

    Returns
    -------
    nchunks : int
        Number of chunks
    """
    tot_m, used_m, free_m = list(map(int, os.popen('free -tm').readlines()[-1].split()[1:]))
    msin_m = float(subprocess.check_output(['du', '-sm', msin]).split()[0]) * fraction
    tot_required_m = msin_m * nsectors * 2.0
    nchunks = max(1, int(np.ceil(tot_required_m / tot_m)))
    return nchunks


def convert_mjd2mvt(mjd_sec):
    """
    Converts MJD to casacore MVTime

    Parameters
    ----------
    mjd_sec : float
        MJD time in seconds

    Returns
    -------
    mvtime : str
        Casacore MVTime string
    """
    t = Time(mjd_sec / 3600 / 24, format='mjd', scale='utc')
    date, hour = t.iso.split(' ')
    year, month, day = date.split('-')
    d = t.datetime
    month = d.ctime().split(' ')[1]

    return '{0}{1}{2}/{3}'.format(day, month, year, hour)


def convert_mvt2mjd(mvt_str):
    """
    Converts casacore MVTime to MJD

    Parameters
    ----------
    mvt_str : str
        MVTime time

    Returns
    -------
    mjdtime : float
        MJD time in seconds
    """
    day_str = mvt_str.split('/')[0].lower()
    months = {'jan':1, 'feb':2, 'mar':3, 'apr':4, 'may':5, 'jun':6,
              'jul':7, 'aug':8, 'sep':9, 'oct':10, 'nov':11, 'dec':12}
    for m in months:
        if m in day_str:
            p = day_str.split(m)
            day_str = '{0}.{1}.{2}'.format(p[1], months[m], p[0])
    time_str = mvt_str.split('/')[1]
    t = dateutil.parser.parse('{0}/{1}'.format(day_str, time_str))

    return Time(t).mjd * 3600 * 24


def main(msin, mapfile_dir, filename, msin_column='DATA', model_column='DATA',
         out_column='DATA', nr_outliers=0, use_compression=False, peel_outliers=False,
         reweight=True, starttime=None, solint_sec=None, solint_hz=None,
         weights_colname="CAL_WEIGHT", gainfile="", uvcut_min=80.0, uvcut_max=1e6,
         phaseonly=True, dirname=None, quiet=True):
    """
    Subtract sector model data

    Parameters
    ----------
    msin : str
        Name of MS file from which subtraction will be done
    mapfile_dir : str
        Path of current pipeline mapfile directory
    filename: str
        Name of mapfile containing model data filenames
    msin_column : str, optional
        Name of input column
    model_column : str, optional
        Name of model column
    out_column : str, optional
        Name of output column
    nr_outliers : int, optional
        Number of outlier sectors. The last nr_outliers files are assumed to be the
        outlier sectors
    use_compression : bool, optional
        If True, use Dysco compression
    peel_outliers : bool, optional
        If True, outliers are peeled before sector models are subtracted
    starttime : str, optional
        Start time in JD seconds
    """
    use_compression = misc.string2bool(use_compression)
    peel_outliers = misc.string2bool(peel_outliers)
    nr_outliers = int(nr_outliers)
    solint_sec = float(solint_sec)
    solint_hz = float(solint_hz)
    uvcut_min = float(uvcut_min)
    uvcut_max = float(uvcut_max)
    uvcut = [uvcut_min, uvcut_max]
    phaseonly = misc.string2bool(phaseonly)

    # Get the model data filenames. We only use files that contain the root of
    # msin, so that models for other observations are not picked up (starttime
    # is also used when a single MS file is used for multiple observations)
    if msin.endswith('_field'):
        msin_root = msin.rstrip('_field')
    else:
        msin_root = msin
    mapfile = os.path.join(mapfile_dir, filename)
    model_map = DataMap.load(mapfile)
    model_list = [item.file for item in model_map if msin_root in item.file]
    if starttime is not None:
        # Filter the list of models to include only ones for the given times
        nrows_list = []
        for msmod in model_list[:]:
            tin = pt.table(msmod, readonly=True, ack=False)
            starttime_chunk = np.min(tin.getcol('TIME'))
            if not misc.approx_equal(starttime_chunk, convert_mvt2mjd(starttime), tol=1.0):
                # Remove files with start times that are not within 1 sec of the
                # specified starttime
                i = model_list.index(msmod)
                model_list.pop(i)
            else:
                nrows_list.append(tin.nrows())
                starttime_exact = convert_mjd2mvt(starttime_chunk)  # store exact value for use later
            tin.close()
        if len(set(nrows_list)) > 1:
            print('subtract_sector_models: Model data files have differing number of rows...')
            sys.exit(1)

    nsectors = len(model_list)
    if nsectors == 1 and nr_outliers == 1:
        # This means we have a single imaging sector and outlier sector, so duplicate
        # the outlier model so that it gets subtracted properly later
        model_list *= 2
    elif nsectors == 0:
        print('subtract_sector_models: No model data found. Exiting...')
        sys.exit(1)
    print('subtract_sector_models: Found {} model data files'.format(nsectors))

    # If starttime is given, figure out startrow and nrows for input MS file
    tin = pt.table(msin, readonly=True, ack=False)
    tarray = tin.getcol("TIME")
    nbl = np.where(tarray == tarray[0])[0].size
    tarray = None
    if starttime is not None:
        times = [convert_mjd2mvt(t) for t in tin.getcol('TIME')]
        startrow_in = times.index(starttime_exact)
        nrows_in = nrows_list[0]
    else:
        startrow_in = 0
        nrows_in = tin.nrows()
    tin.close()

    # If outliers are to be peeled, do that first
    if peel_outliers and nr_outliers > 0:
        # Open input and output table
        tin = pt.table(msin, readonly=True, ack=False)
        msout = '{}_field'.format(msin)
        if not os.path.exists(msout):
            os.system('/bin/cp -r {0} {1}'.format(msin, msout))
        tout = pt.table(msout, readonly=False, ack=False)

        # Define chunks based on available memory
        fraction = float(nrows_in) / float(tin.nrows())
        nchunks = get_nchunks(msin, nr_outliers, fraction)
        nrows_per_chunk = int(nrows_in / nchunks)
        startrows_tin = [startrow_in]
        startrows_tmod = [0]
        nrows = [nrows_per_chunk]
        for i in range(1, nchunks):
            if i == nchunks-1:
                nrow = nrows_in - (nchunks - 1) * nrows_per_chunk
            else:
                nrow = nrows_per_chunk
            nrows.append(nrow)
            startrows_tin.append(startrows_tin[i-1] + nrows[i-1])
            startrows_tmod.append(startrows_tmod[i-1] + nrows[i-1])
        print('subtract_sector_models: Using {} chunk(s) for peeling'.format(nchunks))

        for c, (startrow_tin, startrow_tmod, nrow) in enumerate(zip(startrows_tin, startrows_tmod, nrows)):
            # For each chunk, load data
            datain = tin.getcol(msin_column, startrow=startrow_tin, nrow=nrow)
            if use_compression:
                # Replace flagged values with NaNs before compression
                flags = tin.getcol('FLAG', startrow=startrow_tin, nrow=nrow)
                flagged = np.where(flags)
                datain[flagged] = np.NaN
            datamod_list = []
            for i, msmodel in enumerate(model_list[nsectors-nr_outliers:]):
                tmod = pt.table(msmodel, readonly=True, ack=False)
                datamod_list.append(tmod.getcol(model_column, startrow=startrow_tmod, nrow=nrow))
                tmod.close()

            # For each sector, subtract sum of model data for this chunk
            other_sectors_ind = list(range(nr_outliers))
            datamod_all = None
            for sector_ind in other_sectors_ind:
                if datamod_all is None:
                    datamod_all = datamod_list[sector_ind].copy()
                else:
                    datamod_all += datamod_list[sector_ind]
            tout.putcol(out_column, datain-datamod_all, startrow=startrow_tmod, nrow=nrow)
            tout.flush()
        tout.close()
        tin.close()

        # Now reset things for the imaging sectors
        msin = msout
        model_list = model_list[:-nr_outliers]
        nsectors = len(model_list)
        nr_outliers = 0

    # Open input and output tables
    tin = pt.table(msin, readonly=True, ack=False)
    tout_list = []
    for i, msmod in enumerate(model_list):
        if nr_outliers > 0 and i == len(model_list)-nr_outliers:
            # Break so we don't open output tables for the outliers
            break
        msout = msmod.rstrip('_modeldata')
        if not os.path.exists(msout):
            os.system('/bin/cp -r {0} {1}'.format(msin, msout))
        tout_list.append(pt.table(msout, readonly=False, ack=False))

    # Define chunks based on available memory, making sure each chunk gives a full
    # timeslot (needed for reweighting)
    fraction = float(nrows_in) / float(tin.nrows())
    nchunks = get_nchunks(msin, nsectors, fraction)
    nrows_per_chunk = int(nrows_in / nchunks)
    while nrows_per_chunk % nbl > 0.0:
        nrows_per_chunk -= 1
        if nrows_per_chunk < nbl:
            nrows_per_chunk = nbl
            break
    nchunks = int(np.ceil(nrows_in / nrows_per_chunk))
    startrows_tin = [startrow_in]
    startrows_tmod = [0]
    nrows = [nrows_per_chunk]
    for i in range(1, nchunks):
        if i == nchunks-1:
            nrow = nrows_in - (nchunks - 1) * nrows_per_chunk
        else:
            nrow = nrows_per_chunk
        nrows.append(nrow)
        startrows_tin.append(startrows_tin[i-1] + nrows[i-1])
        startrows_tmod.append(startrows_tmod[i-1] + nrows[i-1])
    print('subtract_sector_models: Using {} chunk(s)'.format(nchunks))

    for c, (startrow_tin, startrow_tmod, nrow) in enumerate(zip(startrows_tin, startrows_tmod, nrows)):
        # For each chunk, load data
        datain = tin.getcol(msin_column, startrow=startrow_tin, nrow=nrow)
        flags = tin.getcol('FLAG', startrow=startrow_tin, nrow=nrow)
        datamod_list = []
        for i, msmodel in enumerate(model_list):
            tmod = pt.table(msmodel, readonly=True, ack=False)
            datamod_list.append(tmod.getcol(model_column, startrow=startrow_tmod, nrow=nrow))
            tmod.close()

        # For each sector, subtract sum of model data for this chunk
        weights = None
        for i, tout in enumerate(tout_list):
            other_sectors_ind = list(range(nsectors))
            other_sectors_ind.pop(i)
            datamod_all = None
            for sector_ind in other_sectors_ind:
                if datamod_all is None:
                    datamod_all = datamod_list[sector_ind].copy()
                else:
                    datamod_all += datamod_list[sector_ind]
            tout.putcol(out_column, datain-datamod_all, startrow=startrow_tmod, nrow=nrow)
            if reweight:
                # Also subtract sector's model to make residual data for reweighting
                if weights is None:
                    datamod_all += datamod_list[i]
                    covweights = CovWeights(model_list[0], solint_sec, solint_hz, startrow_tmod, nrow,
                                            gainfile=gainfile, uvcut=uvcut, phaseonly=phaseonly,
                                            dirname=dirname, quiet=quiet)
                    coefficients = covweights.FindWeights(datain-datamod_all, flags)
                    weights = covweights.calcWeights(coefficients)
                tout.putcol(weights_colname, weights, startrow=startrow_tmod, nrow=nrow)
            tout.flush()
    for tout in tout_list:
        tout.close()
    tin.close()

    # Delete model data
    for msmod in model_list:
        if os.path.exists(msmod):
            os.system('/bin/rm -rf {0}'.format(msmod))


class CovWeights:
    def __init__(self, MSName, solint_sec, solint_hz, startrow, nrow, uvcut=[0, 2000],
                 gainfile=None, phaseonly=False, dirname=None, quiet=True):
        if MSName[-1] == "/":
            self.MSName = MSName[0:-1]
        else:
            self.MSName = MSName
        tab = pt.table(self.MSName, ack=False)
        self.timepersample = tab.getcell('EXPOSURE', 0)
        self.ntSol = max(1, int(round(solint_sec / self.timepersample)))
        tab.close()
        sw = pt.table(self.MSName+'::SPECTRAL_WINDOW', ack=False)
        self.referencefreq = sw.col('REF_FREQUENCY')[0]
        self.channelwidth = sw.col('CHAN_WIDTH')[0][0]
        self.numchannels = sw.col('NUM_CHAN')[0]
        sw.close()
        self.nchanSol = max(1, self.get_nearest_frequstep(solint_hz / self.channelwidth))
        self.uvcut = uvcut
        self.gainfile = gainfile
        self.phaseonly = phaseonly
        self.dirname = dirname
        self.quiet = quiet
        self.startrow = startrow
        self.nrow = nrow

    def FindWeights(self, residualdata, flags):
        ms = pt.table(self.MSName, ack=False)
        ants = pt.table(ms.getkeyword("ANTENNA"), ack=False)
        antnames = ants.getcol("NAME")
        ants.close()
        nAnt = len(antnames)

        u, v, _ = ms.getcol("UVW", startrow=self.startrow, nrow=self.nrow).T
        A0 = ms.getcol("ANTENNA1", startrow=self.startrow, nrow=self.nrow)
        A1 = ms.getcol("ANTENNA2", startrow=self.startrow, nrow=self.nrow)
        tarray = ms.getcol("TIME", startrow=self.startrow, nrow=self.nrow)
        nbl = np.where(tarray == tarray[0])[0].size
        ms.close()

        # apply uvcut
        c_m_s = 2.99792458e8
        uvlen = np.sqrt(u**2 + v**2) / c_m_s * self.referencefreq
        flags[uvlen > self.uvcut[1], :, :] = True
        flags[uvlen < self.uvcut[0], :, :] = True
        residualdata[flags] = np.nan
        residualdata[residualdata == 0] = np.nan

        # initialise
        nChan = residualdata.shape[1]
        nPola = residualdata.shape[2]
        nt = int(residualdata.shape[0] / nbl)
        residualdata = residualdata.reshape((nt, nbl, nChan, nPola))
        A0 = A0.reshape((nt, nbl))[0, :]
        A1 = A1.reshape((nt, nbl))[0, :]

        # make rms and residuals arrays
        rmsarray = np.zeros((nt, nbl, nChan, 2), dtype=np.complex64)
        residuals = np.zeros_like(rmsarray, dtype=np.complex64)
        rmsarray[:, :, :, 0] = residualdata[:, :, :, 1]
        rmsarray[:, :, :, 1] = residualdata[:, :, :, 2]
        residuals[:, :, :, 0] = residualdata[:, :, :, 0]
        residuals[:, :, :, 1] = residualdata[:, :, :, 3]

        # start calculating the weights
        CoeffArray = np.zeros((nt, nChan, nAnt))
        ant1 = np.arange(nAnt)
        tcellsize = self.ntSol
        for t_i in range(0, nt, self.ntSol):
            if (t_i == nt - self.ntSol) and (nt % self.ntSol > 0):
                tcellsize = nt % self.ntSol
            t_e = t_i + tcellsize
            fcellsize = self.nchanSol
            for f_i in range(0, nChan, self.nchanSol):
                if (f_i == nChan - self.nchanSol) and (nChan % self.nchanSol > 0):
                    fcellsize = nChan % self.nchanSol
                f_e = f_i + fcellsize

                # build weights for each antenna in the current time-frequency block
                for ant in ant1:
                    # set of vis for baselines ant-ant_i
                    set1 = np.where(A0 == ant)[0]
                    # set of vis for baselines ant_i-ant
                    set2 = np.where(A1 == ant)[0]
                    CoeffArray[t_i:t_e, f_i:f_e, ant] = np.sqrt(
                        np.nanmean(
                            np.append(residuals[t_i:t_e, set1, f_i:f_e, :],
                                      residuals[t_i:t_e, set2, f_i:f_e, :]) *
                            np.append(residuals[t_i:t_e, set1, f_i:f_e, :],
                                      residuals[t_i:t_e, set2, f_i:f_e, :]).conj()
                            ) -
                        np.nanstd(
                            np.append(rmsarray[t_i:t_e, set1, f_i:f_e, :],
                                      rmsarray[t_i:t_e, set2, f_i:f_e, :])
                            )
                        )

        # get rid of NaNs and low values
        CoeffArray[np.isnan(CoeffArray)] = np.inf
        for i in range(nAnt):
            tempars = CoeffArray[:, :, i]
            thres = 0.25 * np.median(tempars[np.where(np.isfinite(tempars))])
            CoeffArray[:, :, i][tempars < thres] = thres
        return CoeffArray

    def calcWeights(self, CoeffArray):
        ms = pt.table(self.MSName, readonly=False, ack=False)
        ants = pt.table(ms.getkeyword("ANTENNA"), ack=False)
        antnames = ants.getcol("NAME")
        nAnt = len(antnames)
        tarray = ms.getcol("TIME", startrow=self.startrow, nrow=self.nrow)
        darray = ms.getcol("DATA", startrow=self.startrow, nrow=self.nrow)
        tvalues = np.array(sorted(list(set(tarray))))
        nt = tvalues.shape[0]
        nbl = int(tarray.shape[0]/nt)
        nchan = darray.shape[1]
        npol = darray.shape[2]
        A0 = np.array(ms.getcol("ANTENNA1", startrow=self.startrow, nrow=self.nrow).reshape((nt, nbl)))
        A1 = np.array(ms.getcol("ANTENNA2", startrow=self.startrow, nrow=self.nrow).reshape((nt, nbl)))

        # initialize weight array
        w = np.zeros((nt, nbl, nchan, npol))
        A0ind = A0[0, :]
        A1ind = A1[0, :]

        # do gains stuff
        ant1gainarray, ant2gainarray = readGainFile(self.gainfile, ms, nt, nchan, nbl,
                                                    tarray, nAnt, self.MSName, self.phaseonly,
                                                    self.dirname)
        ant1gainarray = ant1gainarray.reshape((nt, nbl, nchan))
        ant2gainarray = ant2gainarray.reshape((nt, nbl, nchan))
        for t in range(nt):
            for i in range(nbl):
                for j in range(nchan):
                    w[t, i, j, :] = 1.0 / (CoeffArray[t, j, A0ind[i]] * ant1gainarray[t, i, j] +
                                           CoeffArray[t, j, A1ind[i]] * ant2gainarray[t, i, j] +
                                           CoeffArray[t, j, A0ind[i]] * CoeffArray[t, j, A1ind[i]] +
                                           0.1)

        # normalize
        w = w.reshape(nt*nbl, nchan, npol)
        w[np.isinf(w)] = np.nan
        w = w / np.nanmean(w)
        w[np.isnan(w)] = 0

        return w

    def get_nearest_frequstep(self, freqstep):
        """
        Gets the nearest frequstep

        Parameters
        ----------
        freqstep : int
            Target frequency step

        Returns
        -------
        optimum_step : int
            Optimum frequency step nearest to target step
        """
        # Generate a list of possible values for freqstep
        if not hasattr(self, 'freq_divisors'):
            tmp_divisors = []
            for step in range(self.numchannels, 0, -1):
                if (self.numchannels % step) == 0:
                    tmp_divisors.append(step)
            self.freq_divisors = np.array(tmp_divisors)

        # Find nearest
        idx = np.argmin(np.abs(self.freq_divisors - freqstep))

        return self.freq_divisors[idx]


def readGainFile(gainfile, ms, nt, nchan, nbl, tarray, nAnt, msname, phaseonly, dirname):
    if phaseonly:
        ant1gainarray1 = np.ones((nt*nbl, nchan))
        ant2gainarray1 = np.ones((nt*nbl, nchan))
    else:
        import losoto.h5parm
        solsetName = "sol000"
        soltabName = "screenamplitude000"
        try:
            gfile = losoto.h5parm.openSoltab(gainfile, solsetName=solsetName, soltabName=soltabName)
        except:
            print("Could not find amplitude gains in h5parm. Assuming gains of 1 everywhere.")
            ant1gainarray1 = np.ones((nt*nbl, nchan))
            ant2gainarray1 = np.ones((nt*nbl, nchan))
            return ant1gainarray1, ant2gainarray1

        freqs = pt.table(msname+"/SPECTRAL_WINDOW").getcol("CHAN_FREQ")
        gains = gfile.val  # axes: times, freqs, ants, dirs, pols
        flagged = np.where(gains == 0.0)
        gains[flagged] = np.nan
        gfreqs = gfile.freq
        times = gfile.time
        dindx = gfile.dir.tolist().index(dirname)
        ant1gainarray = np.zeros((nt*nbl, nchan))
        ant2gainarray = np.zeros((nt*nbl, nchan))
        A0arr = ms.getcol("ANTENNA1", startrow=self.startrow, nrow=self.nrow)
        A1arr = ms.getcol("ANTENNA2", startrow=self.startrow, nrow=self.nrow)
        deltime = (times[1] - times[0]) / 2.0
        delfreq = (gfreqs[1] - gfreqs[0]) / 2.0
        for i in range(len(times)):
            timemask = (tarray >= times[i]-deltime) * (tarray < times[i]+deltime)
            if np.all(~timemask):
                continue
            for j in range(nAnt):
                mask1 = timemask * (A0arr == j)
                mask2 = timemask * (A1arr == j)
                for k in range(nchan):
                    chan_freq = freqs[0, k]
                    freqmask = np.logical_and(gfreqs >= chan_freq-delfreq,
                                              gfreqs < chan_freq+delfreq)
                    if chan_freq < gfreqs[0]:
                        freqmask[0] = True
                    if chan_freq > gfreqs[-1]:
                        freqmask[-1] = True
                    ant1gainarray[mask1, k] = np.nanmean(gains[i, freqmask, j, dindx, :], axis=(0, 1))
                    ant2gainarray[mask2, k] = np.nanmean(gains[i, freqmask, j, dindx, :], axis=(0, 1))
        ant1gainarray1 = ant1gainarray**2
        ant2gainarray1 = ant2gainarray**2

    return ant1gainarray1, ant2gainarray1
