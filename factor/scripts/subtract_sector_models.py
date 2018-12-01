#! /usr/bin/env python
"""
Script to subtract sector model data
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy as np
import sys
import os
import subprocess
from lofarpipe.support.data_map import DataMap


def get_nchunks(msin, nsectors):
    tot_m, used_m, free_m = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
    msin_m = float(subprocess.check_output(['du', '-sh', msin]).split()[0][:-1]) * 1000.0
    tot_required_m = msin_m * nsectors * 2.0
    nchunks = max(1, int(np.ceil(tot_required_m / free_m)))
    return nchunks


def main(msin, mapfile_dir, filename, msin_column='DATA', model_column='DATA',
         out_column='DATA', nr_outliers=0, use_compression=False, peel_outliers=False,
         make_residual_col=True):
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
    make_residual_col : bool, optional
        If True, make a RESIDUAL_DATA column by subtracting all sources
    """
    if isinstance(use_compression, basestring):
        if use_compression.lower() == 'true':
            use_compression = True
        else:
            use_compression = False
    if isinstance(peel_outliers, basestring):
        if peel_outliers.lower() == 'true':
            peel_outliers = True
        else:
            peel_outliers = False
    if isinstance(make_residual_col, basestring):
        if make_residual_col.lower() == 'true':
            make_residual_col = True
        else:
            make_residual_col = False
    nr_outliers = int(nr_outliers)

    # Get the model data filenames. We only use files that contain the root of
    # msin, so that models for other observations are not picked up
    if msin.endswith('_field'):
        msin_root = msin.rstrip('_field')
    else:
        msin_root = msin
    mapfile = os.path.join(mapfile_dir, filename)
    model_map = DataMap.load(mapfile)
    model_list = [item.file for item in model_map if msin_root in item.file]
    nsectors = len(model_list)
    if nsectors == 1 and nr_outliers == 1:
        # This means we have a single imaging sector and outlier sector, so duplicate
        # the outlier model so that it gets subtracted properly later
        model_list *= 2
    elif nsectors == 0:
        print('subtract_sector_models: No model data found. Exiting...')
        sys.exit(1)
    print('subtract_sector_models: Found {} model data files'.format(nsectors))

    # If outliers are to be peeled, do that first
    if peel_outliers:
        # Open input and output table
        tin = pt.table(msin, readonly=True, ack=False)
        msout = '{}_field'.format(msin)
        if not os.path.exists(msout):
            os.system('/bin/cp -r {0} {1}'.format(msin, msout))
        tout = pt.table(msout, readonly=False, ack=False)

        nchunks = get_nchunks(msin, nr_outliers)
        nrows_per_chunk = tin.nrows() / nchunks
        startrows = [0]
        nrows = [nrows_per_chunk]
        for i in range(1, nchunks):
            if i == nchunks-1:
                nrow = tin.nrows() - (nchunks - 1)*nrows_per_chunk
            else:
                nrow = nrows_per_chunk
            nrows.append(nrow)
            startrows.append(startrows[i-1] + nrows[i-1])
        print('subtract_sector_models: Using {} chunks for peeling'.format(nchunks))

        for c, (startrow, nrow) in enumerate(zip(startrows, nrows)):
            # For each chunk, load data
            print('subtract_sector_models: Processing chunk {}...'.format(c))
            datain = tin.getcol(msin_column, startrow=startrow, nrow=nrow)
            if use_compression:
                # Replace flagged values with NaNs before compression
                flags = tin.getcol('FLAG', startrow=startrow, nrow=nrow)
                flagged = np.where(flags)
                datain[flagged] = np.NaN
            datamod_list = []
            for i, msmodel in enumerate(model_list[nsectors-nr_outliers:]):
                print('subtract_sector_models: Loading data for outlier sector {}...'.format(i))
                tmod = pt.table(msmodel, readonly=True, ack=False)
                datamod_list.append(tmod.getcol(model_column, startrow=startrow, nrow=nrow))
                tmod.close()

            # For each sector, subtract sum of model data for this chunk
            print('subtract_sector_models: Peeling data...')
            other_sectors_ind = list(range(nr_outliers))
            datamod_all = None
            for sector_ind in other_sectors_ind:
                if datamod_all is None:
                    datamod_all = datamod_list[sector_ind].copy()
                else:
                    datamod_all += datamod_list[sector_ind]
            tout.putcol(out_column, datain-datamod_all, startrow=startrow, nrow=nrow)
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
            _modeldata
        msout = msmod.rstrip('_modeldata')
        if not os.path.exists(msout):
            os.system('/bin/cp -r {0} {1}'.format(msin, msout))
        tout_list.append(pt.table(msout, readonly=False, ack=False))

    # Define chunks based on available memory
    nchunks = get_nchunks(msin, nsectors)
    nrows_per_chunk = tin.nrows() / nchunks
    startrows = [0]
    nrows = [nrows_per_chunk]
    for i in range(1, nchunks):
        if i == nchunks-1:
            nrow = tin.nrows() - (nchunks - 1)*nrows_per_chunk
        else:
            nrow = nrows_per_chunk
        nrows.append(nrow)
        startrows.append(startrows[i-1] + nrows[i-1])
    print('subtract_sector_models: Using {} chunks'.format(nchunks))

    for c, (startrow, nrow) in enumerate(zip(startrows, nrows)):
        # For each chunk, load data
        print('subtract_sector_models: Processing chunk {}...'.format(c))
        datain = tin.getcol(msin_column, startrow=startrow, nrow=nrow)
        if use_compression:
            # Replace flagged values with NaNs before compression
            flags = tin.getcol('FLAG', startrow=startrow, nrow=nrow)
            flagged = np.where(flags)
            datain[flagged] = np.NaN
        datamod_list = []
        for i, msmodel in enumerate(model_list):
            print('subtract_sector_models: Loading data for sector {}...'.format(i))
            tmod = pt.table(msmodel, readonly=True, ack=False)
            datamod_list.append(tmod.getcol(model_column, startrow=startrow, nrow=nrow))
            tmod.close()

        # For each sector, subtract sum of model data for this chunk
        for i, tout in enumerate(tout_list):
            print('subtract_sector_models: Subtracting data for sector {}...'.format(i))
            other_sectors_ind = list(range(nsectors))
            other_sectors_ind.pop(i)
            datamod_all = None
            for sector_ind in other_sectors_ind:
                if datamod_all is None:
                    datamod_all = datamod_list[sector_ind].copy()
                else:
                    datamod_all += datamod_list[sector_ind]
            tout.putcol(out_column, datain-datamod_all, startrow=startrow, nrow=nrow)
            if make_residual_col:
                # Also subtract sector's model and output to RESIDUAL_DATA column
                datamod_all += datamod_list[i]
                if not 'RESIDUAL_DATA' in tout.colnames():
                    desc = tout.getcoldesc('DATA')
                    desc['name'] = 'RESIDUAL_DATA'
                    coldmi = tout.getdminfo('DATA')
                    coldmi["NAME"] = 'RESIDUAL_DATA'
                    tout.addcols(desc, coldmi)
                tout.putcol('RESIDUAL_DATA', datain-datamod_all, startrow=startrow, nrow=nrow)
            tout.flush()
    for tout in tout_list:
        tout.close()
    tin.close()

    # Delete model data
    for msmod in model_list:
        if os.path.exists(msmod):
            os.system('/bin/rm -rf {0}'.format(msmod))


if __name__ == '__main__':
    descriptiontext = "Add/subtract columns (column_out = column1 +/- column2).\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms1', help='name of MS file 1')
    parser.add_argument('ms2', help='name of MS file 2')
    parser.add_argument('column1', help='name of column 1')
    parser.add_argument('column2', help='name of column 2')
    parser.add_argument('column_out', help='name of the output column (written to ms1)')
    parser.add_argument('op', help='operation: "add" or "subtract"')
    parser.add_argument('in_memory', help='do operation in memory')
    args = parser.parse_args()

    main(args.ms1, args.ms2, args.column1, args.column2, args.column_out, args.op, args.in_memory)
