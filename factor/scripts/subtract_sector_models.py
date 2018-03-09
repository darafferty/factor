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
import glob
import subprocess


def get_nchunks(msin, nsectors):
    tot_m, used_m, free_m = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
    msin_m = float(subprocess.check_output(['du','-sh', msin]).split()[0][:-1]) * 1000.0
    tot_required_m = msin_m * nsectors * 2.0
    nchunks = int(np.ceil(tot_required_m / free_m))
    return nchunks


def main(msin, model_suffix, msin_column='DATA', model_column='DATA',
         out_column='DATA', nr_outliers=0, use_compression=False, peel_outliers=False):
    """
    Subtract sector model data

    Parameters
    ----------
    ms1 : str
        Name of MS file from which subtraction will be done
    ms2 : str or list
        Name of MS file from which column 2 will be taken.
    column1 : str
        Name of column 1
    column2 : str
        Name of column 2
    column_out : str
        Name of output column (written to ms1)
    op : str, optional
        Operation to perform: 'add', 'subtract12', or 'subtract21'
    nr_outliers : int, optional
        Number of outlier sectors. The last nr_outliers files are assumed to be the
        outlier sectors
    use_compression : bool, optional
        If True, use Dysco compression
    peel_outliers : bool, optional
        If True, outliers are peeled before sector models are subtracted

    """
    if type(use_compression) is str:
        if use_compression.lower() == 'true':
            use_compression = True
        else:
            use_compression = False
    nr_outliers = int(nr_outliers)

    # Find the model data files
    i = 0
    model_list = []
    while True:
        matches = glob.glob('{0}{1}_{2}'.format(msin, model_suffix, i))
        if len(matches) == 0:
            break
        model_list.extend(matches)
        i += 1
    nsectors = len(model_list)
    if nsectors == 1 and nr_outliers == 1:
        # This means we have a single imaging sector and outlier sector, so duplicate
        # the outlier model so that it gets subtracted properly later
        model_list *= 2
    elif nsectors == 0:
        print('No model data found. Exiting...')
        sys.exit(1)
    print('subtract_sector_models: found {} model data files'.format(nsectors))

    # If outliers are to be peeled, do that first
    if peel_outliers:
        # Open input and output table
        tin = pt.table(msin, readonly=True, ack=False)
        msout = '{}_peeled'.format(msin)
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
        msin = msout
        nr_outliers = 0
        model_list = model_list[:-nr_outliers]
        nsectors = len(model_list)

    # Open input and output tables
    tin = pt.table(msin, readonly=True, ack=False)
    tout_list = []
    for i, msmod in enumerate(model_list):
        if nr_outliers > 0 and i == len(model_list)-nr_outliers:
            # Break so we don't open output tables for the outliers
            break
        msout = '{}_sub'.format(msmod)
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
            tout.flush()
    for tout in tout_list:
        tout.close()
    tin.close()


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
