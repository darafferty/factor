#! /usr/bin/env python
"""
Script to combine multiple h5parms in time or frequency
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os
from losoto.h5parm import h5parm


def main(input_files, output_file, soltype, solsetname='sol000'):
    """
    Combine multiple h5parms in time or frequency

    Parameters
    ----------
    input_files : str
        Filenames of input h5parms
    output_file : str
        Filename of output h5parm
    soltype : str
        One of 'tecandphase' or 'complexgain'
    solsetname : str
        Name of solset (must be the same in every input file)

    """
    input_files = input_files.strip('[]').split(',')
    input_files  = [f.strip() for f in input_files]

    # Get values
    vals1 = []
    vals2 = []
    times = []
    starttimes = []
    for h5file in input_files:
        H = h5parm(h5file, readonly=True)
        solset = H.getSolset(solsetname)
        if type == 'tecandphase':
            soltab1 = solset.getSoltab('tec000')
            soltab2 = solset.getSoltab('scalarphase000')
        elif type == 'complexgain':
            soltab1 = solset.getSoltab('phase000')
            soltab2 = solset.getSoltab('amplitude000')
        vals1.append(soltab1.vals[:])
        vals2.append(soltab2.vals[:])
        times.append(soltab1.time[:])
        starttimes.append(soltab1.time[0])

    # Sort by time
    sorted_ind = np.argsort(np.array(starttimes))[0]
    vals1 = [vals1[i] for i in sorted_ind]
    vals2 = [vals2[i] for i in sorted_ind]
    times = [times[i] for i in sorted_ind]

    # Combine values and weights
    allvals1 = []
    allweights1 = []
    allvals2 = []
    allweights2 = []
    alltimes = []
    for v1, v2, w1, w2, t in zip(vals1, vals2, weights1, weights2, times):
        allvals1.extend(v1)
        allvals2.extend(v2)
        allweights1.extend(w1)
        allweights2.extend(w2)
        alltimes.extend(t)

    # Get values for other axes
    axesNames = soltab1.getAxesNames()
    axesVals = []
    for axis in axesNames:
        axesVals.append(soltab1.getAxesValues(axis))

    # Create output h5parm
    h5Out = h5parm(output_file, readonly=False)

    # Create output solset
    solset = h5Out.makeSolset('sol000')
    sourceTable = solset._f_get_child('source')
    antennaTable = solset._f_get_child('antenna')
    names = []
    positions = []
    ant_dict = solset.getAnt()
    for k, v in ant_dict.iteritems():
        names.append(k)
        positions.append(v.tolist())
    antennaTable.append(zip(*(ant.keys(), ant.values())))
    names = []
    positions = []
    source_dict = solset.getSou()
    for k, v in source_dict.iteritems():
        names.append(k)
        positions.append(v.tolist())
    sourceTable.append(zip(*(names, positions)))

    # Create output soltabs
    if type == 'tecandphase':
        solset.makeSoltab('sol000', 'tec', 'tec000', axesNames=axesNames,
                         axesVals=axesVals, vals=allvals1, weights=allweights1)
        solset.makeSoltab('sol000', 'scalarphase', 'scalarphase000', axesNames=axesNames,
                         axesVals=axesVals, vals=allvals2, weights=allweights2)
    elif type == 'complexgain':
        solset.makeSoltab('sol000', 'phase', 'phase000', axesNames=axesNames,
                         axesVals=axesVals, vals=allvals1, weights=allweights1)
        solset.makeSoltab('sol000', 'amplitude', 'amplitude000', axesNames=axesNames,
                         axesVals=axesVals, vals=allvals2, weights=allweights2)


if __name__ == '__main__':
    descriptiontext = "Combine h5parms in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_files', help='list of input files')
    parser.add_argument('output_file', help='name of output file')
    parser.add_argument('soltype', help='solution type')
    parser.add_argument('solsetname', help='output solution set name')
    args = parser.parse_args()

    main(args.input_files, args.output_files, args.soltype, args.solsetname)
