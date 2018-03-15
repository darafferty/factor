#!/usr/bin/python
"""
Script to plot solutions
"""
from losoto.h5parm import h5parm
from losoto.operations import plot
import sys, os
import argparse
from argparse import RawTextHelpFormatter


def main(h5file, solset, soltab, root):
    """
    Make various plots
    """
    h = h5parm(h5file)
    ss = h.getSolset(solset)
    st = ss.getSoltab(soltab)

    if st.getType() == 'tec':
        ncol = 1
    else:
        ncol = 0
    if st.getType() == 'amplitude' or st.getType() == 'phase':
        color = 'pol'
    else:
        color = ''

    plot.run(st, ['time'], axisInTable='ant', axisInCol=color, NColFig=ncol)


if __name__ == "__main__":
    descriptiontext = "Plot solutions.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h5file', help="Name of solution h5parm file")
    parser.add_argument('solset', help="Name of solution set")
    parser.add_argument('soltab', help="Name of solution table")
    parser.add_argument('root', help="Root name for output plots")

    args = parser.parse_args()
    main(args.h5file, args.solset, args.soltab, args.root)

