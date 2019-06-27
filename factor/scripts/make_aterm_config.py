#! /usr/bin/env python
"""
Script to make an aterm-config file for WSClean
"""
import sys
from lofarpipe.support.data_map import DataMap


def main(input_file, output_file, tec_mapfile=None, gain_mapfile=None, use_beam=False):
    """
    Make an aterm-config file for WSClean

    Parameters
    ----------
    input_file : str
        Filename of input file (not used but needed for pipeline)
    output_file : str
        Filename of output config file
    tec_mapfile : str, optional
        Full path to mapfile with TEC images
    gain_mapfile : str, optional
        Full path to mapfile with gain images
    use_beam : bool, optional
        If True, use the beam with IDG
    """
    if tec_mapfile is None and gain_mapfile is None:
        print('make_aterm_config: One of tec_mapfile or gain_mapfile must be specified')
        sys.exit(1)

    terms = []
    if tec_mapfile is not None:
        terms.append('tec')
        tec_map = DataMap.load(tec_mapfile)
        tec_images = [item.file for item in tec_map]
    if gain_mapfile is not None:
        terms.append('diagonal')
        gain_map = DataMap.load(gain_mapfile)
        gain_images = [item.file for item in gain_map]
    if use_beam:
        terms.append('beam')
    aterm_str = 'aterms = [{}]\n'.format(', '.join(terms))
    lines = [aterm_str]
    if tec_mapfile is not None:
        lines.append('tec.images = [{}]\n'.format(' '.join(tec_images)))
    if gain_mapfile is not None:
        lines.append('diagonal.images = [{}]\n'.format(' '.join(gain_images)))
    if use_beam:
        lines.append('beam.differential = true\n')
        lines.append('beam.update_interval = 120\n')
        lines.append('beam.usechannelfreq = true\n')

    config_file = open(output_file, 'w')
    config_file.writelines(lines)
    config_file.close()
