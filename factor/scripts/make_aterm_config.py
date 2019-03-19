#! /usr/bin/env python
"""
Script to make an aterm-config file for WSClean
"""


def main(output_file, tec_image=None, gain_image=None, use_beam=False):
    """
    Make an aterm-config file for WSClean

    Parameters
    ----------
    output_file : str
        Filename of output config file
    tec_image : str, optional
        Full path to TEC image
    gain_image : str, optional
        Full path to gain image
    use_beam : bool, optional
        If True, use the beam with IDG
    """
    if tec_image is None and gain_image is None:
        print('make_aterm_config: One of tec_image or gain_image must be specified')
        sys.exit(1)

    terms = []
    if tec_image is not None:
        terms.sppend('tec')
    if gain_image is not None:
        terms.sppend('gain')
    if use_beam is not None:
        terms.sppend('beam')
    aterm_str = 'aterms = [{}]'.format(', '.join(terms))
    lines = [aterm_str]
    if tec_image is not None:
        lines.append('tec.images = [{}]'.format(tec_image))
    if gain_image is not None:
        lines.append('gain.images = [{}]'.format(gain_image))
    if use_beam:
        lines.append('beam.differential = true')
        lines.append('beam.update_interval = 120')
        lines.append('beam.usechannelfreq = true')

    config_file = open(output_file, 'w')
    config_file.writelines(lines)
    config_file.close()
