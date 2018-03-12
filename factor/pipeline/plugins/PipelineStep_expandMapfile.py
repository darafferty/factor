import os
from lofarpipe.support.data_map import DataMap, DataProduct


def plugin_main(args, **kwargs):
    """
    Makes a mapfile by duplicating the entries in an input mapfile

    Note that when mapfile_in has multiple entries, we duplicate each entry N times in
    the output mapfile, where N = len(map_match) / len(map_in), such that the output is:

    [in0, in0, ..., in0, in1, in1, ..., in1, in2, ...]

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap containing single item
    mapfile_to_match : str
        Filename of datamap containing multiple items
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        New parmdb datamap filename
    """
    mapfile_in = kwargs['mapfile_in']
    mapfile_to_match = kwargs['mapfile_to_match']
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_in = DataMap.load(mapfile_in)
    map_match = DataMap.load(mapfile_to_match)
    stride = len(map_match) / len(map_in)

    map_out = DataMap([])
    if 'suffix_to_add' in kwargs:
        suffix_to_add = kwargs['suffix_to_add']
        os.system('cp {0} {0}{1}'.format(map_in[0].file, suffix_to_add))
    else:
        suffix_to_add = ''

    map_match.iterator = DataMap.SkipIterator
    for i, item_in in enumerate(map_in):
        for item in map_match[i:i+stride]:
            map_out.data.append(DataProduct(item.host, '{0}{1}'.format(item_in.file, suffix_to_add), item.skip))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
