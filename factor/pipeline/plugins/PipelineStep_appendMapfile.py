import os
from shutil import copyfile
from lofarpipe.support.data_map import DataMap, DataProduct
from factor.lib import miscellaneous as misc


def plugin_main(args, **kwargs):
    """
    Appends a string to filenames in a mapfile

    Parameters
    ----------
    mapfile_in : str
        Filename of datamap to append to
    append : str
        String to append
    append_index : bool
        If True, append a unique index to each file
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile
    rename: bool
        If True, rename the original filenames to the new ones

    Returns
    -------
    result : dict
        New datamap filename

    """
    mapfile_in = kwargs['mapfile_in']

    if 'rename' in kwargs:
        rename =  misc.string2bool(kwargs['rename'])
    else:
        rename = False
    if 'append_index' in kwargs:
        append_index = misc.string2bool(kwargs['append_index'])
    else:
        append_index = False

    append_str = kwargs['append']
    if append_str == 'None':
        append_str = ''
    mapfile_dir = kwargs['mapfile_dir']
    filename = kwargs['filename']

    map_out = DataMap([])
    map_in = DataMap.load(mapfile_in)
    if append_index:
        # Find number of unique files in input mapfile
        nfiles = len(set([item.file for item in map_in]))

    j = -1
    for i, item in enumerate(map_in):
        if append_index:
            if i % nfiles == 0:
                j += 1
            newfile = item.file+append_str+'_{}'.format(j)
        else:
            newfile = item.file+append_str
        map_out.data.append(DataProduct(item.host, newfile, item.skip))
        if rename:
            if os.path.exists(newfile):
                os.remove(newfile)
            copyfile(item.file, newfile)

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
