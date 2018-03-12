import os
from shutil import copyfile
from lofarpipe.support.data_map import DataMap, DataProduct


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
        rename =  string2bool(kwargs['rename'])
    else:
        rename = False
    if 'append_index' in kwargs:
        append_index = string2bool(kwargs['append_index'])
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
        # Find stride for indexing
        stride = len(map_in) / len(set([f for item.file for item in map_in]))

    j = -1
    for i, item in enumerate(map_in):
        if append_index:
            if i % stride == 0:
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


def string2bool(instring):
    if not isinstance(instring, basestring):
        raise ValueError('string2bool: Input is not a basic string!')
    if instring.upper() == 'TRUE' or instring == '1':
        return True
    elif instring.upper() == 'FALSE' or instring == '0':
        return False
    else:
        raise ValueError('string2bool: Cannot convert string "'+instring+'" to boolean!')
