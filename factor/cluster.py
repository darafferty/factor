"""
Module that holds all compute-cluster-related functions
"""
import os
import logging
import sys
import re
import numpy as np

log = logging.getLogger('factor:cluster')


def make_pbs_clusterdesc():
    """
    Make a cluster description file from the PBS_NODEFILE

    Returns
    -------
    clusterdesc_file
        Filename of resulting cluster description file
    """
    nodes = []
    try:
        filename = os.environ['PBS_NODEFILE']
    except KeyError:
        log.error('PBS_NODEFILE not found. You must have a reservation to '
                  'use clusterdesc = PBS.')
        sys.exit(1)

    with open(filename, 'r') as file:
        for line in file:
            node_name = line.split()[0]
            if node_name not in nodes:
                nodes.append(node_name)

    lines = ['# Clusterdesc file to do parallel processing with PBS / torque\n\n']
    lines.append('ClusterName = PBS\n\n')
    lines.append('# Compute nodes\n')
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(sorted(nodes))))

    clusterdesc_file = 'factor_pbs.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)
    log.info('Using {0} node(s)'.format(len(nodes)))

    return clusterdesc_file


def expand_part(s):
    """Expand a part (e.g. "x[1-2]y[1-3][1-3]") (no outer level commas).

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """
    # Base case: the empty part expand to the singleton list of ""
    if s == "":
        return [""]

    # Split into:
    # 1) prefix string (may be empty)
    # 2) rangelist in brackets (may be missing)
    # 3) the rest

    m = re.match(r'([^,\[]*)(\[[^\]]*\])?(.*)', s)
    (prefix, rangelist, rest) = m.group(1, 2, 3)

    # Expand the rest first (here is where we recurse!)
    rest_expanded = expand_part(rest)

    # Expand our own part
    if not rangelist:
        # If there is no rangelist, our own contribution is the prefix only
        us_expanded = [prefix]
    else:
        # Otherwise expand the rangelist (adding the prefix before)
        us_expanded = expand_rangelist(prefix, rangelist[1:-1])

    return [us_part + rest_part
            for us_part in us_expanded
            for rest_part in rest_expanded]


def expand_range(prefix, range_):
    """ Expand a range (e.g. 1-10 or 14), putting a prefix before.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """
    # Check for a single number first
    m = re.match(r'^[0-9]+$', range_)
    if m:
        return ["%s%s" % (prefix, range_)]

    # Otherwise split low-high
    m = re.match(r'^([0-9]+)-([0-9]+)$', range_)

    (s_low, s_high) = m.group(1, 2)
    low = int(s_low)
    high = int(s_high)
    width = len(s_low)

    results = []
    for i in xrange(low, high+1):
        results.append("%s%0*d" % (prefix, width, i))
    return results


def expand_rangelist(prefix, rangelist):
    """ Expand a rangelist (e.g. "1-10,14"), putting a prefix before.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """
    # Split at commas and expand each range separately
    results = []
    for range_ in rangelist.split(","):
        results.extend(expand_range(prefix, range_))
    return results


def expand_hostlist(hostlist, allow_duplicates=False, sort=False):
    """Expand a hostlist expression string to a Python list.

    Example: expand_hostlist("n[9-11],d[01-02]") ==>
             ['n9', 'n10', 'n11', 'd01', 'd02']

    Unless allow_duplicates is true, duplicates will be purged
    from the results. If sort is true, the output will be sorted.

    Note: Adapted from git://www.nsc.liu.se/~kent/python-hostlist.git
    """
    results = []
    bracket_level = 0
    part = ""

    for c in hostlist + ",":
        if c == "," and bracket_level == 0:
            # Comma at top level, split!
            if part:
                results.extend(expand_part(part))
            part = ""
        else:
            part += c

        if c == "[":
            bracket_level += 1
        elif c == "]":
            bracket_level -= 1

    seen = set()
    results_nodup = []
    for e in results:
        if e not in seen:
            results_nodup.append(e)
            seen.add(e)
    return results_nodup


def make_slurm_clusterdesc():
    """
    Make a cluster description file from the SLURM_JOB_NODELIST

    Returns
    -------
    clusterdesc_file
        Filename of resulting cluster description file
    """
    nodes = []
    try:
        hostlist = os.environ['SLURM_JOB_NODELIST']
    except KeyError:
        log.error('SLURM_JOB_NODELIST not found. You must have a reservation to '
                  'use clusterdesc = SLURM.')
        sys.exit(1)

    nodes = expand_hostlist(hostlist)

    lines = ['# Clusterdesc file to do parallel processing with SLURM\n\n']
    lines.append('ClusterName = SLURM\n\n')
    lines.append('# Compute nodes\n')
    lines.append('Compute.Nodes = [{0}]\n'.format(', '.join(sorted(nodes))))

    clusterdesc_file = 'factor_slurm.clusterdesc'
    with open(clusterdesc_file, 'wb') as file:
        file.writelines(lines)
    log.info('Using {0} node(s)'.format(len(nodes)))

    return clusterdesc_file


def get_compute_nodes(clusterdesc_file):
    """
    Read a cluster description file and return list of nodes

    Parameters
    ----------
    clusterdesc_file : str
        Filename of cluster description file

    Returns
    -------
    result : list
        Sorted list of node names
    """
    from lofarpipe.support import clusterdesc

    cluster = clusterdesc.ClusterDesc(clusterdesc_file)
    return sorted(clusterdesc.get_compute_nodes(cluster))


def find_executables(cluster_parset):
    """
    Adds the paths to required executables to parset dict

    Parameters
    ----------
    cluster_parset : dict
        Cluster-specific parset dictionary
    """
    from distutils import spawn

    executables = {'genericpipeline_executable': ['genericpipeline.py'],
                   'wsclean_executable': ['wsclean'],
                   'losoto_executable': ['losoto'],
                   'h5collector_executable': ['H5parm_collector.py']}
    for key, names in executables.iteritems():
        for name in names:
            path = spawn.find_executable(name)
            if path is not None:
                cluster_parset[key] = path
                break
        if path is None:
            log.error('The path to the {0} executable could not be determined. '
                      'Please make sure it is in your PATH.'.format(name))
            sys.exit(1)


def get_type(cluster_parset):
    """
    Gets the cluster type and sets relevant entries in the input parset

    Parameters
    ----------
    cluster_parset : dict
        Cluster-specific parset dictionary
    """
    if cluster_parset['cluster_type'].lower() == 'pbs':
        log.info('Using cluster setting: "PBS".')
        cluster_parset['clusterdesc'] = make_pbs_clusterdesc()
        cluster_parset['clustertype'] = 'pbs'
    elif cluster_parset['cluster_type'].lower() == 'slurm':
        log.info('Using cluster setting: "SLURM".')
        cluster_parset['clusterdesc'] = make_slurm_clusterdesc()
        cluster_parset['clustertype'] = 'slurm'
    elif cluster_parset['cluster_type'].lower() == 'juropa_slurm':
        log.info('Using cluster setting: "JUROPA_slurm" (Single '
                 'genericpipeline using multiple nodes).')
        # slurm_srun on JUROPA uses the local.clusterdesc
        cluster_parset['clusterdesc'] = os.path.join(cluster_parset['lofarroot'],
                                                     'share', 'local.clusterdesc')
        cluster_parset['clustertype'] = 'juropa_slurm'
        cluster_parset['node_list'] = ['localhost']
    elif cluster_parset['cluster_type'].lower() == 'mpirun':
        log.info('Using cluster setting: "mpirun".')
        # mpirun uses the local.clusterdesc?
        cluster_parset['clusterdesc'] = os.path.join(cluster_parset['lofarroot'],
                                                     'share', 'local.clusterdesc')
        cluster_parset['clustertype'] = 'mpirun'
        cluster_parset['node_list'] = ['localhost']
    else:
        log.info('Using cluster setting: "local" (Single node).')
        cluster_parset['clusterdesc'] = cluster_parset['lofarroot'] + '/share/local.clusterdesc'
        cluster_parset['clustertype'] = 'local'
    if 'node_list' not in cluster_parset:
        cluster_parset['node_list'] = get_compute_nodes(cluster_parset['clusterdesc'])


def check_ulimit(cluster_parset):
    """
    Checks the limit on number of open files

    Parameters
    ----------
    cluster_parset : dict
        Cluster-specific parset dictionary
    """
    try:
        import resource
        nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        if cluster_parset['clustertype'] == 'local' and nof_files_limits[0] < nof_files_limits[1]:
            log.debug('Setting limit for number of open files to: {}.'.format(nof_files_limits[1]))
            resource.setrlimit(resource.RLIMIT_NOFILE, (nof_files_limits[1], nof_files_limits[1]))
            nof_files_limits = resource.getrlimit(resource.RLIMIT_NOFILE)
        log.debug('Active limit for number of open files is {0}, maximum limit '
                  'is {1}.'.format(nof_files_limits[0], nof_files_limits[1]))
        if nof_files_limits[0] < 2048:
            log.warn('The limit for number of open files is small, this could '
                     'result in a "Too many open files" problem when running factor.')
            log.warn('The active limit can be increased to the maximum for the '
                     'user with: "ulimit -Sn <number>" (bash) or "limit descriptors 1024" (csh).')
    except resource.error:
        log.warn('Cannot check limits for number of open files, what kind of system is this?')


def get_total_memory():
    """
    Returns the total memory in GB
    """
    tot_gb, used_gb, free_gb = map(int, os.popen('free -t -g').readlines()[-1].split()[1:])

    return tot_gb


def get_time_chunksize(cluster_parset, timepersample, numsamples, solint_fast_timestep):
    """
    Returns the target chunk size in seconds for an observation

    Parameters
    ----------
    cluster_parset : dict
        Cluster-specific parset dictionary
    timepersample : float
        Time in seconds per time sample
    numsamples : int
        Total number of time samples in the observation
    solint_fast_timestep : int
        Number of time samples in fast-phase solve
    """
    # Try to make at least as many time chunks as there are nodes
    n_nodes = len(cluster_parset['node_list'])
    samplesperchunk = np.ceil(numsamples / solint_fast_timestep / n_nodes)
    target_time_chunksize = timepersample * samplesperchunk

    return target_time_chunksize


def get_frequency_chunksize(cluster_parset, channelwidth, solint_slow_freqstep,
                            solint_slow_timestep, antenna):
    """
    Returns the target chunk size in seconds for an observation

    Parameters
    ----------
    cluster_parset : dict
        Cluster-specific parset dictionary
    channelwidth : float
        Bandwidth in Hz per frequency sample
    solint_slow_freqstep : int
        Number of frequency samples in slow-gain solve
    solint_slow_timestep : int
        Number of time samples in slow-gain solve
    antenna : str
        Antenna type: "HBA" or "LBA"
    """
    # Try to make at least as many time chunks as there are nodes
    n_cpus = cluster_parset['ncpu']
    mem_gb = cluster_parset['fmem'] * get_total_memory()
    if antenna == 'HBA':
        # Memory usage in GB/chan/timeslot/dir of a typical HBA observation
        mem_usage_gb = 1e-3
        ndir = 40
    elif antenna == 'LBA':
        # Memory usage in GB/chan/timeslot/dir of a typical LBA observation
        mem_usage_gb = 2.5e-4
        ndir = 20
    gb_per_solint = mem_usage_gb * solint_slow_freqstep * solint_slow_timestep * ndir
    nsolints = int(round(mem_gb / gb_per_solint))
    channelsperchunk = np.ceil(solint_slow_freqstep * nsolints)
    target_freq_chunksize = channelwidth * channelsperchunk

    return target_freq_chunksize
