"""
Definition of the master Operation class
"""
import os
import logging
import socket
import subprocess
import sys
import uuid
from factor import _logging
from jinja2 import Environment, FileSystemLoader
from lofarpipe.support.utilities import create_directory

DIR = os.path.dirname(os.path.abspath(__file__))
env_parset = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline', 'parsets')))
env_config = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline')))


class Operation(object):
    """
    Generic operation class

    An operation is simply a generic pipeline that performs a part of the
    processing. The corresponding operation object holds the pipeline settings,
    populates the pipeline config and parset templates, and updates the
    direction object with variables needed by later operations.

    Parameters
    ----------
    field : Field object
        Field for this operation
    direction : Direction object, optional
        Direction for this operation
    name : str, optional
        Name of the operation
    """
    def __init__(self, field, direction=None, name=None):
        self.parset = field.parset.copy()
        self.field = field
        self.name = name.lower()
        self.parset['op_name'] = name
        if direction is None:
            self.direction = field
        else:
            self.direction = direction
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor:{0}'.format(self.name))
        self.hostname = socket.gethostname()
        self.node_list = self.parset['cluster_specific']['node_list']

        # Working directory
        self.factor_working_dir = self.parset['dir_working']

        # Pipeline runtime and working dirs (pipeline makes subdir here with
        # name of direction)
        self.pipeline_runtime_dir = os.path.join(self.factor_working_dir,
                                                 'results', self.name)
        self.pipeline_working_dir = self.pipeline_runtime_dir
        create_directory(self.pipeline_runtime_dir)

        # Directory that holds the mapfiles
        self.pipeline_mapfile_dir = os.path.join(self.pipeline_runtime_dir,
                                                 self.direction.name, 'mapfiles')
        create_directory(self.pipeline_mapfile_dir)

        # Directory in the runtime dir that holds parset and config files (also
        # the results of the pipeline)
        self.pipeline_parset_dir = os.path.join(self.pipeline_runtime_dir,
                                                self.direction.name)
        create_directory(self.pipeline_parset_dir)

        # Directory that holds the mapfiles
        self.pipeline_mapfile_dir = os.path.join(self.pipeline_runtime_dir,
                                                 self.direction.name, 'mapfiles')
        create_directory(self.pipeline_mapfile_dir)

        # Local scratch directories and corresponding node recipes
        scratch_subdir = '{0}_{1}'.format(self.direction.name,
                                          str(uuid.uuid4().get_hex()[0:6]))
        if self.parset['cluster_specific']['dir_local'] is None:
            # Not specified
            self.local_scratch_dir = None
            self.local_dir_parent = None
            self.dppp_nodescript = 'executable_args'
        elif self.parset['cluster_specific']['cluster_type'].lower() == 'pbs':
            # PBS: use special DPPP node script
            self.local_scratch_dir = os.path.join(self.parset['cluster_specific']['dir_local'],
                                                  scratch_subdir)
            self.local_dir_parent = self.parset['cluster_specific']['dir_local']
            self.dppp_nodescript = 'dppp_scratch'
        elif self.parset['cluster_specific']['cluster_type'].lower() == 'slurm':
            # SLURM: use special DPPP node script
            self.local_scratch_dir = os.path.join(self.parset['cluster_specific']['dir_local'],
                                                  scratch_subdir)
            self.local_dir_parent = self.parset['cluster_specific']['dir_local']
            self.dppp_nodescript = 'dppp_scratch'
        else:
            # other: use given scratch directory and standard node script
            self.local_scratch_dir = os.path.join(self.parset['cluster_specific']['dir_local'],
                                                  scratch_subdir)
            self.local_dir_parent = self.parset['cluster_specific']['dir_local']
            self.dppp_nodescript = 'executable_args'

        # Directory that holds logs in a convenient place
        self.log_dir = os.path.join(self.factor_working_dir, 'logs', self.name)
        create_directory(self.log_dir)

        # Log name used for logs in log_dir
        self.logbasename = os.path.join(self.log_dir, self.direction.name)

        # Below are paths for scripts, etc. in the Factor install directory
        self.factor_root_dir = os.path.split(DIR)[0]
        self.factor_pipeline_dir = os.path.join(self.factor_root_dir, 'pipeline')
        self.factor_script_dir = os.path.join(self.factor_root_dir, 'scripts')
        self.factor_parset_dir = os.path.join(self.factor_root_dir, 'parsets')

        # Below are the templates and output paths for the pipeline parset and
        # config files. These may need to be re-defined in the subclasses
        # if the operation has non-standard template names
        self.pipeline_parset_template = '{0}_pipeline.parset'.format(self.name)
        self.pipeline_parset_file = os.path.join(self.pipeline_parset_dir,
                                                 'pipeline.parset')
        self.pipeline_config_template = 'pipeline.cfg'
        self.pipeline_config_file = os.path.join(self.pipeline_parset_dir,
                                                 'pipeline.cfg')

        # Define parameters needed for the pipeline config.
        self.cfg_dict = {'lofarroot': self.parset['cluster_specific']['lofarroot'],
                         'pythonpath': self.parset['cluster_specific']['lofarpythonpath'],
                         'factorroot': self.factor_root_dir,
                         'h5collector_executable': self.parset['cluster_specific']['h5collector_executable'],
                         'pipeline_working_dir': self.pipeline_working_dir,
                         'pipeline_runtime_dir': self.pipeline_runtime_dir,
                         'wsclean_executable': self.parset['cluster_specific']['wsclean_executable'],
                         'dppp_nodescript': self.dppp_nodescript,
                         }

        # Define global parameters needed by all pipeline parsets. Other,
        # pipeline-specific, parameters should be defined in the subclasses by
        # updating this dictionary
        self.parms_dict = {'parset_dir': self.factor_parset_dir,
                           'mapfile_dir': self.pipeline_mapfile_dir,
                           'pipeline_dir': self.factor_pipeline_dir,
                           'script_dir': self.factor_script_dir,
                           'local_dir': self.local_scratch_dir,
                           'local_dir_parent': self.local_dir_parent,
                           'pipeline_parset_dir': self.pipeline_parset_dir,
                           'hosts': self.node_list}

        # Add cluster-related info
        if self.parset['cluster_specific']['cluster_type'] == 'localhost':
            self.cfg_dict['remote'] = '[remote]\n' + \
                                      'method = local\n' + \
                                      'max_per_node = {0}\n'.format(self.parset['cluster_specific']['ncpu'])
        elif self.parset['cluster_specific']['cluster_type'] == 'juropa_slurm':
            self.cfg_dict['remote'] = '[remote]\n' + \
                                      'method = slurm_srun\n' + \
                                      'max_per_node = {0}\n'.format(self.parset['cluster_specific']['ncpu'])
        elif self.parset['cluster_specific']['cluster_type'] == 'mpirun':
            self.cfg_dict['remote'] = '[remote]\n' + \
                                      'method = mpirun\n' + \
                                      'max_per_node = {0}\n'.format(self.parset['cluster_specific']['ncpu'])
        elif (self.parset['cluster_specific']['cluster_type'] == 'pbs' or
              self.parset['cluster_specific']['cluster_type'] == 'slurm'):
            self.cfg_dict['remote'] = ''
        else:
            self.log.error('Could not determine the nature of your cluster!')
            sys.exit(1)

    def update_dicts(self):
        """
        Update the dicts used for the pipeline parset templates
        """
        self.cfg_dict.update(self.direction.__dict__)
        self.parms_dict.update(self.direction.__dict__)

    def setup(self):
        """
        Set up this operation

        This involves just filling the pipeline config and parset templates.
        Generally, this does not need to be re-defined in the subclasses
        unless the operation has non-standard template names
        """
        # Update the dictionaries with the attributes of the operation's
        # direction object. Any attributes set in the direction object that are
        # also in the parms_dict will be set to those of the direction object
        # (e.g., 'max_proc_per_node', which is set in the direction object by
        # factor.cluster.divide_nodes() will override the value set above)
        self.update_dicts()

        self.pipeline_parset_template = env_parset.get_template(self.pipeline_parset_template)
        tmp = self.pipeline_parset_template.render(self.parms_dict)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        self.pipeline_config_template = env_config.get_template(self.pipeline_config_template)
        tmp = self.pipeline_config_template.render(self.cfg_dict)
        with open(self.pipeline_config_file, 'w') as f:
            f.write(tmp)

    def finalize(self):
        """
        Finalize this operation

        This should be defined in the subclasses as needed
        """
        pass

    def check_existing_files(self, mapfile):
        """
        Checks if files in input mapfile exist

        Parameters
        ----------
        mapfile : str
            Filename of mapfile to check

        Returns
        -------
        all_exist : bool
            True if all files in mapfile exist, False if not
        """
        from lofarpipe.support.data_map import DataMap

        all_exist = True
        self.log.debug('Checking for existing files...')
        try:
            datamap = DataMap.load(mapfile)
            for item in datamap:
                # Handle case in which item.file is a Python list
                if item.file[0] == '[' and item.file[-1] == ']':
                    files = item.file.strip('[]').split(',')
                else:
                    files = [item.file]
                for f in files:
                    if not os.path.exists(f):
                        all_exist = False
            if all_exist:
                self.log.debug('...all files exist')
            else:
                self.log.debug('...one or more files not found')
            return all_exist
        except IOError:
            self.log.debug('Could not read mapfile {}. Skipping it'.format(mapfile))
            return False

    def can_restart(self):
        """
        Checks the pipeline log for certain conditions that affect auto restarting

        Returns
        -------
        can_restart : bool
            True if pipeline log indicates an error for which auto restart is
            possible
        """
        logfile = self.logbasename + '.out.log'
        can_restart = False
        if os.path.exists(logfile):
            # Read the last few lines and look for 'returncode 123456'
            try:
                with open(logfile, "rb") as f:
                    first = f.readline()        # Read the first line.
                    f.seek(-10000, 2)           # Jump back from end
                    while f.read(1) != b"\n":   # Until EOL is found...
                        f.seek(-2, 1)           # ...jump back the read byte plus one more.
                    last_lines = f.readlines()  # Read last line.

                for line in last_lines:
                    if 'returncode 123456' in line:
                        can_restart = True
                        break
            except IOError:
                can_restart = False

        return can_restart

    def get_steptypes(self):
        """
        Returns the step types of completed pipeline steps

        Returns
        -------
        steptypes : list
            List of step types
        """
        import pickle

        statefile = os.path.join(self.pipeline_parset_dir, 'statefile')
        if os.path.exists(statefile):
            current_state = pickle.load(open(statefile, 'rb'))
            steptypes = [item[0] for item in current_state[1]]
        else:
            steptypes = []

        return steptypes

    def reset_state_to_steptype(self, steptype):
        """
        Resets the pipeline state to before the given steptype

        Steptype is the type of the step as defined in the parset under
        step.control.type

        Parameters
        ----------
        steptype : str
            Step type from which to alter state
        """
        import pickle

        statefile = os.path.join(self.pipeline_parset_dir, 'statefile')
        current_state = pickle.load(open(statefile, 'rb'))
        steptypes = [item[0] for item in current_state[1]]

        # delete steps from the first instance of the given step type
        del_number = steptypes.index(steptype)
        current_state[1] = current_state[1][:del_number]

        pickle.dump(current_state, open(statefile, 'wb'))

    def cleanup(self):
        """
        Cleans up temp files in the scratch directories of each node
        """
        if self.local_scratch_dir is not None:
            for node in self.node_list:
                if node == 'localhost':
                    cmd = ['rm', '-rf', self.local_scratch_dir]
                else:
                    cmd = ['ssh', node, 'rm', '-rf', self.local_scratch_dir]
                tmp = subprocess.call(cmd)

            # Check whether we need to reset the pipeline state to before the sync step
            steptypes = self.get_steptypes()
            if 'sync_files' in steptypes and 'remove_synced_data' not in steptypes:
                self.reset_state_to_steptype('sync_files')
