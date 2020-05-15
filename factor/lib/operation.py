"""
Definition of the master Operation class
"""
import os
import sys
import logging
from factor import _logging
from jinja2 import Environment, FileSystemLoader
from factor.lib import miscellaneous as misc
from toil.leader import FailedJobsException
from toil.cwl import cwltoil
from factor.lib.context import Timer

DIR = os.path.dirname(os.path.abspath(__file__))
env_parset = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline', 'parsets')))


class Operation(object):
    """
    Generic operation class

    An operation is simply a CWL pipeline that performs a part of the
    processing. It holds the pipeline settings, populates the pipeline input and
    parset templates, and runs the pipeline. The field object is passed between
    operations, each of which updates it with variables needed by other, subsequent,
    operations.

    Parameters
    ----------
    field : Field object
        Field for this operation
    name : str, optional
        Name of the operation
    index : int, optional
        Index of the operation
    """
    def __init__(self, field, name=None, index=None):
        self.parset = field.parset.copy()
        self.field = field
        self.rootname = name.lower()
        self.index = index
        if self.index is not None:
            self.name = '{0}_{1}'.format(self.rootname, self.index)
        else:
            self.name = self.rootname
        self.rootname
        self.parset['op_name'] = name
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor:{0}'.format(self.name))

        # Factor working directory
        self.factor_working_dir = self.parset['dir_working']

        # Pipeline working dir
        self.pipeline_working_dir = os.path.join(self.factor_working_dir,
                                                 'pipelines', self.name)
        misc.create_directory(self.pipeline_working_dir)

        # Temp directory (local to the nodes)
        self.temp_dir = self.parset['cluster_specific']['dir_local']

        # Directory that holds the pipeline logs in a convenient place
        self.log_dir = os.path.join(self.factor_working_dir, 'logs')
        misc.create_directory(self.log_dir)

        # Log name used for logs in log_dir
        self.logbasename = os.path.join(self.log_dir, self.name)

        # Paths for scripts, etc. in the Factor install directory
        self.factor_root_dir = os.path.split(DIR)[0]
        self.factor_pipeline_dir = os.path.join(self.factor_root_dir, 'pipeline')
        self.factor_script_dir = os.path.join(self.factor_root_dir, 'scripts')

        # Input template name and output parset and inputs filenames for
        # the pipeline. If the pipeline uses a subworkflow, its template filename must be
        # defined in the subclass by self.subpipeline_parset_template to the right
        # path
        self.pipeline_parset_template = '{0}_pipeline.cwl'.format(self.rootname)
        self.subpipeline_parset_template = None
        self.pipeline_parset_file = os.path.join(self.pipeline_working_dir,
                                                 'pipeline_parset.cwl')
        self.subpipeline_parset_file = os.path.join(self.pipeline_working_dir,
                                                    'subpipeline_parset.cwl')
        self.pipeline_inputs_file = os.path.join(self.pipeline_working_dir,
                                                 'pipeline_inputs.yml')

    def set_parset_parameters(self):
        """
        Define parameters needed for the pipeline parset template

        The dictionary keys must match the jinja template variables used in the
        corresponding pipeline parset.

        The entries are defined in the subclasses as needed
        """
        self.parset_parms = {}

    def set_input_parameters(self):
        """
        Define parameters needed for the pipeline inputs

        The dictionary keys must match the workflow inputs defined in the corresponding
        pipeline parset.

        The entries are defined in the subclasses as needed
        """
        self.input_parms = {}

    def setup(self):
        """
        Set up this operation

        This involves filling the pipeline parset template and writing the inputs file
        """
        # Fill the parset template and save to a file
        self.set_parset_parameters()
        self.pipeline_parset_template = env_parset.get_template(self.pipeline_parset_template)
        tmp = self.pipeline_parset_template.render(self.parset_parms)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)
        if self.subpipeline_parset_template is not None:
            self.subpipeline_parset_template = env_parset.get_template(self.subpipeline_parset_template)
            tmp = self.subpipeline_parset_template.render(self.parset_parms)
            with open(self.subpipeline_parset_file, 'w') as f:
                f.write(tmp)

        # Save the pipeline inputs to a file
        self.set_input_parameters()
        keys = []
        vals = []
        for k, v in self.input_parms.items():
            keys.append(k)
            if type(v) is bool:
                vals.append("'{}'".format(v))
            elif type(v) is list and type(v[0]) is bool:
                vals.append('[{}]'.format(','.join(["'{}'".format(ve) for ve in v])))
            else:
                vals.append(v)
        tmp = '\n'.join(['{0}: {1}'.format(k, v) for k, v in zip(keys, vals)])
        with open(self.pipeline_inputs_file, 'w') as f:
            f.write(tmp)

    def finalize(self):
        """
        Finalize this operation

        This should be defined in the subclasses as needed
        """
        pass

    def call_toil(self):
        """
        Calls Toil to run the operation's pipeline
        """
        max_nodes = self.parset['cluster_specific']['max_nodes']
        scratch_dir = self.parset['cluster_specific']['dir_local']
        batch_system = self.parset['cluster_specific']['batch_system']
        jobstore = os.path.join(self.pipeline_working_dir, 'jobstore')

        # Build the args list
        args = []
        args.extend(['--batchSystem', batch_system])
        if batch_system == 'slurm':
            args.extend(['--disableCaching'])
            args.extend(['--defaultCores', '6'])
            args.extend(['--defaultMemory', '1M'])
            args.extend(['--maxLocalJobs', str(max_nodes)])
        if batch_system == 'singleMachine':
            args.extend(['--maxLocalJobs', str(1)])
        args.extend(['--jobStore', jobstore])
        if os.path.exists(jobstore):
            args.extend(['--restart'])
        args.extend(['--basedir', self.pipeline_working_dir])
        args.extend(['--workDir', self.pipeline_working_dir])
        args.extend(['--outdir', scratch_dir])
        args.extend(['--logFile', self.logbasename+'.log'])
        args.extend(['--preserve-entire-environment'])
        if scratch_dir is not None:
            args.extend(['--tmpdir-prefix', scratch_dir])
            args.extend(['--tmp-outdir-prefix', scratch_dir])
        args.extend(['--logLevel', 'DEBUG'])
        args.extend(['--clean', 'never'])
        args.extend(['--servicePollingInterval', '10'])
        args.append(self.pipeline_parset_file)
        args.append(self.pipeline_inputs_file)

        # Run the pipeline
        self.log.info('<-- Operation {0} started'.format(self.name))
        try:
            status = cwltoil.main(args=args)
            if status == 0:
                self.success = True
            else:
                self.success = False
        except FailedJobsException:
            self.success = False

    def run(self):
        """
        Runs the operation
        """
        self.setup()
        with Timer(self.log):
            self.call_toil()
        if self.success:
            self.log.info('--> Operation {0} completed'.format(self.name))
            self.finalize()
        else:
            self.log.error('Operation {0} failed due to an error'.format(self.name))
            sys.exit(1)
