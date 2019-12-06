"""
Definition of the master Operation class
"""
import os
import logging
from factor import _logging
from jinja2 import Environment, FileSystemLoader
from factor.lib import miscellaneous as misc

DIR = os.path.dirname(os.path.abspath(__file__))
env_parset = Environment(loader=FileSystemLoader(os.path.join(DIR, '..', 'pipeline', 'parsets')))


class Operation(object):
    """
    Generic operation class

    An operation is simply a CWL pipeline that performs a part of the
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
    def __init__(self, field, direction=None, name=None, index=None):
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
        if direction is None:
            self.direction = field
        else:
            self.direction = direction
        _logging.set_level(self.parset['logging_level'])
        self.log = logging.getLogger('factor:{0}'.format(self.name))

        # Factor working directory
        self.factor_working_dir = self.parset['dir_working']

        # Pipeline working dir
        self.pipeline_working_dir = os.path.join(self.factor_working_dir,
                                                 'pipelines', self.name,
                                                 self.direction.name)
        misc.create_directory(self.pipeline_working_dir)

        # Directory that holds the pipeline logs in a convenient place
        self.log_dir = os.path.join(self.factor_working_dir, 'logs', self.name)
        misc.create_directory(self.log_dir)

        # Log name used for logs in log_dir
        self.logbasename = os.path.join(self.log_dir, self.direction.name)

        # Paths for scripts, etc. in the Factor install directory
        self.factor_root_dir = os.path.split(DIR)[0]
        self.factor_pipeline_dir = os.path.join(self.factor_root_dir, 'pipeline')
        self.factor_script_dir = os.path.join(self.factor_root_dir, 'scripts')

        # Input template name and output parset and inputs filenames for
        # the pipeline
        self.pipeline_parset_template = '{0}_pipeline.cwl'.format(self.rootname)
        self.pipeline_parset_file = os.path.join(self.pipeline_working_dir,
                                                 'pipeline_parset.cwl')
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

        This involves filling the pipeline parset template and writing the inputs file.
        Generally, this does not need to be re-defined in the subclasses unless the
        operation has non-standard template name
        """
        # Fill the parset template and save to a file
        self.set_parset_parameters()
        self.pipeline_parset_template = env_parset.get_template(self.pipeline_parset_template)
        tmp = self.pipeline_parset_template.render(self.parset_parms)
        with open(self.pipeline_parset_file, 'w') as f:
            f.write(tmp)

        # Save the pipeline inputs to a file
        self.set_input_parameters()
        keys = []
        vals = []
        for k, v in self.input_parms.items():
            keys.append(k)
            if type(v) is bool:
                vals.append("'{}'".format(v))
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
