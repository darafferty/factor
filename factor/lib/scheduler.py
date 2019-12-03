"""
Module defining the Scheduler class
"""
import logging
import multiprocessing
import os
import sys
from factor.lib.context import Timer
import factor.cluster

log = logging.getLogger('factor:scheduler')


def call_toil(op_name, direction_name, parset, inputs, basedir, dir_local, logbasename,
              batch_system):
    """
    Calls Toil to run the pipeline

    Parameters
    ----------
    op_name : str
        Name of operation
    direction_name : str
        Name of direction
    parset : str
        Name of pipeline parset file
    inputs : str
        Name of pipeline inputs file
    basedir : str
        Path of pipeline working directory
    dir_local : str
        Path of local scratch directory
    logbasename : str
        Log file base name
    batch_system : str
        Type of batch system to use
    """
    from toil.cwl import cwltoil
    from factor.lib.context import RedirectStdStreams

    # Build the args dict
    args = {}
    args['cwltool'] = parset
    args['cwljob'] = inputs
    args['batchSystem'] = batch_system
    if batch_system == 'slurm':
        args['disableCaching'] = True
        args['defaultCores'] = 6
        args['defaultMemory'] = '1M'
    args['jobStore'] = os.path.join(basedir, 'jobstore')
    args['basedir'] = os.path.join(basedir)
    args['outdir'] = os.path.join(basedir)
    args['preserve-environment'] = True
    args['tmpdir-prefix'] = dir_local
    args['tmp-outdir-prefix'] = os.path.join(basedir)

    # Run the pipeline, redirecting screen output to log files.
    if direction_name == 'field':
        log.info('<-- Operation {0} started'.format(op_name))
    else:
        log.info('<-- Operation {0} started (direction: {1})'.format(op_name, direction_name))
    with open("{0}.log".format(logbasename), "w") as out:
        with RedirectStdStreams(stdout=out, stderr=out):
            status = cwltoil.main(args)

    return (op_name, direction_name, status)


class Scheduler(object):
    """
    The scheduler runs all jobs sent to it in parallel

    Parameters
    ----------
    parset : str
        Factor parset
    name : str, optional
        Name of the scheduler
    """
    def __init__(self, parset, name='scheduler'):

        self.parset = parset['cluster_specific'].copy()
        factor.cluster.check_ulimit(self.parset)
        self.nops_simul = self.parset['max_nodes']
        self.scratch_dir = self.parset['dir_local']
        self.batch_system = self.parset['batch_system']
        self.name = name
        self.success = True

    def result_callback(self, result):
        """
        Callback function for apply_async result
        """
        op_name, direction_name, status = result

        # Identify the current operation from the direction name
        try:
            this_op_indx = [op.direction.name for op in self.operation_list].index(direction_name)
            this_op = self.operation_list[this_op_indx]
        except ValueError:
            log.warn('Operation {0} (direction: {1}) not in list of active '
                     'operations. This could indicate a problem with the operation'.
                     format(op_name, direction_name))
            return

        # Finalize the operation
        if status == 0:
            if direction_name == 'field':
                log.info('--> Operation {0} completed'.format(op_name))
            else:
                log.info('--> Operation {0} completed (direction: '
                         '{1})'.format(op_name, direction_name))
            this_op.finalize()
        else:
            self.success = False
            if this_op.can_restart():
                log.warning('Operation {0} failed due to error (direction: '
                            '{1}) but will be automatically resumed'.format(op_name, direction_name))
            else:
                log.error('Operation {0} failed due to an error (direction: '
                          '{1})'.format(op_name, direction_name))

    def run(self, operation_list):
        """
        Runs a list of operations in parallel

        Parameters
        ----------
        operation_list : Operation instance or list of Operation instances
            List of operations to process
        """
        if type(operation_list) != list:
            operation_list = [operation_list]
        self.operation_list = operation_list

        # Run the operation(s)
        while len(self.operation_list) > 0:
            self.allocate_resources()
            with Timer(log, 'operation'):
                pool = multiprocessing.Pool(processes=self.nops_simul)
                self.queued_ops = self.operation_list[self.nops_simul:]
                process_list = []
                for op in self.operation_list:
                    op.setup()
                    process_list.append(
                        pool.apply_async(call_toil, (op.name,
                                         op.direction.name, op.pipeline_parset_file,
                                         op.pipeline_inputs_file, op.pipeline_working_dir,
                                         self.scratch_dir, op.logbasename, self.batch_system),
                                         callback=self.result_callback)
                                        )
                pool.close()
                pool.join()

            # Check for and handle any failed ops
            if not self.success:
                log.error('One or more operations failed due to an error. Exiting...')
                sys.exit(1)
            else:
                self.operation_list = []
def call_toil(op_name, direction_name, parset, inputs, basedir, dir_local, logbasename,
              batch_system):
