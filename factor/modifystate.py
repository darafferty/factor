"""
Module that modifies the state of the pipelines
"""
from factor.parset import parset_read
from factor.strategy import set_strategy
from factor.lib.field import Field
import logging
import glob
import os
import sys

log = logging.getLogger('factor:state')
logging.getLogger('factor:parset').setLevel(logging.CRITICAL)
logging.getLogger('factor:strategy').setLevel(logging.CRITICAL)


def check_operation(operation):
    """
    Returns list of started/completed pipeline names for given operation path

    Parameters
    ----------
    operation : str
        Path of operation output
    """
    pipelines = []
    directions = glob.glob(os.path.join(operation, '*'))
    for d in directions:
        statefile = os.path.join(d, 'statefile')
        if os.path.exists(statefile):
            pipelines.append('{0}/{1}'.format(os.path.basename(operation), os.path.basename(d)))

    return pipelines


def run(parset_file):
    """
    Modifies the state of one or more pipelines

    Parameters
    ----------
    parset_file : str
        Filename of parset containing processing parameters
    """
    # Read parset
    log.info('Reading parset and checking state...')
    parset = parset_read(parset_file, use_log_file=False, skip_cluster=True)

    # Initialize field object
    field = Field(parset, mininmal=True)
    field.outlier_sectors = [None]
    field.imaging_sectors = [None]

    # Get the processing strategy
    strategy_steps = set_strategy(field)

    # Check each operation for started pipelines
    while True:
        pipelines = []
        for iter, step in enumerate(strategy_steps):
            if step['do_calibrate']:
                operation = os.path.join(parset['dir_working'], 'pipelines', 'calibrate_{}'.format(iter+1))
                pipelines.extend(check_operation(operation))
            if step['do_predict']:
                operation = os.path.join(parset['dir_working'], 'pipelines', 'predict_{}'.format(iter+1))
                pipelines.extend(check_operation(operation))
            if step['do_image']:
                operation = os.path.join(parset['dir_working'], 'pipelines', 'image_{}'.format(iter+1))
                pipelines.extend(check_operation(operation))

        # List pipelines and query user
        print('\nCurrent stratgy: {}'.format(field.parset['strategy']))
        print('\nPipelines:')
        i = 0
        if len(pipelines) == 0:
            print('    None')
        else:
            for p in pipelines:
                i += 1
                print('    {0}) {1}'.format(i, p))
        try:
            while(True):
                p_number_raw = raw_input('Enter number of pipeline to reset or press "q" to quit: ')
                if p_number_raw.lower() == "q":
                    sys.exit(0)
                elif int(p_number_raw) > 0 and int(p_number_raw) <= i:
                    break
                else:
                    print("Please enter a number between 1 and {}".format(i))
        except KeyboardInterrupt:
            sys.exit(0)
        pipeline = pipelines[int(p_number_raw)-1]

        # Ask for confirmation
        try:
            while(True):
                answer = raw_input('Reset all pipelines from {} onwards? '.format(pipeline))
                if (answer.lower() == "n" or answer.lower() == "no" or
                    answer.lower() == "y" or answer.lower() == "yes"):
                    break
                else:
                    print('Please enter "y" or "n"')
        except KeyboardInterrupt:
            sys.exit(0)

        # Reset pipeline states as requested
        if answer.lower() == "y" or answer.lower() == "yes":
            print('Reseting state...')
            for pipeline in pipelines[int(p_number_raw)-1:]:
                statefile = os.path.join(parset['dir_working'], 'pipelines', pipeline, 'statefile')
                os.system('rm -f {}'.format(statefile))
            print('Reset complete.')
