"""
Module that holds all the Image class
"""
import os
import logging
from factor.lib.operation import Operation

log = logging.getLogger('factor:image')


class Image(Operation):
    """
    Operation to image the field
    """
    def __init__(self, field, sector, index):
        name = 'Image_{0}'.format(index)
        super(Image, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'image_pipeline.parset'

        # Define extra parameters needed for this operation
        sector_filename = []
        sector_skymodel = []
        sector_patches = []
        for sector in field.sectors:
            for i, obs in enumerate(sector.observations):
                sector_filename.append('{0}.sector_{1}'.format(obs.ms_filename, i))
                sector_skymodel.append(sector.skymodel_file)
                sector_patches.append(sector.patches)
        sector_patches = '[{}]'.format(';'.join(sector_patches))
        obs_filename = []
        for obs in sector.observations:
            if len(field.sectors) > 1:
                # Use the model-subtracted data
                obs_filename.append(obs.ms_subtracted_filename)
            else:
                obs_filename.append(obs.ms_filename)

        self.parms_dict.update({'obs_filename': obs_filename,
                                'h5parm_mapfile': field.h5parm_mapfile})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to field object for later use
        self.field.fast_h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'fast_h5parm.mapfile')
        self.field.slow_h5parm_mapfile = os.path.join(self.pipeline_mapfile_dir,
            'slow_h5parm.mapfile')
