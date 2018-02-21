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
        if len(field.sectors) > 1:
            # Use the model-subtracted data for this sector
            obs_filename = get_obs_parameters('ms_subtracted_filename')
        else:
            obs_filename = get_obs_parameters('ms_filename')
        image_freqstep = sector.get_obs_parameters('image_freqstep')
        image_timestep = sector.get_obs_parameters('image_timestep')

        self.parms_dict.update({'obs_filename': obs_filename,
                                'h5parm_mapfile': field.h5parm_mapfile,
                                'image_freqstep': image_freqstep,
                                'image_timestep': image_timestep})


    def finalize(self):
        """
        Finalize this operation
        """
        # Add output datamaps to field object for later use
        if hasattr(self.field, 'image_mapfile'):
            self.field.image_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
                'sector_image.mapfile')
            self.field.premask_mapfile['facetimage'] = os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')
        else:
            self.field.image_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
                'sector_image.mapfile')}
            self.field.premask_mapfile = {'facetimage': os.path.join(self.pipeline_mapfile_dir,
                'premask.mapfile')}
