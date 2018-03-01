"""
Module that holds all the Image class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap

log = logging.getLogger('factor:image')


class Image(Operation):
    """
    Operation to image a field sector
    """
    def __init__(self, field, sector, index):
        name = 'Image_{0}'.format(index)
        super(Image, self).__init__(field, direction=sector, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'image_pipeline.parset'

        # Define extra parameters needed for this operation
        if len(field.sectors) > 1:
            # Use the model-subtracted data for this sector
            obs_filename = sector.get_obs_parameters('ms_subtracted_filename')
        else:
            obs_filename = sector.get_obs_parameters('ms_filename')
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
        # Save output mapfiles for later use
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'image-MFS-image.fits.mapfile'))
        self.direction.store_output_image_filename(in_map[0].file)
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'image-sources.txt.mapfile'))
        self.direction.store_output_skymodel_filename(in_map[0].file)

        # Create sym links to image files
        dst = os.path.join(self.factor_working_dir, 'images',
                           'field-MFS-image_{}.fits'.format(self.index))
        if os.path.exists(dst):
            os.unlink(dst)
        os.symlink(self.direction.get_output_image_filename(), dst)
