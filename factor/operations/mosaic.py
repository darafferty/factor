"""
Module that holds the Mosaic class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.utilities import create_directory

log = logging.getLogger('factor:mosaic')


class Mosaic(Operation):
    """
    Operation to mosaic sector images
    """
    def __init__(self, field, index):
        name = 'Mosaic_{0}'.format(index)
        super(Mosaic, self).__init__(field, name=name)
        self.index = index

        # Set the pipeline parset to use
        self.pipeline_parset_template = 'mosaic_pipeline.parset'

        # Define extra parameters needed for this operation
        if len(field.sectors) > 1:
            skip_processing = False
        else:
            skip_processing = True
        sector_image_filename = []
        sector_vertices_filename = []
        for sector in field.sectors:
            sector_image_filename.append(sector.I_image_file)
            sector_vertices_filename.append(sector.vertices_file)
        # TODO: make mosaic of resid, model, and QUV?

        self.parms_dict.update({'skip_processing': skip_processing,
                                'sector_image_filename': sector_image_filename,
                                'sector_vertices_filename': sector_vertices_filename})

    def finalize(self):
        """
        Finalize this operation
        """
        # Save the FITS image and model
        dst_dir = os.path.join(self.field.parset['dir_working'], 'images',
                               'image_{}'.format(self.index))
        create_directory(dst_dir)
        self.field.field_image_filename = os.path.join(dst_dir, 'field-MFS-I-image.fits')
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'make_mosaic.mapfile'))
        os.system('cp {0} {1}'.format(in_map[0].file, self.field.field_image_filename))

        # TODO: make mosaic of model + QUV?
#         self.field_model_filename = os.path.join(dst_dir, 'field-MFS-I-model.fits')

        # Delete temp data
        self.cleanup_mapfiles = [os.path.join(self.pipeline_mapfile_dir,
                                 'regrid_image.mapfile'),
                                 os.path.join(self.pipeline_mapfile_dir,
                                 'make_mosaic_template.mapfile')]
        self.cleanup()
