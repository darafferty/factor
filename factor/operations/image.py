"""
Module that holds the Image class
"""
import os
import logging
from factor.lib.operation import Operation
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.utilities import create_directory

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
        starttime = []
        ntimes = []
        for obs in field.observations:
            starttime.append(obs.convert_mjd(obs.starttime))
            ntimes.append(obs.numsamples)

        self.parms_dict.update({'obs_filename': obs_filename,
                                'starttime': starttime,
                                'ntimes': ntimes,
                                'h5parm_mapfile': field.h5parm_mapfile,
                                'fast_aterms_mapfile': field.fast_aterms_mapfile,
                                'slow_aterms_mapfile': field.slow_aterms_mapfile,
                                'image_freqstep': image_freqstep,
                                'image_timestep': image_timestep})

    def finalize(self):
        """
        Finalize this operation
        """
        # Save output mapfiles for later use:
        # The FITS image and model
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'image-MFS-image-pb.fits.mapfile'))
        self.direction.I_image_file = in_map[0].file
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'image-MFS-model-pb.fits.mapfile'))
        self.direction.I_model_file = in_map[0].file
        # NOTE: currently, -save-source-list only works with pol=I -- when it works with other
        # pols, enable IQUV imaging with the following lines
#         in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
#                                            'image-MFS-I-image-pb.fits.mapfile'))
#         self.direction.I_image_file = in_map[0].file
#         in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
#                                            'image-MFS-Q-image-pb.fits.mapfile'))
#         self.direction.Q_image_file = in_map[0].file
#         in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
#                                            'image-MFS-U-image-pb.fits.mapfile'))
#         self.direction.U_image_file = in_map[0].file
#         in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
#                                            'image-MFS-V-image-pb.fits.mapfile'))
#         self.direction.V_image_file = in_map[0].file

        # The sky models, both true sky and apparent sky (the filenames are defined
        # in the factor/scripts/filter_skymodel.py file)
        in_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'filter.mapfile'))
        self.direction.image_skymodel_file_true_sky = in_map[0].file + '.true_sky'
        self.direction.image_skymodel_file_apparent_sky = in_map[0].file + '.apparent_sky'

        # Symlink to datasets and remove old ones
        dst_dir = os.path.join(self.parset['dir_working'], 'datasets', self.direction.name)
        create_directory(dst_dir)
        ms_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'prepare_imaging_data.mapfile'))
        for ms in ms_map:
            dst = os.path.join(dst_dir, os.path.basename(ms.file))
            if os.path.exists(dst):
                os.unlink(dst)
            os.system('ln -s {0} {1}'.format(ms.file, dst))
        if self.index > 1:
            prev_iter_mapfile_dir = self.pipeline_mapfile_dir.replace('image_{}'.format(self.index),
                                                                      'image_{}'.format(self.index-1))
            self.cleanup_mapfiles = [os.path.join(prev_iter_mapfile_dir,
                                     'prepare_imaging_data.mapfile')]
        self.cleanup()
