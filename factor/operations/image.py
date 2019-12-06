"""
Module that holds the Image class
"""
import os
import logging
from factor.lib.operation import Operation
from factor.lib import miscellaneous as misc

log = logging.getLogger('factor:image')


class Image(Operation):
    """
    Operation to image a field sector
    """
    def __init__(self, field, sector, index):
        super(Image, self).__init__(field, direction=sector, name='image', index=index)

    def set_parset_parameters(self):
        """
        Define parameters needed for the pipeline parset template
        """
        self.parset_parms = {'factor_pipeline_dir': self.factor_pipeline_dir,
                             'do_slowgain_solve': self.field.do_slowgain_solve,
                             'use_beam': self.field.use_beam}

    def set_input_parameters(self):
        """
        Define the pipeline inputs
        """
        if len(self.field.sectors) > 1:
            # Use the model-subtracted data for this sector
            obs_filename = self.direction.get_obs_parameters('ms_subtracted_filename')
        else:
            obs_filename = self.direction.get_obs_parameters('ms_filename')
        self.image_root = os.path.join(self.pipeline_working_dir, self.direction.name)
        prepare_filename = [of+'.prep' for of in obs_filename]
        mask_filename = self.image_root + '_mask.fits'
        aterms_config_file = self.image_root + '_aterm.cfg'
        image_freqstep = self.direction.get_obs_parameters('image_freqstep')
        image_timestep = self.direction.get_obs_parameters('image_timestep')
        starttime = []
        ntimes = []
        for obs in self.field.observations:
            starttime.append(obs.convert_mjd(obs.starttime))
            ntimes.append(obs.numsamples)
        phasecenter = '{0} {1}'.format(self.direction.ra, self.direction.dec)
        if self.temp_dir is None:
            local_dir = self.pipeline_working_dir
        else:
            local_dir = self.temp_dir
        aterm_image_filenames = "'[{}]'".format(','.join(self.field.aterm_image_filenames))

        self.input_parms = {'obs_filename': obs_filename,
                            'prepare_filename': prepare_filename,
                            'mask_filename': mask_filename,
                            'aterms_config_file': aterms_config_file,
                            'starttime': starttime,
                            'ntimes': ntimes,
                            'aterm_image_filenames': aterm_image_filenames,
                            'do_slowgain_solve': self.field.do_slowgain_solve,
                            'image_freqstep': image_freqstep,
                            'image_timestep': image_timestep,
                            'channels_out': self.direction.wsclean_nchannels,
                            'phasecenter': phasecenter,
                            'image_name': self.image_root,
                            'ra': self.direction.ra,
                            'dec': self.direction.dec,
                            'wsclean_imsize': self.direction.wsclean_imsize,
                            'vertices_file': self.direction.vertices_file,
                            'region_file': self.direction.region_file,
                            'use_beam': self.direction.use_beam,
                            'wsclean_niter': self.direction.wsclean_niter,
                            'robust': self.direction.robust,
                            'wsclean_image_padding': self.direction.wsclean_image_padding,
                            'cellsize_deg': self.direction.cellsize_deg,
                            'min_uv_lambda': self.direction.min_uv_lambda,
                            'max_uv_lambda': self.direction.max_uv_lambda,
                            'multiscale_scales_pixel': self.direction.multiscale_scales_pixel,
                            'local_dir': local_dir,
                            'taper_arcsec': self.direction.taper_arcsec,
                            'auto_mask': self.direction.auto_mask,
                            'idg_mode': self.direction.idg_mode,
                            'threshisl': self.direction.threshisl,
                            'threshpix': self.direction.threshpix}

    def finalize(self):
        """
        Finalize this operation
        """
        # Save output FITS image and model
        # NOTE: currently, -save-source-list only works with pol=I -- when it works with other
        # pols, save all pols
        self.direction.I_image_file = self.image_root + '-MFS-image-pb.fits'
        self.direction.I_model_file = self.image_root + '-MFS-model-pb.fits'

        # The sky models, both true sky and apparent sky (the filenames are defined
        # in the factor/scripts/filter_skymodel.py file)
        self.direction.image_skymodel_file_true_sky = self.image_root + '.true_sky'
        self.direction.image_skymodel_file_apparent_sky = self.image_root + '.apparent_sky'

        # Symlink to datasets and remove old ones
#         dst_dir = os.path.join(self.parset['dir_working'], 'datasets', self.direction.name)
#         misc.create_directory(dst_dir)
#         ms_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
#                                            'prepare_imaging_data.mapfile'))
#         for ms in ms_map:
#             dst = os.path.join(dst_dir, os.path.basename(ms.file))
#             os.system('ln -fs {0} {1}'.format(ms.file, dst))
#         if self.index > 1:
#             prev_iter_mapfile_dir = self.pipeline_mapfile_dir.replace('image_{}'.format(self.index),
#                                                                       'image_{}'.format(self.index-1))
#             self.cleanup_mapfiles = [os.path.join(prev_iter_mapfile_dir,
#                                      'prepare_imaging_data.mapfile')]
#         self.cleanup()
