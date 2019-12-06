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
        prepare_filename = [of+'.prep' for of in obs_filename]
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

        self.input_parms = {'obs_filename': obs_filename,
                            'prepare_filename': prepare_filename,
                            'starttime': starttime,
                            'ntimes': ntimes,
                            'aterms_mapfile': field.aterms_mapfile,
                            'do_slowgain_solve': field.do_slowgain_solve,
                            'image_freqstep': image_freqstep,
                            'image_timestep': image_timestep,
                            'channels_out': sector.wsclean_nchannels,
                            'phasecenter': phasecenter,
                            'ra': self.direction.ra,
                            'dec': self.direction.dec,
                            'wsclean_imsize': self.direction.wsclean_imsize,
                            'vertices_file': self.direction.vertices_file,
                            'region_file': self.direction.region_file,
                            'use_beam': self.direction.use_beam,
                            'wsclean_niter': self.direction.wsclean_niter,
                            'robust': self.direction.robust,
                            'wsclean_image_padding': self.direction.wsclean_image_padding,
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
        misc.create_directory(dst_dir)
        ms_map = DataMap.load(os.path.join(self.pipeline_mapfile_dir,
                                           'prepare_imaging_data.mapfile'))
        for ms in ms_map:
            dst = os.path.join(dst_dir, os.path.basename(ms.file))
            os.system('ln -fs {0} {1}'.format(ms.file, dst))
        if self.index > 1:
            prev_iter_mapfile_dir = self.pipeline_mapfile_dir.replace('image_{}'.format(self.index),
                                                                      'image_{}'.format(self.index-1))
            self.cleanup_mapfiles = [os.path.join(prev_iter_mapfile_dir,
                                     'prepare_imaging_data.mapfile')]
        self.cleanup()
