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
    def __init__(self, field, index):
        super(Image, self).__init__(field, name='image', index=index)

        # For imaging we use a subworkflow, so we set the template filename for that here
        self.subpipeline_parset_template = '{0}_sector_pipeline.cwl'.format(self.rootname)

    def set_parset_parameters(self):
        """
        Define parameters needed for the pipeline parset template
        """
        self.parset_parms = {'factor_pipeline_dir': self.factor_pipeline_dir,
                             'pipeline_working_dir': self.pipeline_working_dir,
                             'do_slowgain_solve': self.field.do_slowgain_solve,
                             'use_beam': self.field.use_beam}

    def set_input_parameters(self):
        """
        Define the pipeline inputs
        """
        nsectors = len(self.field.sectors)
        obs_filename = []
        prepare_filename = []
        mask_filename = []
        aterms_config_file = []
        starttime = []
        ntimes = []
        aterm_image_filenames = []
        image_freqstep = []
        image_timestep = []
        multiscale_scales_pixel = []
        local_dir = []
        phasecenter = []
        image_root = []
        for sector in self.field.sectors:
            if nsectors > 1:
                # Use the model-subtracted data
                obs_filename.append(sector.get_obs_parameters('ms_subtracted_filename'))
            else:
                obs_filename.append(sector.get_obs_parameters('ms_filename'))

            # Each image job must have its own directory, so we create it here
            image_dir = os.path.join(self.pipeline_working_dir, sector.name)
            misc.create_directory(image_dir)

            image_root.append(os.path.join(image_dir, sector.name))
            prepare_filename.append([os.path.join(image_dir, os.path.basename(of)+'.prep')
                                     for of in obs_filename])
            mask_filename.append(image_root[-1] + '_mask.fits')
            aterms_config_file.append(image_root[-1] + '_aterm.cfg')
            image_freqstep.append(sector.get_obs_parameters('image_freqstep'))
            image_timestep.append(sector.get_obs_parameters('image_timestep'))
            sector_starttime = []
            sector_ntimes = []
            for obs in self.field.observations:
                sector_starttime.append(obs.convert_mjd(obs.starttime))
                sector_ntimes.append(obs.numsamples)
            starttime.append(sector_starttime)
            ntimes.append(sector_ntimes)
            phasecenter.append("'[{0}deg, {1}deg]'".format(sector.ra, sector.dec))
            if self.temp_dir is None:
                local_dir.append(image_dir)
            else:
                local_dir.append(self.temp_dir)
            multiscale_scales_pixel.append("'{}'".format(sector.multiscale_scales_pixel))

            # The following attribute was set by the preceding calibrate operation
            aterm_image_filenames.append("'[{}]'".format(','.join(self.field.aterm_image_filenames)))

        self.input_parms = {'obs_filename': obs_filename,
                            'prepare_filename': prepare_filename,
                            'mask_filename': mask_filename,
                            'aterms_config_file': aterms_config_file,
                            'starttime': starttime,
                            'ntimes': ntimes,
                            'aterm_image_filenames': aterm_image_filenames,
                            'image_freqstep': image_freqstep,
                            'image_timestep': image_timestep,
                            'phasecenter': phasecenter,
                            'image_name': image_root,
                            'multiscale_scales_pixel': multiscale_scales_pixel,
                            'local_dir': local_dir,
                            'do_slowgain_solve': [self.field.do_slowgain_solve] * nsectors,
                            'channels_out': [sector.wsclean_nchannels for sector in self.field.sectors],
                            'ra': [sector.ra for sector in self.field.sectors],
                            'dec': [sector.dec for sector in self.field.sectors],
                            'wsclean_imsize': [sector.imsize for sector in self.field.sectors],
                            'vertices_file': [sector.vertices_file for sector in self.field.sectors],
                            'region_file': [sector.region_file for sector in self.field.sectors],
                            'use_beam': [sector.use_beam for sector in self.field.sectors],
                            'wsclean_niter': [sector.wsclean_niter for sector in self.field.sectors],
                            'robust': [sector.robust for sector in self.field.sectors],
                            'wsclean_image_padding': [sector.wsclean_image_padding for sector in self.field.sectors],
                            'cellsize_deg': [sector.cellsize_deg for sector in self.field.sectors],
                            'min_uv_lambda': [sector.min_uv_lambda for sector in self.field.sectors],
                            'max_uv_lambda': [sector.max_uv_lambda for sector in self.field.sectors],
                            'taper_arcsec': [sector.taper_arcsec for sector in self.field.sectors],
                            'auto_mask': [sector.auto_mask for sector in self.field.sectors],
                            'idg_mode': [sector.idg_mode for sector in self.field.sectors],
                            'threshisl': [sector.threshisl for sector in self.field.sectors],
                            'threshpix': [sector.threshpix for sector in self.field.sectors]}

    def finalize(self):
        """
        Finalize this operation
        """
        # Save output FITS image and model for each sector
        # NOTE: currently, -save-source-list only works with pol=I -- when it works with other
        # pols, save them all
        for sector in self.field.sectors:
            image_root = os.path.join(self.pipeline_working_dir, sector.name)
            sector.I_image_file = image_root + '-MFS-image-pb.fits'
            sector.I_model_file = image_root + '-MFS-model-pb.fits'

            # The sky models, both true sky and apparent sky (the filenames are defined
            # in the factor/scripts/filter_skymodel.py file)
            sector.image_skymodel_file_true_sky = image_root + '.true_sky'
            sector.image_skymodel_file_apparent_sky = image_root + '.apparent_sky'

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
