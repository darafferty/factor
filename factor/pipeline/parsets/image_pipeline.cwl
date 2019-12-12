cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  - id: obs_filename
    type:
      type: array
      items:
        type: array
        items: string
  - id: prepare_filename
    type:
      type: array
      items:
        type: array
        items: string
  - id: starttime
    type:
      type: array
      items:
        type: array
        items: string
  - id: ntimes
    type:
      type: array
      items:
        type: array
        items: int
  - id: image_freqstep
    type:
      type: array
      items:
        type: array
        items: int
  - id: image_timestep
    type:
      type: array
      items:
        type: array
        items: int
  - id: mask_filename
    type: array
    items: string
  - id: phasecenter
    type: array
    items: string
  - id: ra
    type: array
    items: float
  - id: dec
    type: array
    items: float
  - id: image_name
    type: array
    items: string
  - id: cellsize_deg
    type: array
    items: float
  - id: wsclean_imsize
    type:
      type: array
      items:
        type: array
        items: int
  - id: vertices_file
    type: array
    items: string
  - id: region_file
    type: array
    items: string
  - id: aterms_config_file
    type: array
    items: string
  - id: aterm_image_filenames
    type: array
    items: string
  - id: use_beam
    type: array
    items: string
  - id: channels_out
    type: array
    items: int
  - id: wsclean_niter
    type: array
    items: int
  - id: robust
    type: array
    items: float
  - id: wsclean_image_padding
    type: array
    items: float
  - id: min_uv_lambda
    type: array
    items: float
  - id: max_uv_lambda
    type: array
    items: float
  - id: multiscale_scales_pixel
    type: array
    items: string
  - id: local_dir
    type: array
    items: string
  - id: taper_arcsec
    type: array
    items: float
  - id: auto_mask
    type: array
    items: float
  - id: idg_mode
    type: array
    items: string
  - id: threshisl
    type: array
    items: float
  - id: threshpix
    type: array
    items: float

outputs: []

steps:
  - id: image_sector
    label: image_sector
    run: {{ pipeline_working_dir }}/subpipeline_parset.cwl
    in:
      - id: obs_filename
        source: obs_filename
      - id: prepare_filename
        source: prepare_filename
      - id: starttime
        source: starttime
      - id: ntimes
        source: ntimes
      - id: image_freqstep
        source: image_freqstep
      - id: image_timestep
        source: image_timestep
      - id: mask_filename
        source: mask_filename
      - id: phasecenter
        source: phasecenter
      - id: ra
        source: ra
      - id: dec
        source: dec
      - id: image_name
        source: image_name
      - id: cellsize_deg
        source: cellsize_deg
      - id: wsclean_imsize
        source: wsclean_imsize
      - id: vertices_file
        source: vertices_file
      - id: region_file
        source: region_file
      - id: aterms_config_file
        source: aterms_config_file
      - id: aterm_image_filenames
        source: aterm_image_filenames
      - id: use_beam
        source: use_beam
      - id: channels_out
        source: channels_out
      - id: wsclean_niter
        source: wsclean_niter
      - id: robust
        source: robust
      - id: wsclean_image_padding
        source: wsclean_image_padding
      - id: min_uv_lambda
        source: min_uv_lambda
      - id: max_uv_lambda
        source: max_uv_lambda
      - id: multiscale_scales_pixel
        source: multiscale_scales_pixel
      - id: local_dir
        source: local_dir
      - id: taper_arcsec
        source: taper_arcsec
      - id: auto_mask
        source: auto_mask
      - id: idg_mode
        source: idg_mode
      - id: threshisl
        source: threshisl
      - id: threshpix
        source: threshpix
    scatter: [obs_filename, prepare_filename, starttime, ntimes, image_freqstep,
              image_timestep, mask_filename, phasecenter, ra, dec, image_name,
              cellsize_deg, wsclean_imsize, vertices_file, region_file,
              aterms_config_file, aterm_image_filenames, use_beam, channels_out,
              wsclean_niter, robust, wsclean_image_padding, min_uv_lambda, max_uv_lambda,
              multiscale_scales_pixel, local_dir, taper_arcsec, auto_mask, idg_mode,
              threshisl, threshpix]
    scatterMethod: dotproduct
    out: []
