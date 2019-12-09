cwlVersion: v1.0
class: CommandLineTool
baseCommand: [wsclean]
label: "Images a dataset using WSClean+IDG"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - -no-update-model-required
  - -multiscale
  - -fit-beam
  - -reorder
  - -save-source-list
  - -local-rms
  - -join-channels
  - -use-idg
  - -pol I
  - -mgain 0.6
  - -deconvolution-channels 5
  - -fit-spectral-pol 4
  - -multiscale-shape gaussian
  - -weighting-rank-filter 3
  - -auto-threshold 1.0
  - -local-rms-window 50
  - -local-rms-method rms-with-min
  - -aterm-kernel-size 32
  - -nmiter 8

inputs:
  - id: msin
    type: string[]
    inputBinding:
      position: 2
      itemSeparator: " "
  - id: name
    type: string
    inputBinding:
      prefix: -name
      separate: False
  - id: mask
    type: string
    inputBinding:
      prefix: -fits-mask
      separate: False
  - id: config
    type: string
    inputBinding:
      prefix: -aterm-config
      separate: False
  - id: wsclean_imsize
    type: string
    inputBinding:
      prefix: -size
      separate: False
  - id: wsclean_niter
    type: int
    inputBinding:
      prefix: -niter
      separate: False
  - id: robust
    type: float
    inputBinding:
      prefix: -robust
      separate: False
  - id: wsclean_image_padding
    type: float
    inputBinding:
      prefix: -padding
      separate: False
  - id: min_uv_lambda
    type: float
    inputBinding:
      prefix: -minuv-l
      separate: False
  - id: max_uv_lambda
    type: float
    inputBinding:
      prefix: -maxuv-l
      separate: False
  - id: cellsize_deg
    type: float
    inputBinding:
      prefix: -scale
      separate: False
  - id: multiscale_scales_pixel
    type: string
    inputBinding:
      prefix: -multiscale-scales
      separate: False
  - id: local_dir
    type: string
    inputBinding:
      prefix: -temp-dir
      separate: False
  - id: channels_out
    type: int
    inputBinding:
      prefix: -channels-out
      separate: False
  - id: taper_arcsec
    type: float
    inputBinding:
      prefix: -taper-gaussian
      separate: False
  - id: auto_mask
    type: float
    inputBinding:
      prefix: -auto-mask
      separate: False
  - id: idg_mode
    type: string
    inputBinding:
      prefix: -idg-mode
      separate: False

outputs:
  - id: image_nonpb_name
    type: string
    outputBinding:
      outputEval: $(inputs.name)-MFS-image.fits
  - id: skymodel_nonpb
    type: string
    outputBinding:
      outputEval: $(inputs.name)-sources.txt
  - id: skymodel_pb
    type: string
    outputBinding:
      outputEval: $(inputs.name)-sources-pb.txt
