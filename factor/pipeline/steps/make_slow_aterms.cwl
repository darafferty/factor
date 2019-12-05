cwlVersion: v1.0
class: CommandLineTool
baseCommand: [make_aterm_images.py]
label: "Make FITS images of the aterms"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - '--soltabname=gain000'
  - '--smooth_deg=0.1'
  - '--gsize_deg=0.1'
  - '--time_avg_factor=1'

inputs:
  - id: slowh5parm
    type: string
    inputBinding:
      position: 0
  - id: outroot
    type: string
    inputBinding:
      prefix: --outroot=
      separate: false
  - id: fasth5parm
    type: string
    inputBinding:
      prefix: --fasth5parm=
      separate: false
  - id: skymodel
    type: string
    inputBinding:
      prefix: --skymodel=
      separate: false
  - id: sector_bounds_deg
    type: string
    inputBinding:
      prefix: --bounds_deg=
      separate: false
  - id: sector_bounds_mid_deg
    type: string
    inputBinding:
      prefix: --bounds_mid_deg=
      separate: false

outputs:
  - id: aterm_images
    doc: Output aterm images
    type: string
    outputBinding:
      glob: $(inputs.outroot)*fits
  - id: aterm_paths_file
    doc: Output Measurement Set
    type: string
    outputBinding:
      outputEval: $(inputs.outroot).txt
