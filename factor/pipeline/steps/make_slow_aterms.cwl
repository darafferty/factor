cwlVersion: v1.0
class: CommandLineTool
baseCommand: [make_aterm_images.py]
label: "Make FITS images of the aterms"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - '--soltabname=gain000'
  - '--smooth_deg=0.05'
  - '--gsize_deg=0.05'
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

outputs: []
