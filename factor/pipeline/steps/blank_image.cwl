baseCommand: [blank_image.py]
label: "Make a blank FITS image for use as a clean mask"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - '--make_blank_image=True'

inputs:
  - id: imagefile
    type: string
    inputBinding:
      position: 0
  - id: maskfile
    type: string
    inputBinding:
      position: 1
  - id: wsclean_imsize
    type: string
    inputBinding:
      prefix: --imsize=
      separate: false
  - id: vertices_file
    type: string
    inputBinding:
      prefix: --vertices_file=
      separate: false
  - id: ra
    type: float
    inputBinding:
      prefix: --reference_ra_deg=
      separate: false
  - id: dec
    type: float
    inputBinding:
      prefix: --reference_dec_deg=
      separate: false
  - id: cellsize_deg
    type: float
    inputBinding:
      prefix: --cellsize_deg=
      separate: false
  - id: region_file
    type: string
    inputBinding:
      prefix: --region_file=
      separate: false

outputs:
  - id: maskimg
    type: string
    outputBinding:
      outputEval: $(inputs.maskfile)
