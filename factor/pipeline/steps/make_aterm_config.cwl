baseCommand: [make_aterm_config.py]
label: "Make an aterms configuration file for use in imaging"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  - id: msin
    type: string
    inputBinding:
      position: 0
  - id: outfile
    type: string
    inputBinding:
      position: 1
  - id: aterms_mapfile
    type: string
    inputBinding:
      prefix: --gain_mapfile=
      separate: false
  - id: use_beam
    type: string
    inputBinding:
      prefix: --use_beam=
      separate: false

outputs:
  - id: aterms_config
    type: string
    outputBinding:
      outputEval: $(inputs.outfile)
