cwlVersion: v1.0
class: CommandLineTool
baseCommand: [make_aterm_config.py]
label: "Make an aterms configuration file for use in imaging"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  - id: outfile
    type: string
    inputBinding:
      position: 1
  - id: gain_filenames
    type: string
    inputBinding:
      prefix: --gain_filenames=
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
