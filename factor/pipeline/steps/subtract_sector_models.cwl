cwlVersion: v1.0
class: CommandLineTool
baseCommand: [subtract_sector_models.py]
label: "Subtracts sector model data"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - '--reweight=True'
  - '--weights_colname=WEIGHT_SPECTRUM'
  - '--phaseonly=True'

inputs:
  - id: msobs
    type: string
    inputBinding:
      position: 0
  - id: msmod
    type: string
    inputBinding:
      position: 1
      itemSeparator: ","
  - id: obs_starttime
    inputBinding:
      prefix: --starttime=
      separate: False
  - id: solint_sec
    inputBinding:
      prefix: --solint_sec=
      separate: False
  - id: solint_hz
    inputBinding:
      prefix: --solint_hz=
      separate: False
  - id: infix
    inputBinding:
      prefix: --infix=
      separate: False
  - id: min_uv_lambda
    inputBinding:
      prefix: --min_uv_lambda=
      separate: False
  - id: max_uv_lambda
    inputBinding:
      prefix: --max_uv_lambda=
      separate: False
  - id: nr_outliers
    inputBinding:
      prefix: --nr_outliers=
      separate: False
  - id: peel_outliers
    inputBinding:
      prefix: --peel_outliers=
      separate: False

outputs: []
