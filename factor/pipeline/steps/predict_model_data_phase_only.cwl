cwlVersion: v1.0
class: CommandLineTool
baseCommand: [DPPP]
label: "Calibrates a dataset using DDECal"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - numthreads=0
  - msin.datacolumn=DATA
  - msout.overwrite=True
  - msout.writefullresflag=False
  - steps=[predict]
  - predict.type=h5parmpredict
  - predict.operation=replace
  - predict.applycal.correction=phase000
  - predict.applycal.steps=[fastphase]
  - predict.applycal.fastphase.correction=phase000
  - predict.usebeammodel=True
  - predict.beammode=array_factor
  - predict.onebeamperpatch=True
  - msout.storagemanager=Dysco
  - msout.storagemanager.databitrate=0

inputs:
  - id: msin
    type: string
    inputBinding:
      prefix: msin=
      separate: False
  - id: msout
    type: string
    inputBinding:
      prefix: msout=
      separate: False
  - id: starttime
    type: string
    inputBinding:
      prefix: msin.starttime=
      separate: False
  - id: ntimes
    type: int
    inputBinding:
      prefix: msin.ntimes=
      separate: False
  - id: h5parm
    type: string
    inputBinding:
      prefix: predict.applycal.parmdb=
      separate: False
  - id: sourcedb
    type: string
    inputBinding:
      prefix: predict.sourcedb=
      separate: False
  - id: directions
    type: string[]
    inputBinding:
      valueFrom: $('['+self+']')
      prefix: predict.directions=
      itemSeparator: ','
      separate: False

outputs:
  - id: msmod
    type: string
    outputBinding:
      outputEval: $(inputs.msout)
