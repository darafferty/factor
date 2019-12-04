cwlVersion: v1.0
class: CommandLineTool
baseCommand: [DPPP]
label: "Calibrates a dataset using DDECal"

requirements:
  InlineJavascriptRequirement: {}

arguments:
  - numthreads=0
  - msin.datacolumn=DATA
  - msout=.
  - steps=[solve]
  - solve.type=ddecal
  - solve.mode=scalarphase
  - solve.usebeammodel=True
  - solve.beammode=array_factor
  - solve.onebeamperpatch=True

inputs:
  - id: msin
    type: string
    inputBinding:
      prefix: msin=
      separate: False
  - id: msin.starttime
    type: string
    inputBinding:
      prefix: msin.starttime=
      separate: False
  - id: msin.ntimes
    type: int
    inputBinding:
      prefix: msin.ntimes=
      separate: False
  - id: solve.h5parm
    type: string
    inputBinding:
      prefix: solve.h5parm=
      separate: False
  - id: solve.solint
    type: int
    inputBinding:
      prefix: solve.solint=
      separate: False
  - id: solve.nchan
    type: int
    inputBinding:
      prefix: solve.nchan=
      separate: False
  - id: solve.sourcedb
    type: string
    inputBinding:
      prefix: solve.sourcedb=
      separate: False
  - id: solve.maxiter
    type: int
    inputBinding:
      prefix: solve.maxiter=
      separate: False
  - id: solve.propagatesolutions
    type: bool
    inputBinding:
      prefix: solve.propagatesolutions=
      separate: False
  - id: solve.stepsize
    type: float
    inputBinding:
      prefix: solve.stepsize=
      separate: False
  - id: solve.tolerance
    type: float
    inputBinding:
      prefix: solve.tolerance=
      separate: False
  - id: solve.uvlambdamin
    type: float
    inputBinding:
      prefix: solve.uvlambdamin=
      separate: False
  - id: solve.smoothnessconstraint
    type: float
    inputBinding:
      prefix: solve.smoothnessconstraint=
      separate: False

outputs:
  - id: fast_phases_h5parm
    type: string
    outputBinding:
      outputEval: $(inputs.solve.h5parm)
