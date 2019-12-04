cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  - id: timechunk_filename
    type: string[]
  - id: freqchunk_filename
    type: string[]
  - id: starttime
    type: string[]
  - id: ntimes
    type: int[]
  - id: slow_starttime
    type: string[]
  - id: slow_ntimes
    type: int[]
  - id: startchan
    type: int[]
  - id: nchan
    type: int[]
  - id: solint_fast_timestep
    type: int[]
  - id: solint_slow_timestep
    type: int[]
  - id: solint_fast_freqstep
    type: int[]
  - id: solint_slow_freqstep
    type: int[]
  - id: output_fast_h5parm
    type: string[]
  - id: output_slow_h5parm
    type: string[]
  - id: calibration_skymodel_file
    type: string
  - id: calibration_sourcedb
    type: string
  - id: smoothnessconstraint
    type: float
  - id: maxiter
    type: int
  - id: propagatesolutions
    type: bool
  - id: stepsize
    type: float
  - id: tolerance
    type: float
  - id: uvlambdamin
    type: float

#outputs:
#  - id: fast_phases_h5parm
#    outputSource:
#      - solve_fast_phases/fast_phases_h5parm
#    type: string
outputs:
  - id: sdb
    outputSource:
      - make_sourcedb/sourcedb
    type: string

steps:
  - id: make_sourcedb
    label: make_sourcedb
    run: {{ factor_pipeline_dir }}/steps/make_sourcedb.cwl
    in:
      - id: in
        source:
          - calibration_skymodel_file
      - id: out
        source:
          - calibration_sourcedb
    out:
      - id: sourcedb

  - id: solve_fast_phases
    label: solve_fast_phases
    run: {{ factor_pipeline_dir }}/steps/ddecal_solve_scalarphase.cwl
    in:
      - id: msin
        source: timechunk_filename
      - id: msin.datacolumn
        source: data_colname
      - id: msin.starttime
        source: starttime
      - id: msin.ntimes
        source: ntimes
      - id: solve.h5parm
        source: output_fast_h5parm
      - id: solve.solint
        source: solint_fast_timestep
      - id: solve.nchan
        source: solint_fast_freqstep
      - id: solve.sourcedb
        source: make_sourcedb/sourcedb
      - id: solve.mode
        source: mode
      - id: solve.approximatetec
        source: approximatetec
      - id: solve.maxapproxiter
        source: maxapproxiter
      - id: solve.maxiter
        source: maxiter
      - id: solve.propagatesolutions
        source: propagatesolutions
      - id: solve.stepsize
        source: solint_fast_freqstep
      - id: solve.tolerance
        source: tolerance
      - id: solve.uvlambdamin
        source: uvlambdamin
      - id: solve.smoothnessconstraint
        source: smoothnessconstraint
    scatter:
      - msin
    out:
      - id: fast_phases_h5parm
