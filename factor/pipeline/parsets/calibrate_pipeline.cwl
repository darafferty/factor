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
  - id: combined_fast_h5parm
    type: string
  - id: output_slow_h5parm
    type: string[]
  - id: combined_slow_h5parm
    type: string
  - id: calibration_skymodel_file
    type: string
  - id: calibration_sourcedb
    type: string
  - id: smoothnessconstraint
    type: float
  - id: maxiter
    type: int
  - id: propagatesolutions
    type: string
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
      - solve_fast_phases/fast_phases_h5parm
    type: string[]

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
      - id: starttime
        source: starttime
      - id: ntimes
        source: ntimes
      - id: h5parm
        source: output_fast_h5parm
      - id: solint
        source: solint_fast_timestep
      - id: nchan
        source: solint_fast_freqstep
      - id: sourcedb
        source: make_sourcedb/sourcedb
      - id: maxiter
        source: maxiter
      - id: propagatesolutions
        source: propagatesolutions
      - id: stepsize
        source: stepsize
      - id: tolerance
        source: tolerance
      - id: uvlambdamin
        source: uvlambdamin
      - id: smoothnessconstraint
        source: smoothnessconstraint
    scatter: [msin, starttime, ntimes, h5parm, solint, nchan]
    scatterMethod: dotproduct
    out:
      - id: fast_phases_h5parm

  - id: combine_fast_phases
    label: combine_fast_phases
    run: {{ factor_pipeline_dir }}/steps/collect_h5parms.cwl
    in:
      - id: inh5parms
        source: solve_fast_phases/fast_phases_h5parm
      - id: outh5parm
        source: combined_fast_h5parm
    out:
      - id: outh5parm
