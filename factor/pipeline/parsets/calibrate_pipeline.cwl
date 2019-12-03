cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  - id: timechunk_filename
    items: string[]
  - id: freqchunk_filename
    items: string[]
  - id: starttime
    items: string[]
  - id: ntimes
    items: string[]
  - id: slow_starttime
    items: string[]
  - id: slow_ntimes
    items: string[]
  - id: startchan
    items: string[]
  - id: nchan
    items: string[]
  - id: solint_fast_timestep
    items: string[]
  - id: solint_slow_timestep
    items: string[]
  - id: solint_fast_freqstep
    items: string[]
  - id: solint_slow_freqstep
    items: string[]
  - id: output_fast_h5parm
    items: string[]
  - id: output_slow_h5parm
    items: string[]
  - id: calibration_skymodel_file
    type: string

outputs:
  - id: fast_phases_h5parm
    outputSource:
      - solve_fast_phases/fast_phases_h5parm
    type: string

steps:
  - id: make_sourcedb
    label: make_sourcedb
    run: {{ factor_pipeline_dir }}/steps/make_sourcedb.cwl
    in:
      - id: in
        source:
          - calibration_skymodel_file
    out:
      - id: sourcedb

  - id: solve_fast_phases
    label: solve_fast_phases
    run: {{ factor_pipeline_dir }}/steps/ddecal_fast_phases.cwl
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
      scatter:
        - msin
      out:
        - id: fast_phases_h5parm
