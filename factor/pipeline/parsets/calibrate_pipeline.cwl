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
  - id: sector_bounds_deg
    type: string
  - id: sector_bounds_mid_deg
    type: string
  - id: output_aterms_root
    type: string
  - id: combined_h5parms
    type: string

outputs: []

steps:
  - id: make_sourcedb
    label: make_sourcedb
    run: {{ factor_pipeline_dir }}/steps/make_sourcedb.cwl
    in:
      - id: in
        source: calibration_skymodel_file
      - id: out
        source: calibration_sourcedb
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

{% if do_slowgain_solve %}
# Solve for slow gains

  - id: solve_slow_gains
    label: solve_slow_gains
    run: {{ factor_pipeline_dir }}/steps/ddecal_solve_complexgain.cwl
    in:
      - id: msin
        source: freqchunk_filename
      - id: starttime
        source: slow_starttime
      - id: ntimes
        source: slow_ntimes
      - id: fast_h5parm
        source: combine_fast_phases/outh5parm
      - id: h5parm
        source: output_slow_h5parm
      - id: solint
        source: solint_slow_timestep
      - id: nchan
        source: solint_slow_freqstep
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
      - id: slow_gains_h5parm

  - id: combine_slow_gains
    label: combine_slow_gains
    run: {{ factor_pipeline_dir }}/steps/collect_h5parms.cwl
    in:
      - id: inh5parms
        source: solve_slow_gains/slow_gains_h5parm
      - id: outh5parm
        source: combined_slow_h5parm
    out:
      - id: outh5parm

  - id: make_slow_aterms
    label: make_slow_aterms
    run: {{ factor_pipeline_dir }}/steps/make_slow_aterms.cwl
    in:
      - id: slowh5parm
        source: combine_slow_gains/outh5parm
      - id: fasth5parm
        source: combine_fast_phases/outh5parm
      - id: skymodel
        source: calibration_skymodel_file
      - id: outroot
        source: output_aterms_root
      - id: sector_bounds_deg
        source: sector_bounds_deg
      - id: sector_bounds_mid_deg
        source: sector_bounds_mid_deg
    out: []

  - id: combine_fast_and_slow_h5parms
    label: combine_fast_and_slow_h5parms
    run: {{ factor_pipeline_dir }}/steps/combine_h5parms.cwl
    in:
      - id: inh5parm1
        source: combine_slow_gains/outh5parm
      - id: inh5parm2
        source: combine_fast_phases/outh5parm
      - id: outh5parm
        source: combined_h5parms
    out:
      - id: combinedh5parm

{% else %}
# Don't solve for slow gains

  - id: make_fast_aterms
    label: make_fast_aterms
    run: {{ factor_pipeline_dir }}/steps/make_fast_aterms.cwl
    in:
      - id: fasth5parm
        source: combine_fast_gains/outh5parm
      - id: skymodel
        source: calibration_skymodel_file
      - id: outroot
        source: output_aterms_root
      - id: sector_bounds_deg
        source: sector_bounds_deg
      - id: sector_bounds_mid_deg
        source: sector_bounds_mid_deg
    out: []

{% endif %}