cwlVersion: v1.0
class: CommandLineTool
baseCommand: [DPPP, msout=., steps=[solve], solve.type=ddecal]

label: "Calibrates a dataset using DDECal"

inputs:
    msin:
        type: Directory
        inputBinding:
            prefix: "msin="
            separate: False
    msin.datacolumn:
        type: string
        default: "DATA"
        inputBinding:
            prefix: "msin.datacolumn="
            separate: False
    msin.starttime:
        type: string
        inputBinding:
            prefix: "msin.starttime="
            separate: False
    msin.ntimes:
        type: string
        default: "0"
        inputBinding:
            prefix: "msin.ntimes="
            separate: False
    msin.baseline:
        type: string
        default: "*"
        inputBinding:
            prefix: "msin.baseline="
            separate: False
    solve.mode:
        type: string
        default: "tec"
        inputBinding:
            prefix: "solve.mode="
            separate: False
    solve.usebeammodel:
        type: string
        default: "True"
        inputBinding:
            prefix: "solve.usebeammodel="
            separate: False
    solve.beammode:
        type: string
        default: "array_factor"
        inputBinding:
            prefix: "solve.beammode="
            separate: False
    solve.onebeamperpatch:
        type: string
        default: "True"
        inputBinding:
            prefix: "solve.onebeamperpatch="
            separate: False
    solve.h5parm:
        type: string
        inputBinding:
            prefix: "solve.h5parm="
            separate: False
    solve.solint:
        type: string
        default: "1"
        inputBinding:
            prefix: "solve.solint="
            separate: False
    solve.nchan:
        type: string
        default: "1"
        inputBinding:
            prefix: "solve.nchan="
            separate: False
    solve.approximatetec:
        type: string
        default: "True"
        inputBinding:
            prefix: "solve.approximatetec="
            separate: False
    solve.maxapproxiter:
        type: string
        default: "50"
        inputBinding:
            prefix: "solve.maxapproxiter="
            separate: False
    solve.maxiter:
        type: string
        default: "150"
        inputBinding:
            prefix: "solve.maxiter="
            separate: False
    solve.propagatesolutions:
        type: string
        default: "True"
        inputBinding:
            prefix: "solve.propagatesolutions="
            separate: False
    solve.stepsize:
        type: string
        default: "0.2"
        inputBinding:
            prefix: "solve.stepsize="
            separate: False
    solve.tolerance:
        type: string
        default: "0.001"
        inputBinding:
            prefix: "solve.tolerance="
            separate: False
    solve.uvlambdamin:
        type: string
        default: "0"
        inputBinding:
            prefix: "solve.uvlambdamin="
            separate: False
    solve.antennaconstraint:
        type: string
        default: "[]"
        inputBinding:
            prefix: "solve.antennaconstraint="
            separate: False
    solve.sourcedb:
        type: string
        inputBinding:
            prefix: "solve.sourcedb="
            separate: False

outputs:
    msout:
        type: Directory
        outputBinding:
            glob: "$(inputs.msin)"

