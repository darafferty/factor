cwlVersion: v1.0
class: CommandLineTool
baseCommand: [makesourcedb, format=<, append=False, outtype=blob]

label: "Makes a sourcedb file from a sky model file"

inputs:
    in:
        type: File
        inputBinding:
           prefix: "in="
           separate: False


arguments:
    out:
        type: File
        inputBinding:
           prefix: "out="
           valueFrom: "$(inputs.in).sourcedb"

outputs:
    sourcedb:
        type: File
        outputBinding:
            glob: "$(inputs.out)"
