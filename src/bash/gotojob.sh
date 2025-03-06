#!/bin/bash
function gotojob {
    if (( $# == 0 )); then
        echo "gotojob v0.1.0 by Minjoon Hong"
        echo "change current directory to working directory of a slurm job"
        echo "usage: $FUNCNAME [SLURM_JOBID]"
        return
    fi
    
    # Store the output of scontrol in a variable
    scontrol_output=$(scontrol show job $1 2>&1)
    
    # Check if there was an error
    if [[ $scontrol_output == *"error"* ]]; then
        echo "Error: Invalid job ID or job not found"
        return 1
    fi
    workdir=$(scontrol show job $1 | grep Dir | awk -F '=' '{print $2}')
    
    # Change directory only if workdir is not empty
    if [[ -n "$workdir" ]]; then
        cd "$workdir"
    else
        echo "Error: Could not determine working directory for job $1"
        return 1
    fi
}
