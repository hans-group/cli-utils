#!/bin/bash
function gotojob {
    if (( $# == 0 )); then
        echo "gotojob v0.1.0 by Minjoon Hong"
        echo "change current directory to working directory of a slurm job"
        echo "usage: $FUNCNAME [SLURM_JOBID]"
        return
    fi
    workdir=$(scontrol show job $1 | grep Dir | awk -F '=' '{print $2}')
    cd $workdir
}
