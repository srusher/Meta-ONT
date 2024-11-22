#!/bin/bash
source /etc/profile

module load nextflow
module load Python/3.11.1

if [[ -z $2 ]]; then

    SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

else

    SCRIPTDIR=$2

fi


if [[ -z $3 ]]; then

    opts="$1"

else

    opts="$3"

fi

NXF_WORK=/scicomp/scratch/$USER/nextflow


# primary nextflow execution command: alter as needed
nextflow run main.nf -profile singularity,local $opts
