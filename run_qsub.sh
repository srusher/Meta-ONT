#!/bin/bash

# module load nextflow
# module load Python/3.11.1

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

DATETIME=$(date '+%Y-%m-%d_%H-%M-%S')
EMAIL=$USER"@cdc.gov"

qsub \
    -N meta-ont_$DATETIME \
    -M $EMAIL \
    -m abe \
    -q highmem.q \
    -pe smp 8 \
    -l h_vmem=120G \
    -cwd \
    -o $SCRIPTDIR/qsub_logs \
    -e $SCRIPTDIR/qsub_logs/meta-ont_$DATETIME'_error.log' \
    run_local.sh qsub $SCRIPTDIR $1
