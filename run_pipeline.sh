#!/bin/bash

# the nfcore template requires a newer version of python to run its input check scripts 
module load Python/3.11.1

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# primary nextflow execution command: alter as needed
nextflow run main.nf -profile singularity,local

echo -ne "\nWould you like to clear your nextflow \"work\" and \".nextflow\" directories?\nNOTE: You only need to keep these directories around if you intend to troubleshoot the previous nextflow run!\nPlease enter Yes or No: "; read user_input; echo


user_input=$(echo "$user_input" | tr '[:upper:]' '[:lower:]')

if [[ $user_input == "yes" || $user_input == "y" ]]; then

    rm -rf $SCRIPTDIR/work/*
    rm -rf $SCRIPTDIR/.nextflow/*
    
fi
