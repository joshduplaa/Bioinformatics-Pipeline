#!/bin/bash
# filename          : tester.sh
# description       : Executes unit tests for assignment 6 for CS5393
# author            : Eric Rees
# email             : Eric.Rees@ttu.edu
# date              : November 16, 2024
# version           : 2.0
# usage             : ./testerr.sh
# notes             : This script is built to run on login.hpcc.ttu.edu.
# license           : MIT License
#==============================================================================

echo -e "\n\nCS5393 Assignment 6 Testing Script."

filename="viterbi_A6.py"

if [ -e "$filename" ]; then

    #Load the Bioinformatics conda environment.
    source /lustre/work/errees/conda/etc/profile.d/conda.sh
    conda activate cs5393

    #Append the current working directory to the PYTHONPATH.
    export PYTHONPATH="$PYTHONPATH:$(pwd)"

    echo -e "Executing viterbi_main.py\n"
    
    #Execute the program.
    timeout 60 python3 /lustre/work/errees/courses/cs5393/assignment_6/code/viterbi_main.py

    echo -e "\n"

else
    echo "Error: File '$filename' does not exist."
fi

conda deactivate
