#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment4.py
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 3.0
Usage : python3 Joshua_Duplaa_R11583588_final_project.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import argparse
import sys
import os

####Task #1 – Reading and Sorting of FASTA formatted Sequence Data


##Task 2 - sequence clustering
"""
This task will receive the sorted sequence data from Task #1 as well as the scoring matrix given via the command
line argument ‘-s’ which will define the substitution scores that should be utilized by the NW algorithm.
This task will produce the intermediary file <output_file>.clustered.fna where <output_file> is the path and
filename given by the command line argument ‘-o
"""

##Task 3 - Chimera Detection and Filtering
"""
This task will receive the clustered sequences generated during Task #2 as well as the scoring matrix given via
the command line argument ‘-s’ which will define the substitution scores that should be utilized by the SW
algorithm.This task will produce two intermediary files where <output_file> is the path and filename given by the
command line argument ‘-o’:
<output_file>.nonchimeric.fna: A fasta file containing the representative sequences determined likely
to not be chimeras.
<output_file>.chimeric.fna: A fasta file containing the representative sequences determined likely to be
chimeras.
"""

##Task 4 - Multiple sequence Alignment
"""
This task will receive the non-chimeric representative sequences from Task #3 as well as the scoring matrix given
via the command line argument ‘-s’ which will define the substitution scores that should be utilized by the NW
algorithm.
This task will produce the intermediary file <output_file>.msa.fna where <output_file> is the path and filename
given by the command line argument ‘-o’.
"""



#Execute task_2 function, return's clustered sequences and representative sequence list

#Execute task_3 Filter out chimeric sequeces from representative sequence list

#Execute task_4 Perform MSA of representative sequence list.