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
import util.readsortseqs as readsortseqs 

####Task #1 – Sort sequence data from longest to shortest
def Task_1(args):
    sortedSequenceList, labelToSeq, seqToLabel, titleList = readsortseqs.main(args)
    return sortedSequenceList


##Task 2 - sequence clustering
"""
Sequence clustering is grouping high similar sequences together and then use single sequence to represent them all.
Task 2 function "clusters" Sequences by to do the folowing:
    -Take the sorted list of sequences and mark the first(longest) sequence as the first representative sequence.
    -Perform a one to many alignment for the representative sequence to every other sequence.
    -We are adding sequences to a list of representative sequences that have a percent identity >= 97% of the representative sequence
    -percent identity = count_of_aligneNucleotides/len(alignedSequence)
    -return list of representative sequences.
"""
def Task_2():


    return repSequences

##Task 3 - Chimera Detection and Filtering
"""
Conditions for a chimeric sequence::
If a sequence (from a cluster with fewer members) aligns to two or more different sequences (from clusters with a higher member count), it
is likely a chimera. 
A high-level algorithm for this process is as follows:
    1. For each representative sequence, calculate pairwise alignments using the Smith-Waterman algorithm
        with all other representative sequences.
            a. Reject and remove any alignment that fails to meet the following criteria:
                i. Good alignments must reach an alignment score ≥50.
                ii. Good alignments must be at least 60 bases in length.
    2. For each representative sequence (A), determine if there exist good alignments (defined above)
        between it and two other representative sequences (B and C) such that the size of B > A and the size of
        C is >A.
        a. In other words, the number of sequences in clusters B and C are greater than the number in
        cluster A.

    3. If such a match can be found, mark the cluster as possibly chimeric (is a chimera).
        a. Write all chimeric sequences to the output file <output_file>.chimeric.fna defined in the Output
        section below.

    4. If no such match can be found, mark the cluster as likely nonchimeric (not a chimera).
        a. Write all non-chimeric sequences to the output file <output_file>.nonchimeric.fna defined in the
        Output section below.
"""

def Task_3():

    return nonChimericSequences


##Task 4 - Multiple sequence Alignment
"""
Perform an MSA of non chimeric sequences
"""
def Task_4(nonChimericSequences):

    return alignedNonChiSeqs


##############main implementation
def main():
    print("Final_Project :: R#11583588")

    parser = argparse.ArgumentParser(description="sorting Genome Sequences in descending order")
    parser.add_argument("-i", "--input", required=True, type=str, help="File path to input.fna")
    parser.add_argument("-o", "--output", required=True, type=str, help="File path for output.fna")
    parser.add_argument("-s", "--score_matrix", required=True, type=str, help="File path for scoring matrix")
    args = parser.parse_args()

    if args.score_matrix == "input/nucleotide.mtx":
        scoringType = "nucleotide"
    elif args.score_matrix == "input/BLOSUM50.mtx":
        scoringType = "BLOSUM50"
    elif args.score_matrix == "input/BLOSUM62.mtx":
        scoringType = "BLOSUM62"
    
    #read in score matrix
    with open(args.score_matrix, "r") as file:
        score_matrix = file.readlines()
    scoreIndex = 0
    for row in score_matrix:
        row = row.strip("\n")
        score_matrix[scoreIndex] = row
        scoreIndex += 1
    score_matrix[0] = "0" + score_matrix[0]
    score_matrix = [row.split() for row in score_matrix]
    gapPenalty = int(score_matrix[3][len(score_matrix[3])-1])

    #task 1
    sortedSequenceList, labelToSeq, seqToLabel, titleList = Task_1(args)

    #task 2
    #repSequences = Task_2()

    #task 3
    #nonChimericSequences = Task_3()

    #task 4
    #alignedNonChiSeqs = Task_4()



if __name__ == "__main__":
    main()