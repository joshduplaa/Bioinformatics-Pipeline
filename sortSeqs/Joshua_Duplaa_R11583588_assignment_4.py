#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment4.py
Description : This is an implementation Feng Doolittle for progressive MSA using NW, Fitch Margialosh clustering, and sum of pairs.
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 2.0
Usage : python3 Joshua_Duplaa_R11583588_assignment4.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import argparse
import sys
import os
import math 
#modified global aligmment script from assignment 2, returns aligned sequences and score
import globalAlign


def ReadAndStoreSequences(args):
    #sequenceLabel has the title, sequenceDict has title:sequence sequenceDict[title]=sequence
    #validate the input file path
    if not os.path.exists(args.input):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)

    with open(args.input, "r") as file:
        sequenceFile = [line.strip() for line in file]

    #dictionary for quick sequence lookup
    labelToSeq = {}
    sequenceList = []
    titleList = []
    seqIndex = -1
 
    for line in sequenceFile:
        if line[0] == ">":
            #skip line and set label To Seq
            label = line
            sequence = ""
            sequenceList.append("")
            titleList.append(line)
            seqIndex += 1
        else:
            sequence += line
            labelToSeq[label] = sequence
            sequenceList[seqIndex] = sequence
    
    seqToLabel = {value: key for key, value in labelToSeq.items()}
    sortedSequenceList = sorted(sequenceList, key=len, reverse=True)

    return sortedSequenceList, labelToSeq, seqToLabel, titleList

##############main implementation
def main():
    print("Assignment4 :: R#11583588")

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

    #task 1, read in the sequences
    sequenceList, labelToSeq, seqToLabel, titleList = Task_1(args)

    #task 2, generate distance matrix with all v all global alignemnt
    distanceMatrix = Task_2(sequenceList, gapPenalty, score_matrix, scoringType)

    #task 3, using Fitch Margoliash construct guide tree using FM clustering, return the collapse order of the distance matrix
    collapseOrder = Task_3(distanceMatrix)

    #task 4, perform MSA given collapse Order Using feng doolittle
    alignedSequences = Task_4(collapseOrder, sequenceList, gapPenalty, score_matrix)

    #task 5, calculate score of alignedSequences (Sum of Pairs!)
    msaScore = Task_5(alignedSequences,gapPenalty, score_matrix)

    #task 6, output to fasta file
    Task_6(alignedSequences, msaScore, args, titleList)
    return 0


if __name__ == "__main__":
    main()