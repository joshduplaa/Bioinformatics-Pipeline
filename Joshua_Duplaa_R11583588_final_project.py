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
import util.globalAlign as globalAlign
import util.localAlign as localAlign

####Task #1 – Sort sequence data from longest to shortest
def Task_1(args):
    sortedSequenceList, labelToSeq, seqToLabel, titleList = readsortseqs.main(args)
    return sortedSequenceList, labelToSeq, seqToLabel, titleList


##Task 2 - sequence clustering
"""
Sequence clustering is grouping high similar sequences together and then use single sequence to represent them all.
Task 2 function "clusters" Sequences by to do the folowing:
    -Take the sorted list of sequences and mark the first(longest) sequence as the first representative sequence.
    -Perform a global aligmnent for the representative sequence to every other sequence.
    -During alignment, evaluate percent identity, if percentidentity>=97%, add it to the the current repsequence cluster.
    -Create new clusters for sequences that did not get added to a representative sequence list.
    -iterate for every cluster until no new clusters are added.
    -evaluate representative sequences for each new cluster created with 
"""
def Task_2(sortedSequenceList, gapPenalty, score_matrix):
    #clusters are a list of lists, list of individual clusters. Technically a jagged 2d array.
    #first seq is the first representative sequence
    clusters = []
    repSeqs = [sortedSequenceList[0]]
    #the first element of each cluster should be the repSeq, with thr repSeqs hashset index = cluster index 
    #One to all alignment
    sequenceList = [x for x in sortedSequenceList]
    sequenceList.remove(repSeqs[0])

    for i in range(len(sortedSequenceList)):
        repSeq = repSeqs[i]
        repSeqCluster, newCluster = buildCluster(repSeq, sequenceList, gapPenalty, score_matrix)
        repSeqCluster = [repSeq] + repSeqCluster
        clusters.append(repSeqCluster)
        if len(newCluster) == 0:
            break
        repSeqs.append(newCluster[0]) #adds new representative sequence
        sequenceList = newCluster[1:]
        
    return repSeqs, clusters

def buildCluster(repSeq, sequenceList, gapPenalty, score_matrix):
    repSeqCluster = []
    newCluster = []
    for seq in sequenceList:
        sequencePair = [repSeq, seq]
        alignedSequences, alignmentScore = globalAlign.main(sequencePair, gapPenalty, score_matrix)
        percentIdentity = CalculatePercentIdentity(alignedSequences)
        if percentIdentity>=97:
            repSeqCluster.append(seq)
        else:
            newCluster.append(seq)

    return repSeqCluster, newCluster

def CalculatePercentIdentity(alignedSequences):
        identicalAlignments = 0
        for index in range(len(alignedSequences[0])):
            nucleotide1 = alignedSequences[0][index]
            nucleotide2 = alignedSequences[1][index]
            if nucleotide1 == nucleotide2 and nucleotide1 != "-":
                identicalAlignments += 1
        
        percentIdentity = (identicalAlignments/len(alignedSequences[0]))*100
        return percentIdentity
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
    2. For each representative sequence (A), determine: if there exist good alignments (defined above)
        between it and two other representative sequences (B and C) where the cluster size of B > A and the cluster size of
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

def  Task_3(repSeqs, clusters, gapPenalty, score_matrix):
    #goal of the task, classify clusters as nonchimeric or chimeric

    goodAlignments = []
    for i in range(len(repSeqs)-1):
        seq1 = repSeqs[i]
        for j in range(i+1, len(repSeqs)):
            seq2 = repSeqs[j]
            alignedSequences, alignmentScore = localAlign.main([seq1, seq2], gapPenalty, score_matrix)
            alignmentLen = len(alignedSequences[0])
            print(alignmentLen), print(alignmentScore)
            if(alignmentLen<60 or alignmentScore<50):
                continue
            else:
                goodAlignments.append([alignedSequences])
    
                
    
    print(repSeqs)
    print(clusters)


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
    repSeqs, clusters = Task_2(sortedSequenceList, gapPenalty, score_matrix)

    print("pause for debugging")
    #task 3
    nonChimericSequences = Task_3(repSeqs, clusters, gapPenalty, score_matrix)

    #task 4
    #alignedNonChiSeqs = Task_4()



if __name__ == "__main__":
    main()