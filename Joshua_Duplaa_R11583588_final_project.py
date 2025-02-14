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

####Task 1###
def ReadAndStoreSequences(args):
    #sequenceLabel has the title, sequenceDict has title:sequence sequenceDict[title]=sequence
    # Validate the input file path
    if not os.path.exists(args.input):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)

    with open(args.input, "r") as file:
        genome = [line.strip() for line in file]


    #original order sequenceList
    sequenceList = []
    lineToAdd = []

    for line in genome:
        if line[0] == ">":
            if lineToAdd:  # Ensure an empty list isn't appended at the start
                sequenceList.append(''.join(lineToAdd))
            lineToAdd = []
            print("skip")
        else:
            lineToAdd.extend([line])
    
    if lineToAdd:  # Ensure an empty list isn't appended at the start
            sequenceList.append(''.join(lineToAdd))
            


    #sequence Dictionary to store sequence title:index where the same index is the index of the sequence sequenceList 
    sequenceDict = {}
    for sequenceLabel in genome:
        if(">" in sequenceLabel):
            sequenceLength = CountSequence(sequenceLabel, 0, genome)
            sequenceDict[sequenceLabel] = sequenceLength
            print(f"confirmation of sequence length {sequenceDict[sequenceLabel]}")


    labelToSeq = {}
    seqIndex = 0
    for line in genome:
        if line[0] == ">":
            labelToSeq[line] = sequenceList[seqIndex]
            seqIndex += 1

    #Now I need to make a sorted list 
    #found cool python function for sorting the dictionary
    #https://docs.python.org/3/library/functions.html#sorted
    sorted_sequence_data = sorted(sequenceDict.items(), key=lambda x: x[1], reverse=True)
    sortedSequenceList = [""]*len(sorted_sequence_data)
    print(sorted_sequence_data)
    for j in range(len(sorted_sequence_data)):
        sequence = labelToSeq[sorted_sequence_data[j][0]]
        sortedSequenceList[j] += sequence

    print(sortedSequenceList)

    seqToLabel = {value: key for key, value in labelToSeq.items()}

    return sortedSequenceList, labelToSeq, seqToLabel, sorted_sequence_data, genome


def CountSequence(sequenceLabel, sequenceLength, genome):
    #looping through the relavent genome
    #print(genome.index(sequenceLabel))
    for i in range(genome.index(sequenceLabel)+1, len(genome)):
        #check if the next sequence has not been reached. If reached return back to the main function.
        if ">" in genome[i]:
            return sequenceLength
        # Add the last sequence length if we've reached the end of the file
        elif i == len(genome) - 1:
            sequenceLength += len(genome[i])
            return sequenceLength
        else:
            sequenceLength = sequenceLength + len(genome[i])

def outSortedSeqToFile(sorted_sequence_data, genome, labelToSeq):
    fileName = f"{args.output}.sorted.fna"

    #write sorted_sequence_data to file, but only the sequenceLabel in each line. omit the length.

    with open(fileName, "w") as output_file:
        for seq in sorted_sequence_data:
            #sequence header seq[0]
            output_file.write(f"{seq[0]}\n")
            genomeIndex = genome.index(seq[0])
            
            #write sequence data for each header
            for i in range(genomeIndex+1, len(genome)):
                if ">" in genome[i]:
                    break
                # Write each line of the sequence
                output_file.write(f"{genome[i]}\n")

            #special handling for the very last sequence in the file
            if i == len(genome):
                #writing the last line of the last sequence
                output_file.write(genome[i])
            
    return 1


def Task_1(args):
    sortedSequenceList, labelToSeq, seqToLabel, sorted_sequence_data, genome = ReadAndStoreSequences(args)
    outSortedSeqToFile(sorted_sequence_data, genome, labelToSeq)
    return sortedSequenceList, labelToSeq, seqToLabel, sorted_sequence_data, genome

##Task 2 - sequence clustering
"""
I am using the function clusterSequences() to do the folowing:
    -Take the sorted list of sequences and mark the first(longest) sequence as the first representative sequence.
    -Perform a one to many alignment for the representative sequence to every other sequence.
    -We are adding sequences to the representative sequence list that have a percent identity >= 97% of the representative sequence
        -percent identity = count_of_aligneNucleotides/len(alignedSequence)
"""

#function below returns dictionary of sequence clusters, rep_sequence:[sequence1, sequence3, ...]
#also returns representativeSequencesList
def outSequenceCluster(sequenceCluster, seqToLabel):
    #write sorted_sequence_data to file, but only the sequenceLabel in each line. omit the length.
    fileName = f"{args.output}.clustered.fna"

    clustersToRemove = []
        
    with open(fileName, "w") as output_file:
        for cluster in sequenceCluster: 
            title = seqToLabel[cluster]
            size = len(sequenceCluster[cluster])
            if size == 0:
                clustersToRemove.append(cluster)
            else:
                output_file.write(f"{title}; size={size}\n{cluster}\n")

    for item in clustersToRemove:
        del sequenceCluster[item]
            
    return sequenceCluster

def ClusterSequences(repSequence, sortedSequenceList, sequenceClusters):
    #representative sequences are the key in in the sequenceClusters 
    #function builds cluster for every representative sequence, then adds reprentative sequences as new keys in sequence
    queue = [repSequence]
    processed = set()

    while queue:
        currentRep = queue.pop(0)  #take the next sequence to process
        if currentRep in processed:
            continue

        buildCluster(currentRep, sequenceClusters, sortedSequenceList, queue, processed)
        processed.add(currentRep)  #mark as processed


    
    return sequenceClusters

def buildCluster(repSequence, sequenceClusters, sortedSequenceList, queue, processed):
    for sequenceTwo in sortedSequenceList:
        if sequenceTwo == repSequence:
            #print("skip don't align to itself")
            continue
        elif sequenceTwo[0] == ">":
            #print("skip title")
            continue
        elif sequenceTwo in processed: 
            continue  # Skip already classified sequences

        else:
            sequencesToAlign = [repSequence, sequenceTwo]
            #align representative sequence to sequenceTwo
            alignedSequences, alignmentScore = NeedlemanWunsch(sequencesToAlign)
            #calculate percent identity
            percentIdentity = CalculatePercentIdentity(alignedSequences)
            
            if(percentIdentity>97):
                sequenceClusters[repSequence].append(sequenceTwo)
            else:
                sequenceClusters[sequenceTwo] = []
                queue.append(sequenceTwo)
    
        
            #How do I add to the cluster for repSequece Sequence 
    

def CalculatePercentIdentity(alignedSequences):
        identicalAlignments = 0
        for index in range(len(alignedSequences[0])):
            nucleotide1 = alignedSequences[0][index]
            nucleotide2 = alignedSequences[1][index]
            if nucleotide1 == nucleotide2 and nucleotide1 != "-":
                identicalAlignments += 1
        
        percentIdentity = (identicalAlignments/len(alignedSequences[0]))*100
        return percentIdentity

def NeedlemanWunsch(sequencesToAlign):
    #build alignMatrix
    #print(sequencesToAlign)
    alignMatrix = BuildAlignmentMatrix(sequencesToAlign, gapPenalty)
    traceback = BuildTracebackMatrix(sequencesToAlign)
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = CalculateMaxScore(alignMatrix, i, j, gapPenalty, score_matrix)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction
    
    alignmentScore = alignMatrix[len(alignMatrix)-1][len(alignMatrix[0])-1]
    alignedSequences = AlignSequence(traceback)

            
    return alignedSequences, alignmentScore

def BuildAlignmentMatrix(sequencesToAlign, gapPenalty):
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    alignMatrix = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the alignMatrix with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        alignMatrix.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(alignMatrix[0])):
        alignMatrix[0][i] = topSequenceRow[i-2]

    #FWrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in alignMatrix[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal += 1 

    setRowVal = gapPenalty
    #filling in the necessary values -2,-4,-6 etc. for row 2 and beyond in alignMatrix
    for j in range(2,len(alignMatrix[0])):
        alignMatrix[1][j] = setRowVal
        setRowVal += gapPenalty    #decrement -2 everytime

    #filling in the necessary values -2,-4,-6 etc. for col 2 abd beyond in alignMatrix
    setColVal = 0
    for row in alignMatrix[1:]:
        row[1] = setColVal
        setColVal += gapPenalty

    return alignMatrix

def BuildTracebackMatrix(sequencesToAlign):
    #Now I need to build a scoring matrix
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    traceback = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the traceback with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        traceback.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(traceback[0])):
        traceback[0][i] = topSequenceRow[i-2]

    #FWrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in traceback[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal +=1
    
    #filling in the necessary values L (left) for row 2 in traceback
    for j in range(2,len(traceback[0])):
        traceback[1][j] = "L"

    #filling in the necessary values U(up) for col 2 in traceback
    for row in traceback[2:]:
        row[1] = "U"

    return traceback

#Helper function for needleman alg to find alignment score.           
def CalculateMaxScore(alignMatrix, rowIndex, colIndex, gapPenalty, score_matrix):
    #find topLetter and leftLetter so we can find match/mismatch value in score_matrix
    topLetter = alignMatrix[0][colIndex]
    leftLetter = alignMatrix[rowIndex][0]

    scoreIndexTop = score_matrix[0].index(topLetter)
    scoreIndexLeft = score_matrix[0].index(leftLetter)

    #find score value at the score indexes found
    matchValue = int(score_matrix[scoreIndexTop][scoreIndexLeft])

    #Calculate score at the cell in alignMatrix
    diagVal = matchValue + alignMatrix[rowIndex-1][colIndex-1]
    upVal = alignMatrix[rowIndex-1][colIndex] + gapPenalty
    leftVal = alignMatrix[rowIndex][colIndex-1] + gapPenalty

    cellScore = max(diagVal, upVal, leftVal)
    #gap penalty found in scoring file chosen. For now just focus on the nucleatide scoring file (score_matrix)
    direction = ""
    if(cellScore == diagVal):
        direction = "D"
    elif(cellScore == upVal):
        direction = "U"
    elif(cellScore == leftVal):
        direction = "L"

    return cellScore, direction

#Helper function for needleman alg to find alignment alignment between sequences
def AlignSequence(traceback):
    #Finished stepping through matrix when arrived at traceback[1][1], starting point is the bottom left traceback[len(traceback)-1][len(traceback[0])-1]
    startRow = len(traceback)-1
    startCol = len(traceback[0])-1
    i = startRow
    j = startCol

    leftSequence = ""
    topSequence = ""

    while i != 1 and j != 1:
        if(traceback[i][j] == "D"):
            #Write both top and left into sequence listtopSequece
            leftSequence = traceback[i][0]+leftSequence
            topSequence = traceback[0][j]+topSequence 
            i -= 1
            j -= 1
        elif(traceback[i][j] == "U"):
            #Gap up, write left
            topSequence = "-"+topSequence 
            leftSequence = traceback[i][0]+leftSequence
            i -= 1
            
        elif(traceback[i][j] == "L"):
            #Gap Left, write top
            leftSequence = "-"+leftSequence
            topSequence = traceback[0][j]+topSequence 
            j -= 1 

    #topSequence is the first sequence, topSequnce is the second sequence
    alignedSequences = [str(topSequence), str(leftSequence)]

    return alignedSequences

def Task2(sortedSequenceList, genome, labelToSeq, seqToLabe):
    #add first sequence as repSequence sortedSequenceList[0]
    clusterList = {sortedSequenceList[0] : []}
    clusterList = ClusterSequences(sortedSequenceList[0], sortedSequenceList, clusterList)
    clusterList = outSequenceCluster(clusterList, seqToLabel)
    print(clusterList)
    
    return clusterList

##Task 3 - Chimera Detection and Filtering
"""
The task here is to remove chimeric sequences from the representativeSequences list.
I am using the function filterChimeras to detect and filter out chimeras.
it returns a list of non-chimeric and chimeric sequence clusters(with representativeSequences as a key)
it does the following:
    1. Takes in list of representativeSequences from task 2
    2. Perform All to All alignment, SmithWaterman, For each repSequence in the list
    3. If an alignment is at leasat 60 bases in length, and scores >= 50, the sequence is added to the 
    good alignments for it's respective cluster. it's a bad alignment otherwise
    4. After determining the goodness of the alignment scores for every representativeSequence, 
    for each representativeSequence, determine if the representativeSequence has fewer good alignemnts
    in it's cluster than every other representativeSequence in theirs. If so, it's chimeric
    If the representativeSequence doesn't have the fewest good alignments, the representative sequence is
    non-chimeric
    5. Put chimeric and non-chimeric sequences into their respective lists
    
"""


def BuildAlignmentMatrixSW(sequencesToAlign):
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    alignMatrix = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the alignMatrix with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        alignMatrix.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(alignMatrix[0])):
        alignMatrix[0][i] = topSequenceRow[i-2]

    #FWrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in alignMatrix[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal += 1 


    #filling in the necessary null values 0 for row 2 and beyond in alignMatrix
    for j in range(2,len(alignMatrix[0])):
        alignMatrix[1][j] = 0

    #filling in the necessary values null values 0 for col 2 abd beyond in alignMatrix
    for row in alignMatrix[1:]:
        row[1] = 0

    return alignMatrix


def BuildTracebackMatrixSW(sequencesToAlign):
    #Now I need to build a scoring matrix
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    traceback = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the traceback with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        traceback.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(traceback[0])):
        traceback[0][i] = topSequenceRow[i-2]

    #Wrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in traceback[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal +=1

    return traceback

#Helper function for smithwaterman alg to find alignment score.           
def CalculateMaxScoreSW(alignMatrix, rowIndex, colIndex, gap_penalty, score_matrix):
    #find topLetter and leftLetter so we can find match/mismatch value in score_matrix
    topLetter = alignMatrix[0][colIndex]
    leftLetter = alignMatrix[rowIndex][0]

    scoreIndexTop = score_matrix[0].index(topLetter)
    scoreIndexLeft = score_matrix[0].index(leftLetter)

    #find score value at the score indexes found
    matchValue = int(score_matrix[scoreIndexTop][scoreIndexLeft])

    #Calculate score at the cell in alignMatrix
    diagVal = matchValue + alignMatrix[rowIndex-1][colIndex-1]
    upVal = alignMatrix[rowIndex-1][colIndex] + gap_penalty
    leftVal = alignMatrix[rowIndex][colIndex-1] + gap_penalty

    cellScore = max(diagVal, upVal, leftVal)
    #gap penalty found in scoring file chosen. For now just focus on the nucleatide scoring file (score_matrix)
    direction = ""
    if(cellScore == diagVal):
        direction = "D"
    elif(cellScore == upVal):
        direction = "U"
    elif(cellScore == leftVal):
        direction = "L"

    #Smith Waterman condition
    if(cellScore <= 0):
        cellScore = 0
        direction = 0 #0 is the equivalent of x in smithwaterman

    return cellScore, direction
#Helper function for smithwaterman alg to find alignment alignment between sequences
def AlignSequenceSW(traceback, startRow, startCol):
    #Finished stepping through matrix when 0 is reached, starting point is the last occurence of the alignment score
    leftSequence = ""
    topSequence = ""
    i = startRow
    j = startCol

    while traceback[i][j] != 0:
        if(traceback[i][j] == "D"):
            #Write both top and left into sequence listtopSequece
            leftSequence = traceback[i][0]+leftSequence
            topSequence = traceback[0][j]+topSequence 
            i -= 1
            j -= 1
        elif(traceback[i][j] == "U"):
            #Gap top, write left
            topSequence = "-"+topSequence 
            leftSequence = traceback[i][0]+leftSequence
            i -= 1
            
        elif(traceback[i][j] == "L"):
            #Gap Left, write top
            leftSequence = "-"+leftSequence
            topSequence = traceback[0][j]+topSequence 
            j -= 1 
        elif(traceback[i][j] == "0"):
            i = 1
            j = 1

    #topSequence is the first sequence, topSequnce is the second sequence
    alignedSequences = [str(topSequence), str(leftSequence)]

    return alignedSequences

def smithwaterman_algSW(sequencesToAlign):
    alignMatrix = BuildAlignmentMatrixSW(sequencesToAlign)
    traceback = BuildTracebackMatrixSW(sequencesToAlign)
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    cellScoreList = []
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = CalculateMaxScoreSW(alignMatrix, i, j, gapPenalty, score_matrix)
            cellScoreList.append(cellScore)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction

    alignmentScore = max(cellScoreList)
    #find the index of the last occurrence of the cell score to determine where to start in the trace back.
    for rowIndex in range(len(alignMatrix)):  #iterate over row indices
        for colIndex in range(len(alignMatrix[rowIndex])):  #iterate over column indices
            if alignMatrix[rowIndex][colIndex] == alignmentScore:  
                traceRow = rowIndex
                traceCol = colIndex


    alignedSequences = AlignSequenceSW(traceback, traceRow, traceCol) #traceRow and traceCol are starting positions for the traceback
    
    #cellScore is the max value in the alignedMatrix
            
    return alignedSequences, alignmentScore

def alignToAll(seq, repSequences):
    alignmentsToAdd = []
    for item in repSequences:
        if item == seq:
            continue
        else:
            alignedSequences, alignmentScore = smithwaterman_algSW([seq, item])
            if len(alignedSequences[0])>=60 and alignmentScore>=50:
                alignmentsToAdd.append(item)
            else:
                #bad alignment
                continue 



    return alignmentsToAdd

def determineChimeric(seq, repSequences, goodAlignments, clusterList):
    #seq is chimeric if every item in repSequences has more goodAlignments than seq    
    chimericStatus = True
    for item in goodAlignments[seq]:
        if len(clusterList[seq])>len(clusterList[item]):
            #more good alignments in item than seq, so seq is not a chimeric Sequence
            chimericStatus = False
            return chimericStatus

    return chimericStatus
    

def filterChimeras(sequences, clusterList):
    #Perform All to All alignment, SmithWaterman, For each repSequence in the list
    goodAlignments = {}
    for seq in sequences:
        alignmentsToAdd = alignToAll(seq, sequences)
        goodAlignments[seq] = alignmentsToAdd
    
    #determine if the representativeSequence is chimeric, ie it has fewer good alignemnts in it's cluster than every other representativeSequence in theirs. If so, it's chimeric
    chimericSeqs = []
    nonChimericSeqs = []
    for seq in sequences:
        chimericStatus = determineChimeric(seq, sequences, goodAlignments, clusterList)
        if chimericStatus == True:
            chimericSeqs.append(seq)
        elif chimericStatus == False:
            nonChimericSeqs.append(seq)
    print(chimericSeqs, nonChimericSeqs)
    
    return chimericSeqs, nonChimericSeqs

def outChimera(sequences, chimericSequences, seqToLabel, clusterList):
    #write sorted_sequence_data to file, but only the sequenceLabel in each line. omit the length.
    if(chimericSequences == True):
        fileName = f"{args.output}.chimeric.fna"
    if(chimericSequences == False):
        fileName = f"{args.output}.nonchimeric.fna"
        
    with open(fileName, "w") as output_file:
        for repSequence in sequences:
            title = seqToLabel[repSequence]
            size = len(clusterList[repSequence])
            output_file.write(f"{title}; size={size}\n{repSequence}\n")
            
    return 1

def Task3(clusterList, seqToLabel):
    repSequences =[]
    for item in clusterList:
        repSequences.append(item)
    print(repSequences)
    chimericSeqs, nonChimericSeqs = filterChimeras(repSequences, clusterList)
    for item in chimericSeqs:
        chimericSequences = True
        outChimera(chimericSeqs, chimericSequences, seqToLabel, clusterList)
    for item in nonChimericSeqs:
        chimericSequences = False
        outChimera(nonChimericSeqs, chimericSequences, seqToLabel, clusterList)

    return chimericSeqs, nonChimericSeqs
##Task 4 - Multiple sequence Alignment
"""

This task will receive the non-chimeric representative sequences from Task #3 as well as the scoring matrix given
via the command line argument ‘-s’ which will define the substitution scores that should be utilized by the NW
algorithm.
This task will produce the intermediary file <output_file>.msa.fna where <output_file> is the path and filename
given by the command line argument ‘-o’.
"""

###main implementation###

#main implementation
print("Final Project :: R#11583588")
parser = argparse.ArgumentParser(description="sorting Genome Sequences in descending order")
parser.add_argument("-i", "--input", required=True, type=str, help="File path to input.fna")
parser.add_argument("-o", "--output", required=True, type=str, help="File path for output.fna")
parser.add_argument("-s", "--score_matrix", required=True, type=str, help="File path for scoring matrix")
args = parser.parse_args()

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

#Execute task_1 function, return's sorted sequences
sortedSequenceList, labelToSeq, seqToLabel, sorted_sequence_data, genome = Task_1(args)
clusterList = Task2(sortedSequenceList, genome, labelToSeq, seqToLabel)
Task3(clusterList, seqToLabel)


#Execute task_2 function, return's clustered sequences and representative sequence list

#Execute task_3 Filter out chimeric sequeces from representative sequence list

#Execute task_4 Perform MSA of representative sequence list.