#!/usr/bin/env python

#modules
import sys
import os
#make sure you have installed biopython to import Bio. Use: pip install biopython or pip3 install biopython
import Bio
from Bio import SeqIO
import numpy as np

#usage
#seq1 is the sequence 1
#seq2 is sequence 2
#output.txt is where you save the data

#python3 needlemanWunsch.py seq1.fna seq2.fna output.txt


#load files from terminal
def sysfiles(files):
    for input in files:
        print ("Loaded %s" % input)
    with open(sys.argv[1], 'r') as file, open(sys.argv[2], 'r') as file2:
        f1 = file.read()
        f2 = file2.read()
sysfiles(sys.argv[1:3]) #only call the first 2 files

outputFile = open(sys.argv[3], 'w') #this is to create the output file

#read seqs
def openSeq(file):
    sequence = SeqIO.parse(file, 'fasta')
    for info in sequence:
        secInfo = info
        return(secInfo)

#file 1
print("\n")
seq_1 = sys.argv[1]
seq1Info = openSeq(seq_1)
seq1len = int(len(seq1Info.seq))
print('Sequence 1 length %s: ' %seq1Info.id + str(seq1len))
lower = int(input("Insert lower value (e.g., 0): "))
upper = int(input("Insert upper value (e.g., length of sequence 1): "))
seq1 = seq1Info.seq[lower:upper].replace(r'\r\n', '') #removing new lines from the sequences
seq1 = seq1Info.seq[lower:upper].replace(r'\r\n', '')
#print(seq1)
print('Clean sequence')
#print("\n")


#file 2
seq_2 = sys.argv[2]
seq2Info = openSeq(seq_2)
seq2len = int(len(seq2Info.seq))
print('Sequence 2 length %s: ' %seq2Info.id + str(seq2len))
lower2 = int(input("Insert lower value (e.g., 0): "))
upper2 = int(input("Insert upper value (e.g., length of sequence 2): "))
seq2 = seq2Info.seq[lower2:upper2].replace(r'\r\n', '')  #removing \n from sequences
seq2 = seq2Info.seq[lower2:upper2].replace(r'\r\n', '')
#print(seq2)
print('Clean sequence')
print("\n")


print('Working with alignments')
#scores

#
initialGap = int(input("Initial gap penalty: "))
matchr = int(input("Match score: "))
missmatchp = int(input("Mismatch penalty: "))
gap = int(input("Fixed gap penalty: "))

#seq alignment with initial gap penalty
def needlemanW(sequence1, sequence2, initialGap, matchScore, mismatchPenalty, gappenalty):  
    #make arrays, fill and initialize them
    s = np.zeros((len(sequence1)+1,len(sequence2)+1))
    matching = np.zeros((len(sequence1),len(sequence2)))
    for i in range(len(sequence1)):
        for j in range(len(sequence2)):
            if sequence1[i] == sequence2[j]:
                matching[i][j]= matchScore
            else:
                matching[i][j]= mismatchPenalty
    for i in range(len(sequence1)+1):
        s[i][0] = i*gappenalty #vertical
        for j in range(len(sequence2)+1):
            s[0][j] = j*gappenalty
    for i in range(1,len(sequence1)+1):
        for j in range(1,len(sequence2)+1):
            s_0 = s[i-1][j-1]+matching[i-1][j-1]
            s_1 = s[i-1][j]+gappenalty
            s_2 = s[i][j-1]+gappenalty
            s[i][j] = max(initialGap, s_0, s_1, s_2)
    seq1_al = ''
    seq2_al = ''
    l1 = len(sequence1)
    l2 = len(sequence2)
    out_score = s[l1][l2]    
    while(l1 >0 or l2 > 0):
        if (l1 >0 and l2 > 0 and s[l1][l2] == s[l1-1][l2-1]+ matching[l1-1][l2-1]):
            seq1_al = sequence1[l1-1] + seq1_al
            seq2_al = sequence2[l2-1] + seq2_al
            l1 -= 1 
            l2 -= 1
        elif(l1 > 0 and s[l1][l2] == s[l1-1][l2] + gappenalty):
            seq1_al = sequence1[l1-1] + seq1_al
            seq2_al = '_' + seq2_al
            l1 -= 1
        else:
            seq1_al = '_' + seq1_al
            seq2_al = sequence2[l2-1] + seq2_al
            l2 -= 1
    centerRow = ''  #this is the row that will contain the | for match and x for mismatch
    for i in range(len(seq1_al)):
        if seq1_al[i] == seq2_al[i]:
            centerRow += '|'
        elif seq1_al[i] != seq2_al[i]:
            if (seq1_al[i] == '_' or seq2_al[i] == '_'):
                centerRow += ' '
            else:
                centerRow += 'x'
    numberofMatches = 0   #getting the number of matches
    for i in range(len(seq1_al)):
        if seq1_al[i] == seq2_al[i]:
            numberofMatches = numberofMatches + 1
    numberofMismatches = 0  #counting mismatches
    for i in range(len(seq1_al)):
        if seq1_al[i] != seq2_al[i]:
            numberofMismatches = numberofMismatches + 1
    print(out_score,  file = outputFile)
    print('>%s' % seq1Info.id, file = outputFile)
    print (seq1_al,  file = outputFile)
    print(centerRow,  file = outputFile)
    print(seq2_al,  file = outputFile)
    print('>%s' % seq2Info.id, file = outputFile)
    print('Number of matches using initial gap: ', numberofMatches,  file = outputFile)
    print('Number of mismatches: ', numberofMismatches, file = outputFile)
    print('Sequences above were aligned using initial gap penalty', file = outputFile)
    print('\n\n', file = outputFile)

#Add code to save file
needlemanW(seq1, seq2, initialGap, matchr, missmatchp, gap)

print('\n')
print('DONE')
