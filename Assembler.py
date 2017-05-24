#!/usr/bin/python

import argparse
def makeArgs():
    parser=argparse.ArgumentParser(
        description="This script takes a file in fasta format, searches for overlapping sequences and returns a unique sequence")
    parser.add_argument("-inputFile",required=True,
                        help="Input sequence file in fasta format")
    return parser

if __name__ == "__main__":
    arguments=makeArgs()
    arguments=arguments.parse_args()
    fasta=arguments.inputFile

data={}
Id=None

#############################
## Store reads in dictionary
#############################
with open(fasta,"r") as reads:
    for lines in reads:
        line=lines.strip()
        if line[0]==">":
            Id=line
            data[Id]=""
        else:
            data[Id]=data[Id] + line

######################################
## Find overlap between two sequences
######################################
def overlap(seq1,seq2):
    for i in range(len(seq1)):
        ###Strips i from left side of seq1 and right side of seq2 to search for pregressively less overlapping sequence
        if seq1[i:] == seq2[:len(seq1)-i]:
            return seq1[i:]
    return ''

#####################################
##Search reads for overlap 
#####################################
def total_overlap(data):
    overlap_dict={}
    for seq_Id1, read1 in data.items():
        for seq_Id2, read2 in data.items():
            if seq_Id1==seq_Id2:
                continue
            if seq_Id1 not in overlap_dict:
                overlap_dict[seq_Id1]={}
            ###Places length of overlap as value for seq2 which combined are the value for seq1
            overlap_dict[seq_Id1][seq_Id2]=len(overlap(read1,read2))
    return overlap_dict

####################################
##Search overlaps for first read
####################################
def findFirstRead(overlaps):
    for i in overlaps:
        signifOverlaps = False
        ###First read shouldn't significantly overlap anything as the right read
        for j in overlaps[i]:
            if overlaps[j][i] >= 5:
                signifOverlaps = True
        if not signifOverlaps:
            return i
                        
######################################
##Search for largest overlap value
######################################
def findReadLargestOverlap(overlaps):
    maxOverlap = max(overlaps.values())
    for Id in overlaps:
        if overlaps[Id] == maxOverlap:
            return Id

##############################
##Search for order of reads
##############################
def findOrder(first,overlaps):
    order=[first]
    counter=1
    while counter < len(overlaps.keys()):
        if len(order)==1:
            nextRead=findReadLargestOverlap(overlaps[first])
            order.append(nextRead)
            counter += 1
        else:
            nextRead=findReadLargestOverlap(overlaps[order[-1]])
            order.append(nextRead)
            counter += 1
    return order

#########################################
##Check if assembly contains all reads
#########################################
def checkAssembly(order,data):
    idTotal=len(data.keys())
    overlapId=0
    for Id in data.keys():
        if Id in order:
            overlapId += 1
    if idTotal == overlapId:
        return True
    else:
        return False

##############################################
##Assemble reads and output unique sequence
##############################################
def assembleGenome(readOrder, reads, overlaps):
    genome = ''
    for readName in readOrder[:-1]:
        ###Takes length of overlap and removes it from sequence then adds it to genome string
        rightOverlap = max(x for x in overlaps[readName].values() if x >= 5)
        genome += reads[readName][:-(rightOverlap)]
    genome += reads[readOrder[-1]]
    if checkAssembly(readOrder,reads):
        return genome
    else:
        print "Error not all reads contained in assembly"

overlapData=total_overlap(data)
order = findOrder(findFirstRead(overlapData), overlapData)
print assembleGenome(order, data , overlapData)

