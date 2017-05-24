This program takes in a set of reads in fasta format, searches
for overlapping sequence and returns a unique sequence. After reading
in the file the program stores the read IDs and sequences in a
dictionary. It then searches for overlapping sequence, testing a read
to the left and right of every other read, returning a new dictionary
containing the length of all overlap combinations. Next it parses the
overlap dictionary searching for the first read by looking for a read
with no significant overlap in the right position. After finding the
first read the program searches for the order that the reads line up
by looking for the next read with the largest overlap. After finding
the order of the reads, if the assembly order contains all the reads
it will then assembles the sequences by subtracting the overlapping
sequence from each read then adding on the final read to the end. To
test this program simply input a file in fasta format with the option
-inputFile <Path/to/file>. 

