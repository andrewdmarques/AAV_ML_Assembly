#!/usr/bin/env python3

import sys
import datetime
import re

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\nPython computations running...\n" % timer0)

#Setup files.
ref = sys.argv[1]
seqfile = sys.argv[2]
reff = open(ref, 'r')
frem = open(seqfile[:-4] + "_removed.txt", 'w+')
frem.close()
logfile = "log_" + seqfile[:-4] + "_cleaned.txt"
logf = open(logfile, "w+")
logf.write("Start time:			%s\n" % timer0)

#Generate empty variables.
error = 0
pos = 0
nlr = 1
nlc = 1

#Gather sampling sequence data.
with open(seqfile, 'r') as seqf:
	i = 0
	for seq in seqf:
		i += 1

#Define function for removing sequences with errors (sweeping).
def sweep(pos, possibles, file_in, file_out, file_remove, nlr, nlc, error):
	fout = open(file_out, 'w+')
	frem = open(file_remove, 'a')
	with open(file_in) as fin:
		for seq in fin:
			seq = seq.strip()
			if seq[pos] not in possibles:
				error +=1
				if nlr == 1:
					frem.write(seq)
					nlr = 0
				else:
					frem.write('\n' + seq)
			else:
				if nlc == 1:
					fout.write(seq)
					nlc = 0
				else:
					fout.write('\n' + seq)
	fout.close()
	frem.close()
	return error, nlr
	
#Produce list of all possible residues.
for pos_res in reff:
	possibles = []
	pos_res = pos_res.strip()
	for res in pos_res:
		possibles.append(res)
	if pos == 0:
		file_in = sys.argv[2]
	else:
		file_in = seqfile[:-4] + "_cleaned-" + str(pos) + ".txt"
	file_out = seqfile[:-4] + "_cleaned-" + str(pos+1) + ".txt"
	file_remove = seqfile[:-4] + "_removed.txt"
	error, nlr = sweep(pos, possibles, file_in, file_out, file_remove, nlr, nlc, error)
	print('Completed filtering for position ' + str(pos + 1))
	pos += 1
reff.close()

#*****************************Cleaned2Sort*****************************

sys.stderr.write("Removing duplicates...\n\n")

#Open the input file.
filename = file_out
#print('File to be sorted:\t' + file_out)
#Sort the infile.
logs = open(file_out[:-4] + "_sorted.txt", "w+")
def build_index(filename, sort_col):
    index = []
    f = open(filename)
    while True:
        offset = f.tell()
        line = f.readline()
        if not line:
            break
        length = len(line)
        col = line.split('\t')[sort_col].strip()
        index.append((col, offset, length))
    f.close()
    index.sort()
    return index

def print_sorted(filename, col_sort):
    index = build_index(filename, col_sort)
    f = open(filename)
    for col, offset, length in index:
        f.seek(offset)
        sortedseqs = f.read(length).rstrip('\n')
        logs.write(sortedseqs + "\n")

if __name__ == '__main__':
    sort_col = 0
    print_sorted(filename, sort_col)

#Close the infile and the sorted file.
logs.close()

#*****************************Sort2RemoveDuplicates*****************************

#Compresses the matching sequences.
log_unnamed = open(file_out[:-4] + "_sorted.txt", "r")
log_compressed = open(file_out[:-4] + "_no-duplicates.txt", "w+")
previous = ""
repnum = 1
dup_seq = 0
for line in log_unnamed:
	current = line
	if current == previous:
		repnum += 1
		dup_seq += 1
	elif previous != "":
		spaces = ' '*(16-len(str(repnum)))
		log_compressed.write(previous)
		repnum = 1
	previous = current
log_unnamed.close()
log_compressed.close()

#Write statistics.
print('Total number of sequences checked:\t' + str(i))
print('Number of sequences without errors:\t' + str(i - error))
print('Number of sequences with errors:\t' + str(error))
print('Percent error:\t\t\t\t' + str(round((error/i) * 100,2)) + '%')
print("Number of duplicate sequences:\t\t" + str(dup_seq))
print("Percent of duplicate sequences:\t\t" + str(dup_seq/i * 100) + "%\n")
logf.write('\nTotal number of sequences checked:\t' + str(i))
logf.write('\nNumber of sequences without errors:\t' + str(i - error))
logf.write('\nNumber of sequences with errors:\t' + str(error))
logf.write('\nPercent error:\t\t\t\t' + str(round((error/i) * 100,2)) + '%')
logf.write("\nNumber of duplicate sequences:\t" + str(dup_seq))
logf.write("\nPercent of duplicate sequences:\t" + str(dup_seq/i * 100) + "%")

#Report the end time.
end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)
logf.write("\n\nPython Completed:		%s\n\n" % timer_0)
logf.close()

"""
PURPOSE: To remove any sequences that contain residues not in the primer design of this experiment. 

SampleOutput (from 20-03-21):
Number of sequences with errors:	3021748
Number of sequences without errors:	15499047
Total number of sequences checked:	18520795
Percent error:				16.32%


"""
