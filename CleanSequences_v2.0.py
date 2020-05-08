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

#Write statistics.
print('Number of sequences with errors:\t' + str(error))
print('Number of sequences without errors:\t' + str(i - error))
print('Total number of sequences checked:\t' + str(i))
print('Percent error:\t\t\t\t' + str(round((error/i) * 100,2)) + '%')
reff.close()

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)

"""
PURPOSE: To remove any sequences that contain residues not in the primer design of this experiment. 

SampleOutput (from 20-03-21):
Number of sequences with errors:	3021748
Number of sequences without errors:	15499047
Total number of sequences checked:	18520795
Percent error:				16.32%


"""
