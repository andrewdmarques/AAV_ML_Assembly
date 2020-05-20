#!/usr/bin/env python3

import sys
import datetime
import re
import random

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\nPython computations running...\n" % timer0)

#Parameters
num_seqs = 100000
wt_pos = 33
#Open files
infile = sys.argv[1]
outfile = str(num_seqs/1000) + 'k-seqs_wt-pos-' + str(wt_pos) + ".txt"
outf = open(outfile, 'w+')


wt = 'TPSTTQRQKTNNNTKKEEKQGSEKTNVDIEKRR'

#Define the function to generate the randomized sequence.
def rand_seq(wt_pos):
	with open(infile) as f:
		seq = ''
		i = 0
		for possibles in f:
			if i == wt_pos - 1:
				res = wt[i]
			else:
				num_possibles = len(possibles) - 1 #-2 for the '\n' and zero indexing
				rand_pos = random.randrange(0,num_possibles)
				res = possibles[rand_pos]
			seq += res
			i += 1
		return seq

#Generate x number of randomizes sequences.
for i in range(0,num_seqs):
	seq = rand_seq(wt_pos)
	if i == 0:
		outf.write(seq)		
	else:
		outf.write('\n' + seq)	
	
end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)


"""
PURPOSE: The purpose of this script is to read the CL8_possibles script and generate a specified number of sequences containing the desired mutations.

PROCEDURE:
1.
2.
3.
4.




Sample input file:


Sample output file:

"""
