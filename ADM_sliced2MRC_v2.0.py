#!/usr/bin/env python3

import sys
import datetime
import re
import os
import string
import random

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\nPython computations running...\n" % timer0)

inseq = sys.argv[1]

#Assign the sequences from the input file.
x = 0
with open(inseq) as f:
	for line in f:
		if x == 1:
			nt_cap = line
			x = 0
		if '>part2' in line:
			x = 1
x = 0
with open(inseq) as f:
	for line in f:
		if x == 1:
			aa_mutate = line
			x = 0
		if '>part3' in line:
			x = 1

x = 0
with open(inseq) as f:
	for line in f:
		if x == 1:
			aa_template = line
			x = 0
		if '>part4' in line:
			x = 1

#Determine the start of the nt sequence.
start = 0
for nt in nt_cap:
	if nt.isupper():
		break
	start += 1
nt1 = nt_cap[start-15:start-7].upper()

#Determine the end position of the nt sequence.
end = 0
p = 0
for nt in nt_cap:
	if nt.isupper():
		nt_end = p - start + 16
	p += 1
		
#Define the function for determining the beginning of each variable region.
start_vr = 0
detect_vr = 0
end_vr = 0
final_vr = 0
p = 0
previous = 'not X'
vr_pos = []
vr_len = []
for res in aa_mutate:
	#If residue is not a variable residue.
	if res != 'X':
		start_vr += 1
		if detect_vr == 1:
			end_vr += 1
		if previous == 'X':
			final_vr = 0
		final_vr += 1
	#If beginning of new variable region.
	if start_vr >= 9 and res == 'X':
		#print('Position ' + str(p))
		vr_pos.append(p)
		start_vr = 0
		end_vr = 0
	#If variable residue in variable region.
	if res == 'X':
		start_vr = 0
		end_vr += 1
		detect_vr = 1
	#If sequence carries on for 9 residues without a variable residue, end of variable region is detected.
	if start_vr == 9:
		end_vr_record = end_vr
		if detect_vr == 1:
			vr_len.append(end_vr_record - 9)
	#If the final variable region in the sequence.
	if res == '\n':
		vr_len.append(end_vr - final_vr)
	previous = res
	p += 1

#Determine the number of variable regions
vr_num = len(vr_pos)
mutate_vr = []
template_vr = []
query_vr = []


#Assemble template residues and query residues.		
for i in range(0,vr_num):
	start_vr = vr_pos[i]
	end_vr = vr_pos[i] + vr_len[i]
	mutate_vr.append(aa_mutate[start_vr:end_vr])
	template_vr.append(aa_template[start_vr:end_vr])
	query_vr.append(aa_template[start_vr - 5:start_vr])

#Set the log names for the remainder of the script.
now = datetime.datetime.now()
time = now.strftime("%Y-%m-%d_%H-%M-%S")
filelocation = "./Output-MRC_" + time
os.system("mkdir " + filelocation)
logfile = filelocation + "/temp_" + time
seqsfile = filelocation + "/" + time
log = open(filelocation + "/log_" + time + ".txt", "w+")
log.write("Beginning ADM_MRC script.\n\n")
log.write("Start time:			%s\n\n" % timer0)
seqsin = 0
seqsout = 0

#Determine which residues were mutated.
mrs = ''
for i in range(0,vr_num):
	mrs += mutate_vr[i] + ' '*5
mrs += '\n'

#*****************************sliced2Trimmed*****************************
sys.stderr.write("\nTrimming computations running...\n")
log.write("Trimming computations running...\n")

sliced_file = sys.argv[2]
trim_file = logfile + '_01-trimmed.txt'

trim = open(trim_file, "w+")
nl = 0
with open(sliced_file) as f:
	for line in f:
		if nl == 1:
			trim.write('\n')
		if nl == 0:
			nl = 1
		line = line.strip()
		for i in range(0,len(line)):
			if mrs[i] == 'X':
				trim.write(line[i])

trim.close()

#Record the end time.
end_trim = datetime.datetime.now()
timer_trim = str(end_trim)
sys.stderr.write("Trimming Completed:		%s\n\n" % timer_trim)
log.write("Trimming Completed:		%s\n\n" % timer_trim)

#*****************************Trimmed2Filled*****************************
sys.stderr.write("Filling computations running...\n")
log.write("Filling computations running...\n")

trim_file = logfile + '_01-trimmed.txt'
fill_file = logfile + '_02-filled.txt'

#Determine what the wildtype for each variable residue is.
vres_template = ''
for i in range(0, len(aa_mutate)):
	if aa_mutate[i] == 'X':
		vres_template += aa_template[i]

#Reinsert the template residues into the sequences in place of '.'
fill = open(fill_file, "w+")
j = 0
with open(trim_file) as f:
	for line in f:
		line = line.strip()
		if len(line) == len(vres_template):
			if j != 0:
				fill.write('\n')
			for i in range(0, len(vres_template)):
				if line[i] == '.':
					fill.write(vres_template[i])
				else:
					fill.write(line[i])
			j += 1
fill.close()

#Record the end time.
end_fill = datetime.datetime.now()
timer_fill = str(end_fill)
sys.stderr.write("Filling Completed:		%s\n\n" % timer_fill)
log.write("Filling Completed:		%s\n\n" % timer_fill)

#*****************************Filled2Sorted*****************************

sys.stderr.write("Sorting computations running...\n")
log.write("Sorting computations running...\n")

#Open the input file.
filename = logfile + "_02-filled.txt"

#Sort the infile.
logs = open(logfile + "_03-sorted.txt", "w+")
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

#Record the end time.
end_grp = datetime.datetime.now()
timer_grp = str(end_grp)
sys.stderr.write("Sorting Completed:		%s\n\n" % timer_grp)
log.write("Sorting Completed:		%s\n\n" % timer_grp)

#*****************************Sorted2Compressed*****************************

sys.stderr.write("Compressing computations running...\n")
log.write("Compressing computations running...\n")

#Compresses the matching sequences.
log_unnamed = open(logfile + "_03-sorted.txt", "r")
log_compressed = open(seqsfile + "_compressed.txt", "w+")
previous = ""
nl = 0
repnum = 1
for line in log_unnamed:
	current = line.strip()
	if current == previous:
		repnum += 1
	elif previous != "":
		if nl == 1:
			log_compressed.write('\n')
		if nl == 0:
			nl = 1
		log_compressed.write(previous + '\t' + str(repnum))
		repnum = 1
	previous = current
log_unnamed.close()
log_compressed.close()

#Record the end time.
end_comp = datetime.datetime.now()
timer_comp = str(end_comp)
sys.stderr.write("Compressing Completed:		%s\n\n" % timer_comp)
log.write("Compressing Completed:		%s\n\n" % timer_comp)


#compute sequence statistics and complete log.
end = datetime.datetime.now()
timer_end = str(end)
log.write("\nADM_MRC Completed:	%s\n\n" % timer_end)

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("\nPython Completed:		%s\n\n" % timer_0)

log.close()

"""
PURPOSE:

PROCEDURE:
1. To sort by count with highest first, use "sort -k 2 -r compressed.txt > compressedcounts.txt"




Sample input file:


Sample output file:

"""
