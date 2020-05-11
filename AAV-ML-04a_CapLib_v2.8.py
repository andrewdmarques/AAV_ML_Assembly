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
filelocation = "./Output-Caplib_" + time
os.system("mkdir " + filelocation)
logfile = filelocation + "/temp_" + time
seqsfile = filelocation + "/" + time
log = open(filelocation + "/log_" + time + ".txt", "w+")
log.write("Beginning ADM_CapLib script.\n\n")
log.write("Start time:			%s\n\n" % timer0)
seqsin = 0
seqsout = 0

#*****************************FASTQ2FASTA****************************
sys.stderr.write("Fastq2Fasta computations running...\n")
log.write("Fastq2Fasta computations running...\n")

fq = sys.argv[2]
fa = open(logfile + "_01-fasta.txt", "w+")

location = ""
i = 0

tab = str.maketrans("ACTG", "TGAC")

def reverse_complement(seq):    
    return seq.translate(tab)[::-1]

with open(fq) as f:
	for line in f:
		if line[0] == "@":
			location = line[28:len(line) - 25]
			i = 1
		elif i == 1:
			fa.write( ">" + location + "\n")

			s1 = re.search(r'(?<=' + nt1 + ')\w+', line)
			if s1:				
				fa.write(line + "")
			else:
				line_rev = reverse_complement(line)
				fa.write(line_rev.strip() + "\n")
			i = 0
			location = ""
fa.close()

#Record the end time.
end_fq = datetime.datetime.now()
timer_fq = str(end_fq)
sys.stderr.write("Fastq2Fasta Completed:		%s\n\n" % timer_fq)
log.write("Fastq2Fasta Completed:		%s\n\n" % timer_fq)

#*****************************NT2AAFASTA*****************************

sys.stderr.write("Translating computations running...\n")
log.write("Translating computations running...\n")

#Open log for translating nt to AA

ntfile = logfile + "_01-fasta.txt"
logAA = open(logfile + "_02-translated.txt", "w+")

table = { 
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ATN':'X',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'ACN':'X',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AAN':'X',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'AGN':'X',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CTN':'X',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CCN':'X',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CAN':'X',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CGN':'X',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GTN':'X',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GCN':'X',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GAN':'X',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GGN':'X',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTN':'X',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TCN':'X',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TAN':'X',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 'TGN':'X',
        'NCA':'X', 'NCC':'X', 'NCG':'X', 'NCT':'X', 'NTN':'X',
        'NTC':'X', 'NTT':'X', 'NTA':'X', 'NTG':'X', 'NCN':'X',
        'NAC':'X', 'NAT':'X', 'NAA':'X', 'NAG':'X', 'NAN':'X',
        'NGC':'X', 'NGT':'X', 'NGA':'X', 'NGG':'X', 'NGN':'X',
        'ANA':'X', 'ANC':'X', 'ANG':'X', 'ANT':'X', 'ANN':'X',
        'CNA':'X', 'CNC':'X', 'CNG':'X', 'CNT':'X', 'CNN':'X',
        'GNA':'X', 'GNC':'X', 'GNG':'X', 'GNT':'X', 'GNN':'X',
        'TNA':'X', 'TNC':'X', 'TNG':'X', 'TNT':'X', 'TNN':'X',
        'NNA':'X', 'NNC':'X', 'NNG':'X', 'NNT':'X', 'NNN':'X',
    } 

with open(ntfile) as f:
	for line in f:
		if line[0] == ">":
			seqsin += 1
		#Retrieve region to be translated.	
		s1 = re.search(r'(?<=' + nt1 + ')\w+', line)
		if s1:			
			tempseq = nt1 + s1.group(0)
			ntseq = tempseq[0:nt_end]
			if len(ntseq)%3 == 0: 
				for i in range(0, len(ntseq), 3): 
					codon = ntseq[i:i + 3] 
					protein += table[codon] 
				logAA.write(protein + "\n")  
 
		protein =""
logAA.close() 

#Record the end time.
end_nt = datetime.datetime.now()
timer_nt = str(end_nt)
sys.stderr.write("Translating Completed:		%s\n\n" % timer_nt)
log.write("Translating Completed:		%s\n\n" % timer_nt)

#*****************************Fasta2VRSlices*****************************

sys.stderr.write("Slicing computations running...\n")
log.write("Slicing computations running...\n")

#Open the log file for the VR slices.
logslice = open(logfile + "_03-sliced.txt", "w+")

#Open the file to be sliced.
infileslice = logfile + "_02-translated.txt"

#Start timer.
start = datetime.datetime.now()
timer0 = str(start)

#Search for desired VRs within the sequence.
with open(infileslice) as f:
	vr_retained = []
	num_seq = 0
	num_good_seq = 0
	for line in f:
		num_seq += 1
		vr_list = []
		for i in range(0,vr_num):
			vr_string = ''
			q = query_vr[i]
			m = re.search(r'(?<=' + q + ')\w+', line)
			t = template_vr[i]
			l = len(t)
			if m:
				tempseq = m.group(0)
				VR = tempseq[0:l]
				#Removing consensus residues			
				n = 0
				for res in VR:
					if res == t[n]:
						vr_string += '.'
					else:
						vr_string += res
					n = n + 1
			#If the VR in tact, record the VR.
			if vr_string != '':
				vr_list.append(vr_string)
				vr_retained += str(i + 1)
		#If the sequence has all the VRs in tact, record the VRs.
		if len(vr_list) == vr_num:
			num_good_seq += 1
			for i in range(0,vr_num):
				if i < vr_num - 1:
					logslice.write(vr_list[i])
					logslice.write('     ')
				if i == vr_num - 1:
					logslice.write(vr_list[i])
					logslice.write('\n')
logslice.close()
vr_retained.sort()
#Record the end time.
end_sli = datetime.datetime.now()
timer_sli = str(end_sli)
sys.stderr.write("Slicing Completed:		%s\n\n" % timer_sli)
log.write("Slicing Completed:		%s\n\n" % timer_sli)						

#*****************************VRSlices2Group*****************************

sys.stderr.write("Grouping computations running...\n")
log.write("Grouping computations running...\n")

#Open the input file.
filename = logfile + "_03-sliced.txt"

#Sort the infile.
logs = open(logfile + "_04-sorted.txt", "w+")
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
sys.stderr.write("Grouping Completed:		%s\n\n" % timer_grp)
log.write("Grouping Completed:		%s\n\n" % timer_grp)

#*****************************Grouped2Compressed*****************************

sys.stderr.write("Compressing computations running...\n")
log.write("Compressing computations running...\n")

#Compresses the matching sequences.
log_unnamed = open(logfile + "_04-sorted.txt", "r")
log_compressed = open(logfile + "_05-compressed.txt", "w+")
previous = ""
repnum = 1
for line in log_unnamed:
	current = line
	if current == previous:
		repnum += 1
	elif previous != "":
		spaces = ' '*(16-len(str(repnum)))
		log_compressed.write(str(repnum) + spaces + previous)
		repnum = 1
	previous = current
log_unnamed.close()
log_compressed.close()

#Record the end time.
end_comp = datetime.datetime.now()
timer_comp = str(end_comp)
sys.stderr.write("Compressing Completed:		%s\n\n" % timer_comp)
log.write("Compressing Completed:		%s\n\n" % timer_comp)

#*****************************Compressed2Comsort*****************************

sys.stderr.write("Comsort computations running...\n")
log.write("Comsort computations running...\n")

compressed_infile = logfile + "_05-compressed.txt"
log_comsort = seqsfile + "_comsort.txt"

#Determine which residues were mutated.
mrs = 'Template:' + ' '*2
for i in range(0,vr_num):
	mrs += ' '*5 + mutate_vr[i]
mrs += '\n'

#Insert the header line for comsot file.
f = open(log_comsort, "w+")
f.write(mrs)
f.close()

#Use linux command to sort the sequences by copy number.
os.system('sort -n -r ' + compressed_infile + ' >> ' + log_comsort)

#Record the end time.
end = datetime.datetime.now()
timer_end = str(end)
sys.stderr.write("Comsort Completed:		%s\n\n" % timer_end)
log.write("Comsort Completed:		%s\n\n" % timer_end)

#*****************************Comsort2Trim*****************************
sys.stderr.write("Trimming computations running...\n")
log.write("Trimming computations running...\n")

comsort_file = seqsfile + '_comsort.txt'
trim_file = logfile + '_06-trimmed.txt'

trim = open(trim_file, "w+")
with open(comsort_file) as f:
	next(f)
	for line in f:
		for i in range(0,len(line)):
			if mrs[i] == 'X':
				trim.write(line[i])
		trim.write('\n')
trim.close()

#Record the end time.
end_trim = datetime.datetime.now()
timer_trim = str(end_trim)
sys.stderr.write("Trimming Completed:		%s\n\n" % timer_trim)
log.write("Trimming Completed:		%s\n\n" % timer_trim)

#*****************************Trim2pre-ML*****************************
sys.stderr.write("Pre-ML file computations running...\n")
log.write("Pre-ML file computations running...\n")

trim_file = logfile + '_06-trimmed.txt'
ml_file = seqsfile + '_pre-ml.txt'

#Determine what the wildtype for each variable residue is.
vres_template = ''
for i in range(0, len(aa_mutate)):
	if aa_mutate[i] == 'X':
		vres_template += aa_template[i]

#Reinsert the template residues into the sequences in place of '.'
pre_ml = open(ml_file, "w+")
with open(trim_file) as f:
	for line in f:
		if '_' not in line and len(line) - 1 == len(vres_template):
			for i in range(0, len(line)):
				if line[i] == '.' or line[i] == '_' or line[i] == 'X':
					pre_ml.write(vres_template[i])
				else:
					pre_ml.write(line[i])		
pre_ml.close()

#Record the end time.
end_ml = datetime.datetime.now()
timer_ml = str(end_ml)
sys.stderr.write("Pre-ML file Completed:		%s\n\n" % timer_ml)
log.write("Pre-ML file Completed:		%s\n\n" % timer_ml)

#*****************************pre-ML2pre-MLshuf*****************************
sys.stderr.write("Shuffling computations running...\n")
log.write("Shuffling computations running...\n")

pre_ml_file = seqsfile + '_pre-ml.txt'
shuf_file = seqsfile + '_pre-ml_shuffled.txt'

#Shuffle the pre-ml data
with open(pre_ml_file,'r') as source:
    data = [ (random.random(), line) for line in source ]
data.sort()
with open(shuf_file,'w') as target:
    for _, line in data:
        target.write( line )

#Record the end time.
end_ml = datetime.datetime.now()
timer_ml = str(end_ml)
sys.stderr.write("Shuffling Completed:		%s\n\n" % timer_ml)
log.write("Shuffling Completed:		%s\n\n" % timer_ml)


#compute sequence statistics and complete log.
end = datetime.datetime.now()
timer_end = str(end)
log.write("\nADM_CapLib_v3 Completed:	%s\n\n" % timer_end)
log.write("Number of sequences read in:\t\t" + str( seqsin ) + "\n")
seq_analysis = 'Number of sequences retained:\t\t' + str(num_good_seq) + '\t' + str(round(num_good_seq / num_seq * 100,2)) + '%\n' 
log.write(seq_analysis)
previous = ''
counter = 1
for i in range(0,len(vr_retained)):
	current = vr_retained[i] 
	if current == previous:
		counter += 1
	if previous != '':
		if current != previous or i == len(vr_retained) - 1:
			vr_analysis = 'Number of sequences containing VR' + str(previous) + ':\t' + str(counter) + '\t' + str(round(counter / num_seq * 100,2)) + '%\n'
			sys.stdout.write(vr_analysis)
			log.write(vr_analysis)
			counter = 1
	previous = current
sys.stdout.write(seq_analysis)

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("\nPython Completed:		%s\n\n" % timer_0)

#Write to log the sequence data.
log.write('\nNumber of VRs detected: %s\n\n' % (str(vr_num)))
vrs = ''
qrs = ''
for i in range(0,vr_num):
	vrs += str(i+1) + '.' + template_vr[i] + '\t'
	qrs += str(i+1) + '.' + query_vr[i] + '\t'
log.write('Variable regions:\n' + vrs + '\n')
log.write('\nNucleotide query sequence:\n' + nt1 + '\n')
log.write('\nAmmino acid query sequences:\n' + qrs + '\n')

log.close()

"""
PURPOSE: The purpose of this script is to go from merged fastq files to (1) a compressed format containing only the positions which were designed to have mutations nad (2) a readable format that the user can easily interpret.

PROCEDURE:
1. ./AAV-ML-04a_CapLib_v2.7.py ref-04a_sequence-templates.txt sample-04a_merged-reads.fastq
"""
