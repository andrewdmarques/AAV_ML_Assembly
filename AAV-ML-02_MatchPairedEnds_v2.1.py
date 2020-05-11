#!/usr/bin/env python3
from __future__ import division
import sys
import os
import datetime
import re

#Set the output names.
now = datetime.datetime.now()
time = now.strftime("%Y-%m-%d_%H-%M-%S")
filelocation = "./Output-Matched_" + time
os.system("mkdir " + filelocation)
outfile = filelocation + "/seqs_" + time

start = datetime.datetime.now()
timer0 = str(start)
a = "\n\n\n\n\nStart time:			%s\n\nMatching computations running...\n" % timer0
sys.stderr.write(a)

seqsin_R1 = 0
seqsin_R2 = 0
matches = 0

#Add R1 and R2 to the combined file, transforming them to tab deliniated reads.
out1 = open(outfile + "_temp01-combined.tab", "w+")
R1 = open(sys.argv[1],"r")
R2 = open(sys.argv[2],"r")

linenum = 1
for line in R1:
	if linenum == 4:
		out1.write(line.strip() + "\n")
		linenum = 1
		seqsin_R1 += 1
	else:
		out1.write(line.strip() + "	")
		linenum += 1
linenum = 1
for line in R2:
	if linenum == 4:
		out1.write(line.strip() + "\n")
		linenum = 1
		seqsin_R2 += 1
	else:
		out1.write(line.strip() + "	")
		linenum += 1

out1.close()
R1.close()
R2.close()

#Sort the tab seperated files.
combined = outfile + "_temp01-combined.tab"
sortd = outfile + "_temp02-sortd.tab"

os.system('sort -n ' + combined + ' > ' + sortd)

#Break up the sorted files into two matching files
out2 = open(outfile + "_temp02-sortd.tab","r")
out3 = open(outfile + "_matchesR1.fastq", "w+")
out4 = open(outfile + "_matchesR2.fastq", "w+")

l = 0
turn = "A"
turnA = ""
turnB = ""
headB = ""

for line in out2:
	if turn == "A":
		for word in line:
			l += 1
			if word == " ":
				break
		seqA = line
		headA = line[0:l-1]
		turn = "B"
		l = 0
	elif turn == "B":
		for word in line:
			l += 1
			if word == " ":
				break
		seqB = line
		headB = line[0:l-1]
		turn = "A"
		l = 0
	if headA == headB:
		seqA_fq = seqA.replace("	", "\n")
		seqB_fq = seqB.replace("	", "\n")
		out3.write(seqA_fq)
		out4.write(seqB_fq)
		matches += 1

out2.close()
out3.close()
out4.close()

#Output Statistics
end1 = datetime.datetime.now()
timer1 = str(end1)
b = "Matching Completed:		%s\n\n" % timer1
c = "\nSequences in from R1:	" + str(seqsin_R1)
d = "\nSequences in from R2:	" + str(seqsin_R2)
e = "\nMatching sequences:	" + str(matches)
percent = round(((2 * matches) / (seqsin_R1 + seqsin_R2)) * 100,2)
f = "\nPercent of sequences with matches: " + str(percent) + "%"

logfile = filelocation + "/log_" + time
log = open(logfile + ".txt", "w+")
log.write(a + b + c + d + e + f)
log.close()

sys.stderr.write(b + c + d + e + f)
sys.stderr.write("\n\nYou may want to delete the residual files temp01 and temp02. Files matchesR1 and matchesR2 are the new R1 and R2 files.\n\n")

"""
PURPOSE:
The purpose of this script is to read in two paired end FASTQ files with unequal reads and output two files which can be used for merging (they are in the same order and each well location has a match).

PROCEDURE:
1. Combine R1 and R2 in a single file using a tab delineated format.
2. Sort the sequence headers in the file so that matching sequences are placed together.
3. Break the sequences into two matching files of R1 and R2.
4. Convert the .tab files back to .fastq files.

$ ./AAV-ML-02_MatchPairedEnds_v2.0.py sample-02_R1-ForwardReads_100-reads.fq sample-02_R2-ReverseReads_100-reads.fq

Sample input file:

R1:
@A00589:65:HCKLHDRXX:2:2101:1316:1016 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACCGTCCTTCCGGAAGTGAAACGATGTCAGCACTTGAATTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAAAGCAAGACGGCGAGAATCAGCACACCGACTTCTCGTGGACAGGAGCTACCACTTACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FF::F:FFFFFFFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFFFF,F:F:FF:,
@A00589:65:HCKLHDRXX:2:2101:3974:1016 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACTGTTGTTCCGGATCGGCTACGCCTTCACCCCTTGTATTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAACAAAAGACGGCGAGCAGAACAATACCGACTTCTCGTGGACAGGAGCTACCACATACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFF,FFF:FF:
@A00589:65:HCKLHDRXX:2:2101:1542:1031 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACAGCGGTGGTGGACACCCTACGATTTCAGGGCTTCAGTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAACACAAGACGGCGAGAATCACAACACCGACTTCTCGTGGACTGGAGCTACCACCTACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFFFFF::FFFFFFFFF:,FFFFF::FFFF,FFFFFFFFFFFFFFF:::FFFFFFFFFF::FFFFFFF,FFFFFFFFF:F:FF:FFFFFFFFFFFFFFFF,F::FF:,FFFFFFFFF:FFF:FFFFF,FF,F:FF
@A00589:65:HCKLHDRXX:2:2101:2465:1031 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACTGCTCTAGTGGAACAATTACGCTCTCACAGCTTTGTTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAACCAAAGACGGCGAGAATCACCCGACCGACTTCTCGTGGACGGGAGCTACCACCTACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFF,FFFFFFFF:FFF:F,FFFF:,FFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFF,FFFFFFF:FFFFF:FF:FFFFFFFFFFFF:FFFFFFFF:FFF,FFFFFF,FFFFFFFFFFFFF,:FFFF,
@A00589:65:HCKLHDRXX:2:2101:3477:1031 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACTCCGGATCGGGAAATCCGACGATATCACAACTTTGGTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAAAAACCGACGGCGAGACTCAGCACACCGACTTCTCGTGGACTGGAGCTACCACGTACCACCTCACTGGCAGAGACTCTCTGGTGAATCCGGGA
+
F,F:FF,:FFFFFFF:FFF:FFF,F:FFFFFFFFFFFF:FFFF::FFFF:FFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:,FFFFF:FFFFFFFFFF,FFFFFFFFF,FFFFFFFFFF:FFF,FFFFFFFFFFFFFF,FF:FFFFF:FFFF,FFFFF:FFF:FFFFFF:FFF,FFFFFFFFFFFF,FFFF,F:FFFF:FF,FF:F::FFFFFF:
@A00589:65:HCKLHDRXX:2:2101:4580:1031 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACAGTACAAGTGGAAGCATTACGCATTCACCCCTTCTGTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAAAACAAGACGGCGAGAATCAACAAACCGACTTCTCGAGGACAGGAGCTACCAAGTACCACCTCAATGGCAGAGACTCTCTGGTGAATCAGGGA
+
:FFFFFFFFFFFF:,:FF:FFF:FF:,FFFFF:FFFFFFF:FFFF,FFFFFFFF,:F:FFF,:FFFFFFFFF:FFFFF,FFFF:FFFFFFFFFFFFFFF,:FFFF:F:FFFFFF,FFF,F,,FF:FF:F::FFFFFFFFF,FFF:::FFFFFFFFF,FFFF:FF,FF:FFFFFFFFFFFFFFF,,F,F:F:FFF,F:::FFFFF,FFF:FF::F,FFF:,FFF,FFFFFFFFF,FFFFFFFFFFFF:F,FF


R2:
@A00589:65:HCKLHDRXX:2:2101:1316:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGGCCGTTGCCGCTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACGCTCTCTTCTTCAGCGTCTTTCGCCGGAGTATTCTCCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAATTGTTCTTCATCGTCACGGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFF,:FFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF:FFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFF:FF,F,FFFFFFFFFFFF:,:F,FFFFFFFFF,FFFFFF:F,F:FFFF:FF:::FFFFFFFFFFF::FF:F:FF:FFFFFFFFFFFFFF,FFFF:FF,FFF,F::FFFFFFF,F:FFF:FF:F,FFFF::FFFFFFFF:F
@A00589:65:HCKLHDRXX:2:2101:3974:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCGGCGGTAGCTGCTTGGGCGTTGCCTTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACCTTTTCAATGTCCACATTTGTTTTCTCTGAGCCTTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAACCGATCTTCATCGTCTTTGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFF
@A00589:65:HCKLHDRXX:2:2101:6415:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGACGGTTGCCTTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACTTGGTCGATCGCGACATTCTTCCCATGCGCACCCTCCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAAGCTGTCGTCATCGTCTCTGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF,FFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFF::FFF:FFFFFFFFFFFFFFF,FFFF::FFFFFF:
@A00589:65:HCKLHDRXX:2:2101:6632:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGCCGGTTGCCCCCCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACTCGATCCATATCCACATTGTTTTGTCCTGTCCATTCCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAATAACTTTCCATCATCGTCCTTGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
:FFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF,FFFFFFFFFFFFFFF,,FFFFFFF:FFF,FFFF:FFFFFFFFFF,:FFFFFFFFFF,FF:F,FFFFFFFFFFFFF
@A00589:65:HCKLHDRXX:2:2101:8404:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGTTCGTTGCCTTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCTTCTGTAATCATTACCCGTTCATCTGCTACGTCCTCCCCGTCTCTGCTTTCCTTCCCAAAGATGAGAAACCCGCTCTGAGGAAAAAAACGGTCCTCATCGTCAGTGTGGCGTGCCCTGGCCGGTCCCGGAGTCACCCGAGAGTC
+
FFFFF,FFF:FFF,FFFFFF::,,:FFFFFFFF,,,,FF,F:FF::F:FF:,FFF,FFFFFFF:FFF:,FFF::F,FFFF:,F:F,::F:,FF,:FFFFF,,,FFF:,FFFFF,FFFF:,FF:F,:FFF::::FF,:,F,:,FFFFFF:FFFFFFF,F:FF,:FF,F,FF:FFF,,FFF,FFFF,,FFF,FFFFF:,,,FF::FF::,FF,,,FF,F,FF:F,,,,FFF,FFF:,:F,:F,F,,FFFFF:F
@A00589:65:HCKLHDRXX:2:2101:12038:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGGCTGTTGCCTGCCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTGATCATGACGCCTTCGATCTCAGTGTTTTTGCCTGTGGTGCACTCCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAATTGATCTGCATCGTCACTGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFF,:FFF:FFFFFFFFFFFFF



Sample output file:

R1 out:
@A00589:65:HCKLHDRXX:2:2101:1316:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCTGCGGTAGCTGCTTGGCCGTTGCCGCTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACGCTCTCTTCTTCAGCGTCTTTCGCCGGAGTATTCTCCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAATTGTTCTTCATCGTCACGGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFF,:FFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF:FFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFF:FF,F,FFFFFFFFFFFF:,:F,FFFFFFFFF,FFFFFF:F,F:FFFF:FF:::FFFFFFFFFFF::FF:F:FF:FFFFFFFFFFFFFF,FFFF:FF,FFF,F::FFFFFFF,F:FFF:FF:F,FFFF::FFFFFFFF:F
@A00589:65:HCKLHDRXX:2:2101:3974:1016 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACTGTTGTTCCGGATCGGCTACGCCTTCACCCCTTGTATTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAACAAAAGACGGCGAGCAGAACAATACCGACTTCTCGTGGACAGGAGCTACCACATACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFF,FFF:FF:

R2 out:
@A00589:65:HCKLHDRXX:2:2101:1316:1016 1:N:0:TCCGCGAA+GTACTTAC
CGCTCAGTTCGTACCTGTATTTCTTGAGCAGAACAAACCGTCCTTCCGGAAGTGAAACGATGTCAGCACTTGAATTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAAAGCAAGACGGCGAGAATCAGCACACCGACTTCTCGTGGACAGGAGCTACCACTTACCACCTCAATGGCAGAGACTCTCTGGTGAATCCGGGA
+
FFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FF::F:FFFFFFFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFFFF,F:F:FF:,
@A00589:65:HCKLHDRXX:2:2101:3974:1016 2:N:0:TCCGCGAA+GTACTTAC
TCGTGGAGCGCGGCGGTAGCTGCTTGGGCGTTGCCTTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACCTTTTCAATGTCCACATTTGTTTTCTCTGAGCCTTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAACCGATCTTCATCGTCTTTGTGGCTTGCCATGGCCGGTCCCGGATTCACCAGAGAGTC
+
FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFF


"""
