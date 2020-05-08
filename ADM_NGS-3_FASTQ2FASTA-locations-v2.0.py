import sys
import datetime
import re

nt1 = "TTGAGCAG"

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\nFastq2Fasta computations running...\n" % timer0)

inf = sys.argv[1]
outf = open(inf + "_locations.fasta", "w+")

location = ""
i = 0

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases


with open(inf) as f:
	for line in f:
		if line[0] == "@":
			location = line[28:len(line) - 25]
			i = 1
		elif i == 1:
			outf.write( ">" + location + "\n")

			s1 = re.search(r'(?<=' + nt1 + ')\w+', line)
			if s1:				
				outf.write(line + "")
			else:
				line_rev = reverse_complement(line)
				outf.write(line_rev.strip() + "\n")
			i = 0
			location = ""


outf.close()

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Fastq2Fasta Completed:		%s\n\n" % timer_0)

"""
PURPOSE:
The purpose of this file is to take NGS FastQ files and generate FASTA files where the title of each sequence is the location of the sequence on the lane.
UPDATED: Now this script orients the reads so that they are all oriented the same way.

PROCEDURE:
Run this script with the next argument being the file to be read. 







Sample input file:

@A00589:65:HCKLHDRXX:2:2101:1045:1016 1:N:0:TCCGCGAA+GTACTTAC
TATAGTAGCTATCTGCGGTAGCTGCTTGTTCGTTGCCATGCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACATCTTCTATTGCTATATCTTTGTTTCGAGACCCTTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAACTTTTCTTCATCGTCCTTGTGGCTTGCCATGGCCGGTCCCGGATTCACCCGAGAG
+
:FFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:F:FF:FFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFF:FFFFFFFFFFFF,FFFFFFFFF,FFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:F:FFFFF,FFFFFFFF:FFFFFFFFF,,:F:FF
@A00589:65:HCKLHDRXX:2:2101:1081:1016 1:N:0:TCCGCGAA+GTACTTAC
TATAGTAGCTATCTGCGGTAGCTGCTTGCGTGTTGCCGTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACGTCGTCCAGCTCGACATTACAACTCGGCCCAGACTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAAGAGCCCGTCATCGTCACGGTGGCTTGCCATGGCCGGTCCCGGATTCACCCGAGAG
+
FFF:FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFF:,FFFFFFFFF:FFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,,FFFFFF,F,FFFF,F,FF,F,F,FF,FF,,,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F,FF,F:FFFFFF,F,F:FFF:FFFFF:FFFFFFFFFFF::F:F:FF,F:FF:
@A00589:65:HCKLHDRXX:2:2101:1099:1016 1:N:0:TCCGCGAA+GTACTTAC
TCATAGATTGTGTACTCTCATCGACCAGTACCTGTATTTCTTGAGCAGAACAAACGGCCAGAGCGGACTCAACACGCTCTCAGACCTTCGGTTTTCTCAGGCCGGAGCGAGTGACATTCGGGACCAGTCTAGGAACTGGCTTCCTGGACCCTGTTACCGCCAGCAGCGAGTATCAACAACACCCGCCGACAACAACAACAGTGACTTCTCGTGGCCTGGAGCTACCAAGTACCACCTCAATGGCAGAGACT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFFF:FFF:FF,FFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFF:F,FFFF::FFFF:FFFFFFFFFFF,F:F,FFFFFFF:FFFFFFFFFFFFFFF:F:F:FF,FF:FFF,:F,F:F,:FF:FF


Sample output file:

>1045:1016
TATAGTAGCTATCTGCGGTAGCTGCTTGTTCGTTGCCATGCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACATCTTCTATTGCTATATCTTTGTTTCGAGACCCTTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAACTTTTCTTCATCGTCCTTGTGGCTTGCCATGGCCGGTCCCGGATTCACCCGAGAG

>1081:1016
TATAGTAGCTATCTGCGGTAGCTGCTTGCGTGTTGCCGTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACGTCGTCCAGCTCGACATTACAACTCGGCCCAGACTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAAGAGCCCGTCATCGTCACGGTGGCTTGCCATGGCCGGTCCCGGATTCACCCGAGAG

>1082:1016
TATAGTAGCTATCTGCGGTAGCTGCTTGCGTGTTGCCGTTCTGGAGGTTGGTAGATACAGAACCATACTGCTCCGTGGCCACGGGATTGGTTGTCCTGATTTCCTCTTCGTCTGTAATCATGACGTCGTCCAGCTCGACATTACAACTCGGCCCAGACTGCTTCCCAAAGATGAGAACCCCGCTCTGAGGAAAAAAGAGCCCGTCATCGTCACGGTGGCTTGCCATGGCCGGTCCCGGATTCACCCGAGAG

"""