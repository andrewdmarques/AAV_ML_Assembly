#!/usr/bin/env python3

import numpy as np
import sys
import datetime
import re
import os

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\n" % timer0)

#Count the number of input sequences:
pre_ml = sys.argv[1]
num_seq = 0
with open(pre_ml,'r') as f:
	for seq in f:
		num_seq += 1


#Determine the composition of the dataset.
sys.stderr.write('\nIntroduction:\n' + str(num_seq) + "input sequences have been detected. At the prompt below, specify the number of assigned\npositive and negative datapoints in the set. The input sequences must be arranged with\npositive sequences at the start of the file and negative sequences at the end of the file.\nPress 'ENTER' at the prompt to auto-generate equal positive and negative results.\n\n")

num_pos = ''
num_neg = ''

num_pos = input('Number of positive datapoints: ')
if num_pos == '':
	print('\tAuto-generating results.')
	num_pos = round(num_seq/2)
	num_neg = num_seq - round(num_seq/2)
	print('\tPositive results: ' + str(num_pos))
	print('\tNegative results: ' + str(num_neg))
	
else:
	num_neg = input('\nNumber of negative datapoints: ')

#*****************************Pre-ml2Matrix****************************
sys.stderr.write("\nMatrix computations running...\n")

def matrix(residue):
	if residue == 'A':
		b = '0.17\t0.70\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'C':
		b = '0.29\t0.78\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0'
	elif residue == 'D':
		b = '0.30\t0.11\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0'
	elif residue == 'E':
		b = '0.47\t0.11\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0'
	elif residue == 'F':
		b = '0.77\t0.81\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0'
	elif residue == 'G':
		b = '0.00\t0.46\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'H':
		b = '0.56\t0.14\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0'
	elif residue == 'I':
		b = '0.64\t1.00\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'K':
		b = '0.65\t0.07\t1\t1\t0\t1\t0\t0\t0\t0\t0\t0'
	elif residue == 'L':
		b = '0.64\t0.92\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'M':
		b = '0.61\t0.71\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0'
	elif residue == 'N':
		b = '0.32\t0.11\t1\t1\t1\t0\t0\t0\t0\t0\t0\t1'
	elif residue == 'P':
		b = '0.31\t0.32\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'Q':
		b = '0.50\t0.11\t1\t1\t1\t0\t0\t0\t0\t0\t0\t1'
	elif residue == 'R':
		b = '0.68\t0.00\t1\t1\t0\t1\t0\t0\t0\t0\t0\t0'
	elif residue == 'S':
		b = '0.17\t0.41\t1\t1\t1\t0\t0\t0\t0\t0\t1\t0'
	elif residue == 'T':
		b = '0.33\t0.42\t1\t1\t1\t0\t0\t0\t0\t0\t1\t0'
	elif residue == 'V':
		b = '0.48\t0.97\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0'
	elif residue == 'W':
		b = '1.00\t0.40\t0\t1\t0\t0\t0\t0\t1\t0\t0\t0'
	elif residue == 'Y':
		b = '0.80\t0.36\t1\t1\t1\t0\t0\t0\t1\t0\t0\t0'
	return b

matrix_file = pre_ml[:-4] + '_PROPERTIES.csv'

matrix_out = open(matrix_file, 'w+')

#Write the headers into the CSV file
matrix_out.write('Volume\tHydropathy\tPolar\tH-Don\tH-Acc\tChargePos\tChargeNeg\tAliphatic\tAromatic\tSulfur\tHydroxyl\tAmide\t'*33 + 'Assembled\n' )

#Write the sequence date into the CSV file
seq_count = 0
with open(pre_ml,'r') as f:
	for seq in f:
		seq_count += 1
		i = 0
		len_seq = len(seq)
		for res in seq:
			if res == '\n':
				break
			res = res.strip()
			if i == 1:
				matrix_out.write('\t')
			i = 1
			bi = matrix(res)
			matrix_out.write(bi)
		if seq_count <= int(num_pos):
			matrix_out.write('\t1')
		elif seq_count > int(num_pos):
			matrix_out.write('\t0')
		matrix_out.write('\n')
matrix_out.close()

#Record the end time.
end_bi = datetime.datetime.now()
timer_bi = str(end_bi)
sys.stderr.write("Matrix Completed:		%s\n\n" % timer_bi)
'''
#*****************************Matrix2Transpose****************************

sys.stderr.write("Transposing computations running...\n")

#Record the end time.
end_tra = datetime.datetime.now()
timer_tra = str(end_tra)
sys.stderr.write("Transposing Completed:		%s\n\n" % timer_tra)
print('Saving file...')

#Name the numpy file.
npy_file = matrix_file[5:-11] + '_x'

#Create npy file.
x = np.loadtxt(matrix_file)
x = x.transpose()
np.save(npy_file, x)
os.remove(matrix_file)

#Record the end time.
end_tra = datetime.datetime.now()
timer_tra = str(end_tra)
sys.stderr.write("File saved:			%s\n\n" % timer_tra)

#*****************************X_npy2y_npy****************************

if num_pos != '':
	num_pos = int(num_pos)
	num_neg = int(num_neg)

# Define functions.
def load_x_data(infile):
	x = np.load(infile)
	return x

def load_y_data(x, num_pos, num_neg):
	if num_pos == '':
		num_pos = int(x.shape[1]/2)
		num_neg = int(x.shape[1]/2)
	if num_pos + num_neg == x.shape[1] - 1:
		num_neg += 1
	if num_pos + num_neg != x.shape[1]:
		print('****ERROR! Discrepency between number of x samples and y results****')
		print('****' + str(x.shape[1]) + ' input samples were detected and ' + str(num_pos + num_neg) + ' results were selected****')
	iterations = num_pos + num_neg
	outfile = 'temp_y-' + str(iterations) + '.csv'
	outf1 = open(outfile, 'w+')
	for x in range(0,iterations):
		if x == 0:
			outf1.write('1')
		elif x < num_pos:
			outf1.write('	1')
		else:
			outf1.write('	0')
	outf1.close()
	y_0dim = np.loadtxt(outfile)
	y = np.expand_dims(y_0dim, axis=0)
	os.remove(outfile)
	return y

# Run script.
if num_pos == '':
	print('Auto-generating equal positive and negative results')
elif type(num_pos) == int:
	print('Generating %s positive results and %s Negative results...' % (str(num_pos), str(num_neg)))
else: 
	print('****Error: the inputs for positive and negative results are not integers')


npy = npy_file + '.npy'
x = load_x_data(npy)
y = load_y_data(x, num_pos, num_neg)
np.save(npy[:-6] + '_y', y)
'''
#Record the end time.
end_y = datetime.datetime.now()
timer_y = str(end_y)
sys.stderr.write("Assigned results generated:	%s\n\n" % timer_y)
		
end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)

