#!/usr/bin/env python3

import numpy as np
import sys
import datetime
import re
import os

start = datetime.datetime.now()
timer0 = str(start)
sys.stderr.write("Start time:			%s\n\n" % timer0)
#Determine the composition of the dataset.
sys.stderr.write("\nIntroduction:\nPositive datapoints will be assigned to the sequences in the beginning of the dataset. If\nyour data is not organized with positive results first, please format your data in this way.\nBelow, you will be prompted to type the number of positive and negative datapoints. The\ndefault will be set to equal number of positive and negative datasets. This will be applied\nif you hit 'ENTER' at the prompt.\n\n")

num_pos = ''
num_neg = ''

num_pos = input('Number of positive datapoints: ')
if num_pos == '':
	print('Default setting selected')
else:
	num_neg = input('\nNumber of negative datapoints: ')

#*****************************Pre-ml2Binary****************************
sys.stderr.write("\nBinary computations running...\n")

pre_ml = sys.argv[1]

def binary(residue):
  R = "ACDEFGHIKLMNPQRSTVWY"
  r = R.index(residue)
  return "0\t"*r+"1"+"\t0"*(len(R)-r-1)

binary_file = 'temp_' + pre_ml[:-4] + '_binary.csv'

binary_out = open(binary_file, 'w+')

with open(pre_ml,'r') as f:
	for seq in f:
		i = 0
		len_seq = len(seq)
		for res in seq:
			if res == '\n':
				break
			res = res.strip()
			if i == 1:
				binary_out.write(',')
			i = 1
			bi = binary(res)
			binary_out.write(bi)
		binary_out.write('\n')
binary_out.close()

#Record the end time.
end_bi = datetime.datetime.now()
timer_bi = str(end_bi)
sys.stderr.write("Binary Completed:		%s\n\n" % timer_bi)

#*****************************Binary2LinearTransformed****************************

sys.stderr.write("Transforming computations running...\n")

# Defining Functions.
def CSV2DNN(infile):
	file1 = infile
	file2 = infile[:-4] + '.DNN'
	outf1 = open(file2, 'w+')
	num_of_items = 660
	# Data transformation - switch rows and columns
	def retrieve_position(file,i):
		it = 0
		with open(file,'r') as f:
			out = ''
			for line in f:
				if it == 0:
					out += str(line[i])
					it = 1
				else:
					out += 	'	' + str(line[i])
			return out
	position = ''
	for x in range(0,num_of_items):
		#print(x + 1)
		outf1.write('')
		position += ''
		position += retrieve_position(file1,2*x)
		position += '\n'
	outf1.write(position)
	#Record the end time.
	end_tra = datetime.datetime.now()
	timer_tra = str(end_tra)
	sys.stderr.write("Transforming Completed:		%s\n\n" % timer_tra)
	print('Saving file...')
	outf1.close()

# Script to run the functions.
DNN_file = binary_file[:-4] + '.DNN'
npy_file = DNN_file[5:-11] + '_x'

# Create DNN file.
CSV2DNN(binary_file)

# Create npy file.
x = np.loadtxt(DNN_file)
np.save(npy_file, x)
os.remove(DNN_file)

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

#Record the end time.
end_y = datetime.datetime.now()
timer_y = str(end_y)
sys.stderr.write("Assigned results generated:	%s\n\n" % timer_y)
		
end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0)

