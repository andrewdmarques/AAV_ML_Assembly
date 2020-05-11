#!/usr/bin/env python3

# Support Vector Machine 
# Source: https://labs.cognitiveclass.ai/tools/jupyterlab/lab/tree/labs/coursera/ML0101EN/ML0101EN-Clas-SVM-cancer-py-v1.ipynb
# Link for data file: https://s3-api.us-geo.objectstorage.softlayer.net/cf-courses-data/CognitiveClass/ML0101ENv3/labs/cell_samples.csv

import pandas as pd
import pylab as pl
import numpy as np
from numpy import asarray
from numpy import savetxt
import scipy.optimize as opt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import sys
import datetime
import pickle
from sklearn.metrics import classification_report, confusion_matrix
import itertools

start = datetime.datetime.now()
timer1 = str(start)
sys.stderr.write("Start time:			%s\n\nPython computations running...\n" % timer1)

#Open the sequence data.
test_df = pd.read_csv(sys.argv[2],sep='\t')

#Create the NumPy array, excluding the results column.
test_x = np.asarray(test_df.drop('Assembled', axis = 1))

#Create the NumPy array for the results column.
test_y = np.asarray(test_df['Assembled'])

#Retrieve the parameters data.
para_f_name = sys.argv[1]
with open(para_f_name, 'rb') as f:
	clf = pickle.load(f)

#Predict using the fitted model. 
yhat = clf.predict(test_x)
savetxt('rerun_' + sys.argv[1][:-18] + '_01_' + sys.argv[2][:-4] + '_predictions.csv', yhat, delimiter=',') 
#print('Sample of positive predictions:\t' , end = '')
#print(yhat[:5])
#print('Sample of negative predictions:\t' , end = '')
#print(yhat[-5:])

#Evaluation of SVM.
pred_a = 0
pred_u = 0
tot    = 0
for pred in yhat:
	tot += 1
	if pred == 1:
		pred_a += 1
	elif pred == 0:
		pred_u += 1
perc_pred_a = str(round(pred_a/tot * 100,2)) + '%'
#perc_pred_u = str(round(pred_u/tot * 100,2)) + '%'
print('\nPercent Assembled:\t' + perc_pred_a)
#stat = classification_report(test_y, yhat)
#print(stat)

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0) 


