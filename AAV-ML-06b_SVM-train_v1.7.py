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
train_df = pd.read_csv(sys.argv[1],sep='\t')
test_df = pd.read_csv(sys.argv[2],sep='\t')

#Create the NumPy array, excluding the results column.
train_x = np.asarray(train_df.drop('Assembled', axis = 1))
test_x = np.asarray(test_df.drop('Assembled', axis = 1))

#Create the NumPy array for the results column.
train_y = np.asarray(train_df['Assembled'])
test_y = np.asarray(test_df['Assembled'])

#SVM modeling.
from sklearn import svm
clf = svm.SVC(kernel='rbf')
clf.fit(train_x, train_y)

#Save Parameters for Trained Model.
para_f_name = 'svm_' + sys.argv[1][:-4] + '_01-parameters.txt'
with open(para_f_name, 'wb') as f:
	pickle.dump(clf, f)

#Predict using the fitted model. 
yhat = clf.predict(test_x)
savetxt('svm_' + sys.argv[1][:-4] + '_02_' + sys.argv[2][:-4] + '_predictions.csv', yhat, delimiter=',') 
print('Sample of positive predictions:\t' , end = '')
print(yhat[:5])
print('Sample of negative predictions:\t' , end = '')
print(yhat[-5:])

#Evaluation of SVM.
stat = classification_report(test_y, yhat)
print(stat)

#Alternative Statistics: determining the F1 score.
from sklearn.metrics import f1_score
f1 = str(f1_score(test_y, yhat, average='weighted'))
print('F1 score: ' + f1 + '\n')

end_0 = datetime.datetime.now()
timer_0 = str(end_0)
sys.stderr.write("Python Completed:		%s\n\n" % timer_0) 

#Save Stats to log file.
log_f_name = 'svm_' + sys.argv[1][:-4] + '_03_log.txt'
with open(log_f_name, 'w+') as log:
	log.write("Start time:			%s\n\nPython computations running...\n\n" % timer1)
	log.write('Train file:\t' + sys.argv[1] + '\nTest file:\t' + sys.argv[2] + '\n\nClassification Report:\n\n')
	log.write(stat)
	log.write('\nF1 score: ' + f1 + '\n\n')
	log.write("Python Completed:		%s" % timer_0)

