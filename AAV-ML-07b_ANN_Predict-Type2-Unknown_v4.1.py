#!/usr/bin/env python3
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import pickle

# **********************************Defining basic Function**********************************
def predict(X, parameters):
    pos = 0
    neg = 0
    m = X.shape[1]
    n = len(parameters) // 2
    p = np.zeros((1,m))
    probas, caches = L_model_forward(X, parameters)
    for i in range(0, probas.shape[1]):
        if probas[0,i] > 0.5:
            p[0,i] = 1
            pos += 1
        else:
            p[0,i] = 0
            neg += 1
    per_pos = 100*(pos/(pos+neg))
    print('\nPercent positive results:	' + str(round(per_pos,2)) + '%\n')
    return p, per_pos
def L_model_forward(X, parameters):
    caches = []
    A = X
    L = len(parameters) // 2
    for l in range(1, L):
        A_prev = A 
        A, cache = linear_activation_forward(A_prev, parameters['W' + str(l)], parameters['b' + str(l)], activation = "relu")
        caches.append(cache)
    AL, cache = linear_activation_forward(A, parameters['W' + str(L)], parameters['b' + str(L)], activation = "sigmoid")
    caches.append(cache)
    assert(AL.shape == (1,X.shape[1]))
    return AL, caches
def linear_activation_forward(A_prev, W, b, activation):
    if activation == "sigmoid":
        Z, linear_cache = linear_forward(A_prev, W, b)
        A, activation_cache = sigmoid(Z)    
    elif activation == "relu":
        Z, linear_cache = linear_forward(A_prev, W, b)
        A, activation_cache = relu(Z)    
    assert (A.shape == (W.shape[0], A_prev.shape[1]))
    cache = (linear_cache, activation_cache)
    return A, cache
def linear_forward(A, W, b):
    Z = W.dot(A) + b
    assert(Z.shape == (W.shape[0], A.shape[1]))
    cache = (A, W, b)
    return Z, cache
def relu(Z):
    A = np.maximum(0,Z)    
    assert(A.shape == Z.shape)    
    cache = Z 
    return A, cache
def sigmoid(Z):    
    A = 1/(1+np.exp(-Z))
    cache = Z    
    return A, cache
# **********************************Script**********************************
# Load the parameters.
infile = sys.argv[1]
with open(infile, "rb") as pFile:
    parameters = pickle.load(pFile)

# Load the x values.
x = np.load(sys.argv[2])

# Predict the outcome.
p, per_pos = predict(x, parameters)


"""
PURPOSE: The purpose of this script is to predict the outcome of a known dataset using known parameters. 

PROCEDURE:
    Run: './DNN_predict_v3.0.py data06_200kst_5000it_0.1lr_660-440-1_parameters.txt data13_20copies-assembled_' 
"""
