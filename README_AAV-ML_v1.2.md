# AAV ML

These series of scripts are intended for classifying outcomes of large protein libraries. Specifically, the original application of these scripts is to begin with paired end next generation sequencing data and convert the sequences into a format that the machine learning (ML) algorithms can read. Once converted, one of two algorithms (a support vector machine (SVM) or artificial neural network (ANN)) may be used to learn the datasets of two classifications. In the original application, sequences were classified as assembled or not assembled virus-like particles.

## Getting Started

To begin, move your sequence data, in the form of DNA sequences, to the desired directory containing the scripts. Once this is done, open a command line terminal. 

### Prerequisites

To use the scripts associated with this project, the following criteria must be met:
```
*A Linux-based computer.
*Python 3.7 or newer. 
*Datetime Python package.
*Itertools Python package.
*Matplotlib Python package.
*Numpy Python package.
*Os Python package.
*Pandas Python package.
*Pickle Python package.
*Pylab Python package.
*Random Python package.
*Re Python package.
*Scipy Python package.
*Sklearn Python package.
*String Python package.
*Python package.Python package.
*Python package.
*Python package.
```

## Usage

###Matching paired ends

AAV-ML-02_MatchPairedEnds.py is a Python 3 script written to process the NGS-gathered sequence data of libraries. Specifically, it ensures that paired end library reads have a matching forward and reverse read. After this script is used, the sequences will then be ready to be merged using any sequence merging software. We recommend using Flash by Magoc and Salzberg. 

```
./AAV-ML-02_MatchPairedEnds_v2.0.py R1-ForwardReads.fq R2-ReverseReads.fastq
```

###Extracting capsid library regions

AAV-ML-04a_CapLib.py is a Python 3 script written to process NGS-gathered sequence data of libraries. Specifically, it takes long sequences of libraries and extracts the variable residues, removing the constant residues from the sequences. Additionally, it condenses reads to ensure that only one unique variant is displayed. Optionally, the copy number of each unique variant can be included. 
```
./AAV-ML-04a_CapLib.py ref-04a.txt merged-reads.fastq
```

###Cleaning Sequences

AAV-ML-04b_CleanSequences.py is a Python 3 script written to process NGS-gathered sequence data of libraries. Specifically, it ensures that reads from AAV-ML-04a_CapLib.py result only in sequences that were intended to be synthesized by the degenerate primers.

```
./AAV-ML-04b_CleanSequences.py ref_4b.txt caplib_pre-ml.txt
```

###Pre-selection sequences

AAV-ML-04c_GeneratePre-Selection-Sequences.py is a Python 3 script written to produce sequences that represent the pre-selection distributions of residues that are proportionally accurate for residue distributions by using the degenerate codons used for library construction. After this step, the user must curate files for the training and testing sets with the positive and negative samples of sequences. To curate your sequences, have the sequences ready from the positive sample (those assembled believed to assemble into virus-like particles) and the sequences from the negative sample (those sequences believed to not assemble into virus-like particles). The curated file can have any ration of positive to negative sequences, but sequences that are positive must be grouped together and appear before the negative sequences. For example, a curated file may contain 200,000 lines, with each line being a sequence of variable amino acid residues where the first 100,000 lines represent a sampling of positive sequences and the last 100,000 lines represent a sampling of negative sequences.

```
./AAV-ML-04c_GeneratePre-Selection-Sequences.py ref-04c.txt
```

###Hypothetical library sequences
AAV-ML-04d_GenerateHypotheticalCapSeq.py is a Python 3 script written to produce sequences that represent the hypothetical mutation patterns for capsid libraries. To alter the number of sequences generated, set num_seqs equal to the desired number of sequences. To alter which amino acid position remains wild-type, change wt_pos to equal the position desired to be left not mutated. For CapLib8, use the same reference text file to generate these hypothetical mutations.

```
AAV-ML-04d_GenerateHypotheticalCapSeq.py ref-04c.txt
```

###Preparation of sequences for artificial neural network algorithm

AAV-ML-05a_MLprep_ANN.py is a Python 3 script written to process NGS-gathered sequence data of libraries. Specifically, sequences from AAV-ML-04b_CleanSequences_v2.0.py are converted into a binary matrix format that may be read by the artificial neural network. Once the script is run, indicate the number of positive and negative samples in the curated file. To choose the default setting where an equal number of positive and negative sequences will be used, then leave the prompt black and just press “ENTER”. 


```
./AAV-ML-05a_MLprep_ANN.py curated-seqs.txt
```

###Preparation of sequences for support vector machine using binary residues representation

 AAV-ML-05b_MLprep_SVMRESIDUES.py is a Python 3 script written to process NGS-gathered sequence data of libraries. Specifically, sequences from AAV-ML-04b_CleanSequences_v2.0.py are converted into a binary matrix format that may be read by the support vector machine. Once the script is run, indicate the number of positive and negative samples in the curated file. To choose the default setting where an equal number of positive and negative sequences will be used, then leave the prompt black and just press “ENTER”. 


```
./AAV-ML-05b_MLprep_SVMRESIDUES.py curated-seqs.txt
```

###Preparation of sequences for support vector machine using residues’ properties representation

AAV-ML-05c_MLprep_SVMPROPERTIES.py is a Python 3 script written to process NGS-gathered sequence data of libraries. Specifically, sequences from AAV-ML-04b_CleanSequences_v2.0.py are converted into a properties-based matrix format that may be read by the support vector machine. Once the script is run, indicate the number of positive and negative samples in the curated file. To choose the default setting where an equal number of positive and negative sequences will be used, then leave the prompt black and just press “ENTER”. 


```
./AAV-ML-05c_MLprep_SVMPROPERTIES.py curated-seqs.txty
```

###Training the artificial neural network

AAV-ML-06a_ANN-train.py is a Python 3 script written to train an artificial neural network. four files must be used to train the artificial neural network. File 1 is the features for the set that will train the algorithm. File 2 contains the assigned results (positive or negative) for the outcome of the samples features of file 1. File 3 is the features for the set that will be used to test/validate the algorithm. File 4 contains the assigned results (positive or negative) for the outcome of the samples features of file 3. Files 1 and 2, and files 3 and 4 must contain the same stem of the file, but end is x.csv or y.csv (i.e. train_x.npy, train_y.npy, test_x.npy, test_y.npy). For the input, all files must be in the same directory and use only the stem of the file names. 

The architecture of the neural network may be changed by increasing or decreasing the number of hidden layers and nodes within each layer using the layers_dims list. To add additional layers, add more numbers to the list. The number for each item in the list determines the number of nodes in the list. The first item in the list is for the input layer. The last item in the list is for the output layer. 

To alter the number of learning iterations, set num_iteration equal to the desired number.

To alter the learning rate, set learning_rate equal to the desired value. 

```
./AAV-ML-06a_ANN-train.py train_ test_
```


###Training the support vector machine

AAV-ML-06b_SVM-train.py is a Python3 is a Python 3 script written to train a support vector machine.Two files must be used to train the support vector machine. File 1 contains the features and known classification outcome for the sample that will be used to train the algorithm. File 2 contains the features and outcomes that will be used to test/validate the algorithm. Both the training and the testing sets must be the same type of representation.

To alter the type of kernel used in the support vector machine, set the kernel variable equal to the designed kernel (either ‘linear’, ‘poly’, ‘rbf’, or ‘sigmoid’).

```
./AAV-ML-06b_SVM-train.py train_RESIDUES.csv test_RESIDUES.csv
```

###Predictions for known outcomes using a trained artificial neural network

AAV-ML-07a_ANN-Predict-Type1-Known.py is a Python 3 script written to use the trained artificial neural network parameters to predict samples of known outcomes. 3 files are required to run this script. File 1 contains the trained parameters of the neural network. File 2 contains the features of the sequences which will be predicted. File 3 contains the known classification outcomes for the features that will be predicted. Files 2 and 3 should contain the same stem of the file, but end is x.csv or y.csv. Only the stem of files 2 and 3 need to be submitted to the script, the script will then open both of these files automatically. For example seqs_x.npy and seqs_y.npy will be called by typing “seqs_” one time.

```
./AAV-ML-07a_ANN-Predict-Type1-Known.py parameters.txt seqs_
```

###Predictions for unknown outcomes using a trained artificial neural network

AAV-ML-07b_ANN_Predict-Type2-Unknown.py is a Python 3 script written to use the trained artificial neural network parameters to predict samples of unknown outcomes. 2 files are required to run this script. File 1 contains the trained parameters of the neural network. File 2 contains the features of the sequences which will be predicted.

```
./AAV-ML-07b_ANN_Predict-Type2-Unknown.py parameters.txt seqs_x.npy
```


###Predictions for known outcomes using a support vector machine

AAV-ML-07c_SVM_Predict-Type1-Known.py is a Python 3 script written to use the train support vector machine parameters to predict samples of known outcomes. 2 files are required to run this script. File 1 contains the trained parameters of the support vector machine. File 2 contains the features and known outcomes of the sequences which will be predicted.

```
AAV-ML-07c_SVM_Predict-Type1-Known.py parameters.txt seqs.npy
```

###Predictions for unknown outcomes using a support vector machine

AAV-ML-07d_SVM-Predict-Type2-Unknown.py is a Python 3 script written to use the train support vector machine parameters to predict samples of unknown outcomes. 2 files are required to run this script. File 1 contains the trained parameters of the support vector machine. File 2 contains the features of the sequences which will be predicted.

```
AAV-ML-07d_SVM-Predict-Type2-Unknown.py parameters.txt seqs.npy
```

## Authors

* **Andrew David Marques** - *Initial work* - [University of Florida](https://www.linkedin.com/in/andrew-marques-290a29164/)

* **Michael Kummer** - *Review* - [University of Florida](https://www.linkedin.com/in/michael-kummer-034b57194/)

## License

This project is licensed under the UF License - see the LICENSE.md file for details




