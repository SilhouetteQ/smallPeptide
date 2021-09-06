import sys
import pickle
import numpy as np
from numpy import array
import pandas as pd
from pickle import dump
from tensorflow.keras.utils import to_categorical
from download_fasta import seq2ngram
from keras.preprocessing.text import Tokenizer
from sklearn.model_selection import train_test_split, GridSearchCV
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Embedding

# load doc into memory
def load_doc(filename):
	# open the file as read only
	file = open(filename, 'r')
	# read all text
	text = file.read()
	# close the file
	file.close()
	return text

# infile = './fasta/liverribo15_seqonly_Uniq.fa'
infile1 = './fasta/breastribo15_seqonly_Uniq.fa'
infile2 = './fasta/breastrna15_seqonly_Uniq.fa'
def seqProcess(infile1, infile2):
	raw_text1 = load_doc(infile1)
	raw_text2 = load_doc(infile2)
	raw_text =  raw_text1 + raw_text2
	t1 = len(raw_text1.split('\n')[:-1])
	t2 = len(raw_text2.split('\n')[:-1])
	print(t1, t2)
	y = np.array([1] * t1 + [0] * t2)
	lines = raw_text.split('\n')
	kmer = seq2ngram(lines, 4)
	kmer_seq = [line[:-1] for line in kmer.splitlines()]
	tokenizer = Tokenizer()
	tokenizer.fit_on_texts(kmer_seq)
	sequences = tokenizer.texts_to_sequences(kmer_seq)
	sequences = sequences[:-1]
	xt = sequences[0]
	for a in sequences[1:]:
		xt = np.vstack((xt, a))
	X = np.array(xt)
	print(X)
	sys.exit()
	vocab_size = len(tokenizer.word_index) + 1
	return vocab_size, X, y

vocab_size, X, y = seqProcess(infile1, infile2)

# def load_obj(name):
#     with open('../model/' + name + '.pkl', 'rb') as f:
#         return pickle.load(f)
#
# from keras.preprocessing.text import Tokenizer
# from keras.preprocessing.sequence import pad_sequences
#
# pos = load_obj('liver' + 'CUT' + 'ribo')
# neg = load_obj('liver' + 'CUT' + 'rna')
# overlap = [value for value in list(pos) if value in list(neg)]
# overlap = set(overlap)
# if len(overlap) > 0:
# 	print('overwhelming ' + str(len(overlap)) + ' overlaps')
# 	print(overlap)
# 	# print(pos)
# 	# print(neg)
# 	for item in overlap:
# 		# pos.remove(str(item))
# 		neg.remove(str(item))
# print('need to n-gram %d seqs' % len(pos))
# print('need to n-gram %d seqs' % len(neg))
# posK = seq2ngram(pos, 5)
# negK = seq2ngram(neg, 5)
# pos_seqs = [line[:-1] for line in posK.splitlines()]
# neg_seqs = [line[:-1] for line in negK.splitlines()]
# seqs = pos_seqs + neg_seqs
# tokenizer = Tokenizer()
# tokenizer.fit_on_texts(seqs)
# sequences = tokenizer.texts_to_sequences(seqs)
# sequences = np.array(sequences)
#
# vocab_size = len(tokenizer.word_index) + 1
#
# X = sequences
# y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
#
#
# # infile1 = './fasta/liverribo15_seqonly_Uniq.fa'
# # infile2 = './fasta/liverrna15_seqonly_Uniq.fa'
# # pos = seqProcess(infile1)
# # neg = seqProcess(infile2)
# # X = np.append(pos, neg)
# # y = np.array([1] * len(pos) + [0] * len(neg))
#

from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE

X_train, X_test, y_train, y_test = train_test_split(
	X, y, test_size=0.33, random_state=42)
print(X.shape)
print(X_train.shape)
print(X_test.shape)

over = SMOTE(sampling_strategy=0.5)
under = RandomUnderSampler(sampling_strategy=0.5)
# X_train, y_train = under.fit_resample(X_train, y_train)

print(X.shape)
print(X_train.shape)
print(X_test.shape)
seq_length = X_train.shape[1]
from keras import callbacks
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Flatten
from keras.layers import Embedding
from keras.layers.convolutional import Conv1D
from keras.layers.convolutional import MaxPooling1D
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, f1_score

model = Sequential()
model.add(Embedding(vocab_size, 50, input_length=seq_length))
model.add(Conv1D(filters=32, kernel_size=8, activation='relu'))
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(10, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
print(model.summary())
# # compile model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['AUC'])
# # fit model
model.fit(X_train, y_train, batch_size=128, epochs=60,
		  callbacks = callbacks.EarlyStopping(monitor='loss', patience=3))

# loss, accuracy = model.evaluate(padded_docs, labels, verbose=0)


# # save the model to file
# model.save('model.h5')
# # save the tokenizer
# dump(tokenizer, open('tokenizer.pkl', 'wb'))
q = model.evaluate(X_test, y_test)
print(q)
y_pred = model.predict_classes(X_test, verbose=0)
from sklearn.metrics.cluster import contingency_matrix
print(contingency_matrix(y_test, y_pred))
print(roc_auc_score(y_test, y_pred))
print(f1_score(y_test, y_pred))
