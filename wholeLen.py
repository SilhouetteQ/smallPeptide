import re
import pickle
import sys
from pickle import dump
import pandas as pd
import numpy as np
from download_fasta import seq2ngram

def save_obj(obj, name):
    with open('../model/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def load_obj(name):
    with open('../model/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
def threeNTCount(sam):
    # threeNTPath = '../../../data/' + sam + '/ribo/final_RiboCode_ORFs_result.txt'
    threeNTPath = '../../../data/' + sam + '/ribo/RiboCode_ORFs_result.txt'
    df = pd.read_table(threeNTPath, low_memory = False)
    df1 = df.loc[df['gene_type'] == 'lncRNA', ['transcript_id', 'ORF_tstart', 'ORF_tstop']]
    dID = []
    seqID = df1['transcript_id'].to_list()
    dID += [item[:15] for item in seqID]
    # print('unique 3-nt transcripts number for ' + sam + ': ', len(set(dID)))
    # dID = []
    # threeNTPath = '../../../data/' + sam + '/ribo/RiboCode_ORFs_result.txt'
    # df = pd.read_table(threeNTPath, low_memory = False)
    # df1 = df.loc[df['gene_type'] == 'lncRNA', ['transcript_id', 'ORF_tstart', 'ORF_tstop']]
    # seqID = df1['transcript_id'].to_list()
    # dID += [item[:15] for item in seqID]
    # print(len(set(dID)))
    return df1
def findORF(re2, orfCutoff, cutStart):
    substring = re.findall('(ATG)((\w\w\w)*?)((TAA)|(TGA)|(TAG))', re2)
    sindex = [m.start() + 1 for m in re.finditer('(ATG)((\w\w\w)*?)((TAA)|(TGA)|(TAG))', re2)]
    # print(len(substring[0][0] + substring[0][1] + substring[0][3]))
    # print(sindex)
    ssDict = {}
    # each key has corresponding values at most 3 and at least 1
    for num, sstr in enumerate(sindex):
        # store ORF start and corresponding sequence in a dictionary
        ssDict[sstr] = substring[num][0] + substring[num][1] + substring[num][3]
    for item in list(ssDict.keys()):
        if len(ssDict[item]) < orfCutoff: ### setting ORF length
            del ssDict[item]
            sindex.remove(item)
            continue
        if item < cutStart + 1: ### setting starting cutting sites
            del ssDict[item]
            sindex.remove(item)
    return sindex, ssDict

def seqUpDown(sam, type, cutLen, cutStart, orfCutoff, cutRank):
    if type == 'ribo':
        seqDict = load_obj('threeNTFasta')
        seqInfo = threeNTCount(sam)
        seqInfo['t_id'] = [item[:15] for item in seqInfo['transcript_id']]
        emptySeq = ['ENST00000675070', 'ENST00000408966', 'ENST00000536039', 'ENST00000235290']
        t_id = []
        t_id1 = []
        wholeLen = []
        for i, id in enumerate(seqInfo['t_id']):
            if id in emptySeq: # this is the ENST transcript that can not be downloaded
                continue
            seq = seqDict[id] # retrive sequence
            cut = seqInfo['ORF_tstart'].iloc[i] # retrive transcript start
            cut = int(cut)
            if re.search('ATG', seq) == None:  # remove transcripts that have no ATG
                continue
            if cut < cutStart + 1: # remove transcripts that having 5`UTR length shorter than cutLen
                continue
            re2 = ''.join(re.split('\n', seq)[1:])
            sindex, qualify = findORF(re2, orfCutoff, cutStart) # find qualified ORF
            if len(re2) > 15000:
                continue
            if len(qualify) < cutRank+1:
                continue
            t_id1.append(id)
            if sindex[cutRank] == cut:
                t_id.append(id)
                wholeLen.append(re2)
        print('All possible unique transcript ID: ', len(set(t_id1)))
        print('Unique transcript ID on site ', cutRank, ' : ', len(set(t_id)))
        save_obj(set(wholeLen), sam + 'wCUT' + type)
        return wholeLen, t_id
    if type == 'rna':
        seqDict = load_obj('rnaCountFasta')
        emptySeq = load_obj('missingrnaCountFasta')
        a = load_obj('rnaTranscript')
        samID = a[sam].to_list()
        odID = []
        odID = [item[:15] for item in samID] # transcript overlap with ribo 3nt should be removed
        threeNTdID = []
        # transcript overlap with ribo-seq should be removed
        for s in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
            seqInfo = threeNTCount(s)
            seqID = seqInfo['transcript_id'].to_list()
            threeNTdID += [item[:15] for item in seqID]
        dID = [value for value in odID if value not in threeNTdID]
        wholeLen = []
        t_id = []
        t_id1 = []
        for i, id in enumerate(dID):
            if id in emptySeq:
                continue
            seq = seqDict[id]
            re2 = ''.join(re.split('\n', seq)[1:])
            sindex, qualify = findORF(re2, orfCutoff, cutStart)
            if len(re2) > 15000:
                continue
            if len(qualify) < cutRank + 1:
                continue
            wholeLen.append(re2)
            t_id.append(id)
        save_obj(set(wholeLen), sam + 'wCUT' + type)
        return wholeLen, t_id
# seqUpDown('liver', 'ribo', 10, 10, 50, 0)
# seqUpDown('liver', 'rna', 10, 10, 50, 0)

pos = load_obj('liver' + 'wCUT' + 'ribo')
neg = load_obj('liver' + 'wCUT' + 'rna')
## overlap detect
overlap = [value for value in list(pos) if value in list(neg)]
overlap = set(overlap)
if len(overlap) > 0:
	print('overwhelming ' + str(len(overlap)) + ' overlaps')
	print(overlap)
	for item in overlap:
		# pos.remove(str(item))
		neg.remove(str(item))
print('need to n-gram %d seqs' % len(pos))
print('need to n-gram %d seqs' % len(neg))
def upseq2ngram(X, k):
    XK = []
    for line in X:
        kmer = ''
        line = line.lower()
        l = len(line)
        s = 1
        for i in range(0, l, s):
            if i + k >= l + 1:
                break
            kTmp = ''.join(line[i:i + k])
            kmer += kTmp
            kmer += ' '
        kmer += '\n'
        XK.append(kmer)
    return XK

posK = upseq2ngram(pos, 4)
negK = upseq2ngram(neg, 4)
# pos_seqs = posK.split('\n')
pos_seqs = [line[:-2] for line in posK]
neg_seqs = [line[:-1] for line in negK]

seqs = pos_seqs + neg_seqs
max_length = max([len(s.split()) for s in seqs])

from keras.preprocessing.text import Tokenizer
from sklearn.model_selection import train_test_split, GridSearchCV
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Embedding
from keras.preprocessing.sequence import pad_sequences
from keras import callbacks


tokenizer = Tokenizer()
tokenizer.fit_on_texts(seqs)
sequences = tokenizer.texts_to_sequences(seqs)
vocab_size = len(tokenizer.word_index) + 1
pad_seq = pad_sequences(sequences, maxlen=max_length, padding='post')

X = pad_seq
y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))

seq_length = X.shape[1]

X_train, X_test, y_train, y_test = train_test_split(
	X, y, test_size=0.33, random_state=42)

model = Sequential()
model.add(Embedding(vocab_size, 50, input_length=seq_length))
model.add(LSTM(100, return_sequences=True))
model.add(LSTM(100))
model.add(Dense(100, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
# model.add(Dense(vocab_size, activation='softmax'))
print(model.summary())
# compile model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['AUC'])
# fit model
model.fit(X_train, y_train, batch_size=128, epochs=60,
		  callbacks = callbacks.EarlyStopping(monitor='loss', patience=3))
y_pred = model.predict_classes(X_test, verbose=0)

from sklearn.metrics.cluster import contingency_matrix
print(contingency_matrix(y_test, y_pred))
print(roc_auc_score(y_test, y_pred))
print(f1_score(y_test, y_pred))


