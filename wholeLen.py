### python 3.8
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
    threeNTPath = '../../../data/' + sam + '/ribo/final_RiboCode_ORFs_result.txt'
    # threeNTPath = '../../../data/' + sam + '/ribo/RiboCode_ORFs_result.txt'
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

def seqUpDown(sam, type, cutLen, cutStart, orfCutoff, cutRank, thresh = 1000):
    if type == 'ribo':
        seqDict = load_obj('threeNTFasta')
        seqInfo = threeNTCount(sam)
        seqInfo['t_id'] = [item[:15] for item in seqInfo['transcript_id']]
        emptySeq = ['ENST00000675070', 'ENST00000408966', 'ENST00000536039', 'ENST00000235290']
        t_id = []
        t_id1 = []
        pre = []
        post = []
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
            seqLen = len(re2)
            sindex, qualify = findORF(re2, orfCutoff, cutStart) # find qualified ORF

            if len(qualify) < cutRank+1:
                continue
            if sindex[cutRank] == cut:
                ## b & c
                t_id.append(id)
                ## b
                # wholeLen.append(re2)
                if cut < thresh:
                    pre.append(re2[:cut-1])
                    if seqLen - cut < thresh + 3:
                        post.append(re2[(cut+2):])
                    else:
                        post.append(re2[(cut+2):(cut+2+thresh)])
                else:
                    pre.append(re2[(cut-thresh-1):(cut-1)])
                    if seqLen - cut < thresh + 3:
                        post.append(re2[(cut + 2):])
                    else:
                        post.append(re2[(cut+2):(cut+2+thresh)])
            ## a
            # t_id.append(id)
            # wholeLen.append(re2)
        padPre = [line.rjust(thresh, '0') for line in pre]
        padPost = [line.ljust(thresh, '0') for line in post]
        ## c
        wholeLen = list(map(lambda x, y: x + 'ATG' + y, padPre, padPost))
        print(len(set(wholeLen)))
        print(len(set(t_id)))
        print(len(set(t_id1)))
        print('All possible unique transcript ID: ', len(set(t_id1)))
        print('Unique transcript ID on site ', cutRank, ' : ', len(set(t_id)))
        # save_obj(set(wholeLen), sam + 'wCUT' + type)
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
        pre = []
        post = []
        t_id = []
        wholeLen = []
        for i, id in enumerate(dID):
            if id in emptySeq:
                continue
            seq = seqDict[id]
            re2 = ''.join(re.split('\n', seq)[1:])
            seqLen = len(re2)
            sindex, qualify = findORF(re2, orfCutoff, cutStart)
            # if len(re2) > 15000:
            #     continue
            if len(qualify) < cutRank + 1:
                continue
            t_id.append(id)
            # a & b
            # wholeLen.append(re2)
            cut = sindex[cutRank]
            if cut < thresh:
                pre.append(re2[:cut-1])
                if seqLen - cut < thresh + 3:
                    post.append(re2[(cut+2):])
                else:
                    post.append(re2[(cut+2):(cut+2+thresh)])
            else:
                pre.append(re2[(cut-thresh-1):(cut-1)])
                if seqLen - cut < thresh + 3:
                    post.append(re2[(cut + 2):])
                else:
                    post.append(re2[(cut+2):(cut+2+thresh)])
        padPre = [line.rjust(thresh, '0') for line in pre]
        padPost = [line.ljust(thresh, '0') for line in post]
        ## c
        wholeLen = list(map(lambda x, y: x + 'ATG' + y, padPre, padPost))
        print(len(set(wholeLen)))
        print(len(set(t_id)))
        # save_obj(set(wholeLen), sam + 'wCUT' + type)
        return wholeLen, t_id
# seqUpDown(sam, 'ribo', 15, 15, 50, 0)
# seqUpDown(sam, 'rna', 15, 15, 50, 0)

def writeFasta(sam, type, cutLen, cutStart, orfCutoff, cutRank, enst=True):
    cutSeq, t_id = seqUpDown(sam, type, cutLen, cutStart, orfCutoff, cutRank)
    if enst is False:
        seqPath = './fasta/' + sam + type + str(cutLen) + '_seqonly.fa'
        f = open(seqPath, 'w')
        for item in set(cutSeq):
            f.write(item)
            f.write('\n')
        f.close()
    else:
        forMotif = './fasta/' + sam + type + str(cutLen)+ '_enst.fa'
        f1 = open(forMotif, 'w')
        rep = []
        print(cutSeq[12])
        for m, item in enumerate(cutSeq):
            if item in rep:
                continue
            rep.append(item)
            f1.write('>')
            f1.write(t_id[m])
            f1.write('\n')
            f1.write(item)
            f1.write('\n')
            # f1.write(item)
        f1.close()

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
            if '0' in kTmp:
                kTmp = '0' * k
            kmer += kTmp
            kmer += ' '
        kmer += '\n'
        XK.append(kmer)
    return XK

def kmer_solution(sam, k):
    pos = load_obj(sam + 'wCUT' + 'ribo')
    neg = load_obj(sam + 'wCUT' + 'rna')
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
    posK = upseq2ngram(pos,k)
    negK = upseq2ngram(neg,k)
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
    tokenizer.word_index['0'* k] = 0
    pad_seq = pad_sequences(sequences, maxlen=max_length, padding='post')
    seqlen = [len(line) for line in sequences]

    X = pad_seq
    y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))

    seq_length = X.shape[1]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42)
    model = Sequential()
    model.add(Embedding(vocab_size, 50, input_length=seq_length))
    # model.add(Conv1D(filters=32, kernel_size=3, padding='same', activation='relu'))
    # model.add(MaxPooling1D(pool_size=2))
    model.add(LSTM(100, return_sequences=True))
    model.add(LSTM(100))
    # model.add(LSTM(100, dropout=0.2, recurrent_dropout=0.2))
    model.add(Dense(100, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    # model.add(Dense(vocab_size, activation='softmax'))
    print(model.summary())
    # compile model
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    # fit model
    model.fit(X_train, y_train, batch_size=128, epochs=60,
              callbacks = callbacks.EarlyStopping(monitor='loss', patience=3))

    scores = model.evaluate(X_test, y_test, verbose=0)
    print("Accuracy: %.2f%%" % (scores[1] * 100))

    sys.exit()
    from keras.wrappers.scikit_learn import KerasClassifier
    from sklearn.model_selection import cross_val_score
    from sklearn.model_selection import StratifiedKFold

    # def create_baseline():
    #     # create model
    #     model = Sequential()
    #     model.add(Embedding(vocab_size, 50, input_length=seq_length))
    #     model.add(LSTM(100, return_sequences=True))
    #     model.add(LSTM(100))
    #     model.add(Dense(100, activation='relu'))
    #     model.add(Dense(1, activation='sigmoid'))
    #     # Compile model
    #     model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['AUC'])
    #     return model
    # estimator = KerasClassifier(build_fn=create_baseline, epochs=60, batch_size=5, verbose=0)
    # kfold = StratifiedKFold(n_splits=5, shuffle=True)
    # results = cross_val_score(estimator, X, encoded_Y, cv=kfold)
    # print("Baseline: %.2f%% (%.2f%%)" % (results.mean() * 100, results.std() * 100))

    from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, f1_score
    from sklearn.metrics.cluster import contingency_matrix
    y_pred = model.predict(X_test)
    seq_predictions = np.transpose(y_pred)[0]  # transformation to get (n,)
    # Applying transformation to get binary values predictions with 0.5 as thresold
    y_pred = list(map(lambda x: 0 if x < 0.8 else 1, seq_predictions))
    print(contingency_matrix(y_test, y_pred))
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
    y_pred = list(map(lambda x: 0 if x < 0.6 else 1, seq_predictions))
    print(contingency_matrix(y_test, y_pred))
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
    y_pred = list(map(lambda x: 0 if x < 0.4 else 1, seq_predictions))
    print(contingency_matrix(y_test, y_pred))
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
    y_pred = list(map(lambda x: 0 if x < 0.2 else 1, seq_predictions))
    print(contingency_matrix(y_test, y_pred))
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))

    y_pred = list(map(lambda x: 0 if x < 0.1 else 1, seq_predictions))
    print(contingency_matrix(y_test, y_pred))
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
# kmer_solution(sam,4)

def knn(X_train, y_train, X_test, y_test, cv):
    from sklearn.neighbors import KNeighborsClassifier  # The k-nearest neighbor classifier
    from sklearn.feature_selection import VarianceThreshold  # Feature selector
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import Normalizer, StandardScaler, MinMaxScaler, PowerTransformer, MaxAbsScaler
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('selector', VarianceThreshold()),
        ('classifier', KNeighborsClassifier())
    ])
    parameters = {'scaler': [StandardScaler(), MinMaxScaler(),
                             Normalizer(), MaxAbsScaler()],
                  'selector__threshold': [0, 0.001],
                  'classifier__n_neighbors': [1, 3, 5, 7, 10],
                  'classifier__p': [1, 2],
                  'classifier__leaf_size': [1, 5, 10, 15]
                  }
    grid = GridSearchCV(pipe, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid
def rf(X_train, y_train, X_test, y_test, cv):
    from imblearn.ensemble import BalancedRandomForestClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.feature_selection import VarianceThreshold  # Feature selector
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import Normalizer, StandardScaler, MinMaxScaler, PowerTransformer, MaxAbsScaler
    print(X_train.shape)
    print(X_test.shape)
    model = BalancedRandomForestClassifier()
    parameters = {
        'bootstrap': [True],
        'max_depth': [60, 70, 80],
        'max_features': ['auto'],
        'min_samples_leaf': [1, 3],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [30, 40, 50]
    }
    # parameters = {'max_depth': [3, 5, 10], 'min_samples_split': [2, 5, 10]}
    grid = GridSearchCV(model, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid
def bag(X_train, y_train, X_test, y_test, cv):
    from imblearn.ensemble import BalancedBaggingClassifier
    from sklearn.feature_selection import VarianceThreshold  # Feature selector
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import Normalizer, StandardScaler, MinMaxScaler, PowerTransformer, MaxAbsScaler
    print(X_train.shape)
    print(X_test.shape)
    model = BalancedBaggingClassifier()
    parameters = {
        'bootstrap': [True],
        'n_estimators': [30, 40, 50]
    }
    grid = GridSearchCV(model, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid
def rus(X_train, y_train, X_test, y_test, cv):
    from imblearn.ensemble import RUSBoostClassifier
    from sklearn.model_selection import GridSearchCV
    print(X_train.shape)
    print(X_test.shape)
    model = RUSBoostClassifier()
    parameters = {
        'n_estimators': [20, 50, 80]
    }
    grid = GridSearchCV(model, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid
def ens(X_train, y_train, X_test, y_test, cv):
    from imblearn.ensemble import EasyEnsembleClassifier
    from sklearn.model_selection import GridSearchCV
    print(X_train.shape)
    print(X_test.shape)
    model = EasyEnsembleClassifier()
    parameters = {
        'n_estimators': [70, 80, 90]
    }
    grid = GridSearchCV(model, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid
def svm(X_train, y_train, X_test, y_test, cv):
    from sklearn.feature_selection import VarianceThreshold  # Feature selector
    from sklearn.svm import SVC
    from sklearn.model_selection import GridSearchCV
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.preprocessing import Normalizer, StandardScaler, MinMaxScaler, PowerTransformer, MaxAbsScaler
    parameters = {'kernel': ['rbf'], 'gamma': [0.1, 0.01, 1e-3, 1e-4],
                  'C': [1, 10, 100, 1000], 'class_weight': ['balanced', None]}
    from sklearn.utils.fixes import loguniform
    parameters = {'C': loguniform(1e0, 1e3), 'gamma': loguniform(1e-4, 1e-3),
                  'kernel': ['rbf'], 'class_weight': ['balanced', None]}
    model = SVC()
    grid = RandomizedSearchCV(model, parameters, cv=cv).fit(X_train, y_train)
    print('Training set score: ' + str(grid.score(X_train, y_train)))
    print('Test set score: ' + str(grid.score(X_test, y_test)))
    return grid

def curated():
    import sys
    from numpy import mean, arange
    from sklearn.model_selection import train_test_split, GridSearchCV, RepeatedStratifiedKFold
    from sklearn.metrics.cluster import contingency_matrix
    from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, f1_score
    sys.path.append("../../../../lncRNAIdentification")
    from sequence_attributes_orf import SequenceAttributes
    # from sequence_attributes_kmer import SequenceAttributes
    pos = SequenceAttributes('./fasta/'+sam+'ribo15_enst.fa', 1)
    neg = SequenceAttributes('./fasta/'+sam+'rna15_enst.fa', 0)
    # pos = SequenceAttributes('./fasta/breastfibroblastHEK293_2liverribo30_enst_Uniq.fa', 1)
    # neg = SequenceAttributes('./fasta/breastfibroblastHEK293_2liverrna30_enst_Uniq.fa', 0)
    pos_df = pd.DataFrame(pos.process())
    neg_df = pd.DataFrame(neg.process())
    # df.to_csv('trainingfile.txt', index=False)
    # X = pd.concat([pos_df, neg_df.iloc[:pos_df.shape[0]]])
    X = pd.concat([pos_df, neg_df])
    y = X['class']
    X = X.drop(columns=['id', 'class'])
    # X = X.drop(columns=['id', 'class', 'length', 'fl', 'fp', 'll', 'lp',
    #                     'getoentropy3', 'getoentropy4', 'getoentropy5',
    #                     'toentropy3', 'toentropy4', 'toentropy5'])

    # X, y = prepareData(sam, 5, cutLen)
    # X[:, 0] = y
    # unique_rows, indices, occurrence_count = np.unique(
    #     X, axis = 0, return_counts = True, return_index = True)
    # X = unique_rows
    # y =  np.copy(X[:, 0])
    # X[:, 0] = 0
    # X = np.array(X)
    # y = np.array(y)
    print('dimension of positive data', pos_df.shape)
    print('dimension of negative data', neg_df.shape)
    print('dimension of all data', X.shape)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=10, stratify=y)

    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=5, random_state=1)
    grid = knn(X_train, y_train, X_test, y_test, cv)
    grid = rf(X_train, y_train, X_test, y_test, cv)
    grid = bag(X_train, y_train, X_test, y_test, cv)
    grid = rus(X_train, y_train, X_test, y_test, cv)
    grid = ens(X_train, y_train, X_test, y_test, cv)
    grid = svm(X_train, y_train, X_test, y_test, cv)

    # # Access the best set of parameters
    # best_params = grid.best_params_
    # print(best_params)
    # # Stores the optimum model in best_pipe
    # best_pipe = grid.best_estimator_
    # print(best_pipe)

    # result_df = pd.DataFrame.from_dict(grid.cv_results_, orient='columns')
    # print(result_df.columns)

# curated()
for sam in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
    writeFasta(sam, 'ribo', 15, 15, 50, 0)
    writeFasta(sam, 'rna', 15, 15, 50, 0)
    print(sam)
    curated()

### remove repeats in X_train and X_test
# X1 = np.vstack((X_train,X_test))
# X1[:, 0] = np.append(y_train, y_test)
#
# unique_rows, indices, occurrence_count = np.unique(
#     X1, axis = 0, return_counts = True, return_index = True)
#
# uX_test = unique_rows[indices > X_train.shape[0]-1]
# # print(uX_test[:, 0])
# y_test =  np.copy(uX_test[:, 0])
# uX_test[:, 0] = 0
# print(X_train.shape)
# print(uX_test.shape)
