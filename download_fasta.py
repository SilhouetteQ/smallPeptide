#!/usr/bin/env python
# encoding: utf-8
import sys
import os

import matplotlib.pyplot as plt
import pandas as pd
from gtfparse import read_gtf
import pickle
import numpy as np
import requests, sys
import re
from Bio import SeqIO

import pylab
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable

print(os.getcwd())
def save_obj(obj, name):
    with open('../model/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('../model/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def annoENST():
    annoFILE = 'H:/Project/smallPeptide/data/gencode.v35.annotation.gtf'
    df = read_gtf(annoFILE)
    df_lnc = df[df["gene_type"] == "lncRNA"]
    df_lnc['transcript_id'].replace('', np.nan, inplace=True)
    df_lnc1 = df_lnc.loc[df_lnc['transcript_id'].drop_duplicates().index]
    df_lnc2 = df_lnc1.loc[df_lnc1['transcript_id'].dropna().index]
    save_obj(df_lnc2[['gene_id', 'transcript_id']], 'gencode.v35.annotation')
# annoENST()

def readCount(sam, type):
    countPath = '../../../data/'+ sam + '/' + type + '/counts.txt'
    df = pd.read_table(countPath, header= None)
    df.drop(df.tail(5).index,inplace=True)
    return df[df[1] > 0][0]

def readGTF(sam, type):
    countPath = '../../../data/'+ sam + '/' + type + '/final_count.gtf'
    df = read_gtf(countPath, column_converters={"TPM": float, "FPKM": float})
    save_obj(df, sam+"count")
    # df = load_obj(sam+"count")
    ### set TPM>0.2 for RNA-Seq threshold
    # a = df[df.TPM > 0.2]
    # annoLNC = load_obj('gencode.v35.annotation')
    # a['reference_id'].replace('', np.nan, inplace=True)
    # a1 = a.loc[a['reference_id'].dropna().index]
    # af = a1[a1.reference_id.isin(annoLNC.transcript_id)]
    # return af['reference_id']
    return df

def saveSamDict():
    samRiboTranscript = {}
    sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
    for count, name in enumerate(sam):
        a = readGTF(name, 'ribo')
        samRiboTranscript[name] = a
        print(name, ' ribo ', len(a))
    save_obj(samRiboTranscript, 'riboTranscript')
    samRNATranscript = {}
    for count, name in enumerate(sam):
        a = readGTF(name, 'rna')
        samRNATranscript[name] = a
        print(name, ' rna ', len(a))
    save_obj(samRNATranscript, 'rnaTranscript')
# saveSamDict()

def rnaRiboOverlap():
    rna = load_obj('rnaTranscript')
    ribo = load_obj('riboTranscript')
    seqDict = load_obj('threeNTFasta')
    rnt = rna['breast'].append(rna['fibroblast']).append(rna['HEK293_2']).append(rna['liver'])
    rnt = [item[:15] for item in rnt]
    rit = ribo['breast'].append(ribo['fibroblast']).append(ribo['HEK293_2']).append(ribo['liver'])
    rit = [item[:15] for item in rit]
    dID = [value for value in rnt if value in rit]

# from matplotlib_venn import venn2, venn3
# venn2((len(set(rit)), len(set(rnt)), len(set(dID))), ['Ribo-Seq', 'RNA-Seq'])
# plt.show()

alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}
def venn4(data=None, names=None, output=None, fill="number", show_names=True, show_plot=True, **kwds):
    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize)   # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, labels['1000'], **alignment)
    pylab.text(280, 200, labels['0100'], **alignment)
    pylab.text(155, 250, labels['0010'], **alignment)
    pylab.text(245, 250, labels['0001'], **alignment)
    # 2
    pylab.text(200, 115, labels['1100'], **alignment)
    pylab.text(140, 225, labels['1010'], **alignment)
    pylab.text(145, 155, labels['1001'], **alignment)
    pylab.text(255, 155, labels['0110'], **alignment)
    pylab.text(260, 225, labels['0101'], **alignment)
    pylab.text(200, 240, labels['0011'], **alignment)
    # 3
    pylab.text(235, 205, labels['0111'], **alignment)
    pylab.text(165, 205, labels['1011'], **alignment)
    pylab.text(225, 135, labels['1101'], **alignment)
    pylab.text(175, 135, labels['1110'], **alignment)
    # 4
    pylab.text(200, 175, labels['1111'], **alignment)
    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    pylab.savefig(output + '.png')
    if show_plot:
        pylab.show()


def get_labels(data, fill="number"):
    """
    to get a dict of labels for groups in data
    input
      data: data to get label for
      fill = ["number"|"logic"|"both"], fill with number, logic label, or both
    return
      labels: a dict of labels for different sets
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
    Out[12]:
    {'001': '001: 0',
     '010': '010: 5',
     '011': '011: 0',
     '100': '100: 3',
     '101': '101: 2',
     '110': '110: 2',
     '111': '111: 3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    if fill == "number":
        labels = {k: len(set_collections[k]) for k in set_collections}
    elif fill == "logic":
        labels = {k: k for k in set_collections}
    elif fill == "both":
        labels = {k: ("%s: %d" % (k, len(set_collections[k]))) for k in set_collections}
    else:  # invalid value
        raise Exception("invalid value for fill")

    return labels

def rnaDownload(sam):
    a = load_obj('rnaTranscript')
    seqID = []
    dID = []
    for b in sam:
        seqID += a[b].to_list()
    dID = [item[:15] for item in seqID]
    # print(len(set(dID)))
    seqDict = load_obj('rnaCountFasta')

    # missing = load_obj('missingrnaCountFasta')
    # print(len(seqDict.values()))
    # print(missing)
    # seqDict = {}
    # missing = []
    # server = "https://rest.ensembl.org"
    # for seq in set(dID):
    #     if seq in seqDict.keys():
    #         continue
    #     ext = "/sequence/id/" + seq + "?type=cdna"
    #     r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    #     if not r.ok:
    #         print(seq)
    #         missing.append(seq)
    #         continue
    #     seqDict[seq] = r.text
    #     save_obj(seqDict, 'rnaCountFasta')
    #     save_obj(missing, 'missingrnaCountFasta')
# rnaDownload(['breast', 'fibroblast', 'HEK293_2', 'liver'])
# seqDict = load_obj('rnaCountFasta')
# missing = load_obj('missingrnaCountFasta')

def riboCountLength():
    sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
    a = load_obj('riboTranscript')
    for b in sam:
        print(len(a[b]))
# riboCountLength()

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
# threeNTCount('liver')

### construct a fasta dictionary

def threeNTDownload(sam):
    dID = []
    seqDict = {}
    missing = []
    for s in sam:
        seqInfo = threeNTCount(s)
        seqID = seqInfo['transcript_id'].to_list()
        dID += [item[:15] for item in seqID]
    # print(len(dID))
    # print(len(set(dID)))
    # seqInfo = threeNTCount(sam)
    # seqID = seqInfo['transcript_id'].to_list()
    # dID += [item[:15] for item in seqID]
    seqDict = load_obj('threeNTFasta')
    missing = load_obj('missingthreeNTFasta')
    # print(len(set(seqDict.keys())))
    # server = "https://rest.ensembl.org"
    # for seq in set(dID):
    #     if seq in seqDict.keys():
    #         continue
    #     ext = "/sequence/id/" + seq + "?type=cdna"
    #     r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    #     if not r.ok:
    #         print(seq)
    #         missing.append(seq)
    #         continue
    #     seqDict[seq] = r.text
    #     save_obj(seqDict, 'threeNTFasta')
    #     save_obj(missing, 'missingthreeNTFasta')
# threeNTDownload('liver')
# threeNTDownload(['breast', 'fibroblast', 'HEK293_2', 'liver'])
# empty seq = ['ENST00000675070', 'ENST00000408966', 'ENST00000536039', 'ENST00000235290']

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

def riboStat(sam, type, cutRank, orfCutoff, cutStart, cutSeq1, cutSeq, t_id1, t_id):
    # seqpath = 'riboStat.txt'
    seqpath = type + 'Stat.txt'
    f = open(seqpath, 'a')
    f.write(sam)
    f.write('\t')
    f.write(str(len(set(cutSeq1))))
    f.write('\t')
    f.write(str(len(set(cutSeq))))
    f.write('\t')
    f.write(str(len(set(t_id1))))
    f.write('\t')
    f.write(str(len(set(t_id))))
    f.write('\t')
    f.write(str(cutRank))
    f.write('\t')
    f.write(str(orfCutoff))
    f.write('\t')
    f.write(str(cutStart))
    f.write('\n')
    f.close()

def seqUpDown(sam, type, cutLen, cutStart, orfCutoff, cutRank):
    if type == 'ribo':
        seqDict = load_obj('threeNTFasta')
        seqInfo = threeNTCount(sam)
        # for a in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
        #     seqInfo = threeNTCount(a)
        #     print(a)
        #     print(seqInfo.shape)
        seqInfo['t_id'] = [item[:15] for item in seqInfo['transcript_id']]
        cutSeq = []
        cutSeq1 = []
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
            if len(qualify) < cutRank+1:
                continue
            t_id1.append(id)
            cutSeq1.append(re2[(sindex[cutRank] - cutLen - 1):(sindex[cutRank] + cutLen + 2)])
            # print(qualify[cutRank], cut)
            if sindex[cutRank] == cut:
                # f.write(seqInfo['t_id'].iloc[i] + ' ' + str(seqInfo['ORF_tstart'].iloc[i]) + ' ' + str(len(re2)))
                # f.write('\n')
                # f.write(re2[(cut - cutLen - 1):(cut + cutLen + 2)])
                # f.write('\n')
                cutSeq.append(re2[(cut - cutLen - 1):(cut + cutLen + 2)])
                t_id.append(id)
                wholeLen.append(seq)
                # print(re2)
                # print(re2[(cut - cutLen - 1):(cut + cutLen + 2)])
        # riboStat(sam, type, cutRank, orfCutoff, cutStart,
        #          cutSeq1, cutSeq, t_id1, t_id)
        print(sam)
        print('All possible cut seq: ', len(cutSeq1))
        print('All possible cut seq on site ', cutRank, ' : ', len(cutSeq))
        print('All possible unique cut seq: ', len(set(cutSeq1)))
        print('Unique cut seq on site ', cutRank, ' : ', len(set(cutSeq)))
        print('All possible unique transcript ID: ', len(set(t_id1)))
        print('Unique transcript ID on site ', cutRank, ' : ', len(set(t_id)))
        # return len(set(t_id1)), len(set(t_id)), len(set(cutSeq1)), len(set(cutSeq))
        # return t_id
        save_obj(set(cutSeq), sam + 'CUT' + type)
        # save_obj(set(wholeLen), sam + 'CUT' + type)
        # return wholeLen, t_id
        return cutSeq, t_id
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
        cutSeq = []
        cutSeq1 = []
        wholeLen = []
        t_id = []
        t_id1 = []
        for i, id in enumerate(dID):
            if id in emptySeq:
                continue
            seq = seqDict[id]
            re2 = ''.join(re.split('\n', seq)[1:])
            sindex, qualify = findORF(re2, orfCutoff, cutStart)
            if len(qualify) < cutRank + 1:
                continue
            cutSeq.append(re2[(sindex[cutRank] - cutLen - 1):(sindex[cutRank] + cutLen + 2)])
            wholeLen.append(seq)
            t_id.append(id)
        save_obj(set(cutSeq), sam + 'CUT' + type)
        # save_obj(set(cutSeq), sam + 'sCUT' + type)
        # save_obj(set(wholeLen), sam + 'CUT' + type)
        return cutSeq, t_id
        # return wholeLen, t_id
        # riboStat(sam, type, cutRank, orfCutoff, cutStart, cutSeq1, cutSeq, t_id1, t_id)

# seqUpDown('liver', 'ribo', 15, 15, 50, 0)
# seqUpDown('liver', 'rna', 15, 15, 50, 0)

def writeFasta(sam, type, cutLen, cutStart, orfCutoff, cutRank, enst=True, rep = True):
    cutSeq, t_id = seqUpDown(sam, type, cutLen, cutStart, orfCutoff, cutRank)
    if enst is False:
        if rep is False:
            seqPath = './fasta/' + sam + type + str(cutLen)+ '_seqonly_Uniq.fa'
        else:
            seqPath = './fasta/' + sam + type + str(cutLen) + '_seqonly_Dupli.fa'
        f = open(seqPath, 'w')
        for item in set(cutSeq):
            f.write(item)
            f.write('\n')
        f.close()
    else:
        if rep is False:
            forMotif = './fasta/' + sam + type + str(cutLen)+ '_enst_Uniq.fa'
        else:
            forMotif = './fasta/' + sam + type + str(cutLen)+ '_enst_Dupli.fa'
        f1 = open(forMotif, 'w')
        rep = []
        for m, item in enumerate(cutSeq):
            if item in rep:
                pass
            else:
                # rep.append(item)
                # f1.write('>')
                # f1.write(t_id[m])
                # f1.write('\n')
                f1.write(item)
                # f1.write('\n')
        f1.close()


# t = 'breast'+'fibroblast'+'HEK293_2'+'liver'
# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# cutSeq0 = []
# t_id0 = []
# type = 'rna'
# cutLen = 30
# seqpath = './fasta/' + t + type + str(cutLen)+ '_enst_Uniq.fa'
# f1 = open(seqpath, 'w')
# for b in sam:
#     cutSeq, t_id = seqUpDown(b, type, cutLen, cutLen, 50, 0)
#     cutSeq0 += cutSeq
#     t_id0 += t_id
# rep = []
# for m, item in enumerate(cutSeq0):
#     if item in rep:
#         pass
#     else:
#         rep.append(item)
#         f1.write('>')
#         f1.write(t_id0[m])
#         f1.write('\n')
#         f1.write(item)
#         f1.write('\n')

# t = 'breast'+'fibroblast'+'HEK293_2'+'liver'

# writeFasta('liver', 'ribo', 15, 15, 50, 0, True, False)
# writeFasta('liver', 'rna', 15, 15, 50, 0, True, False)


def intersecT(sam):
    t1 = []
    for sam in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
        # globals()[sam] = list(seqUpDown(sam, 'ribo', 15, 15, 50, 0))
        globals()[sam] = list(seqUpDown(sam, 'rna', 15, 15, 50, 0))
        t1 += sam
    # venn4([set(breast), set(fibroblast),set(HEK293_2), set(liver)],
    #       ['breast', 'fibroblast', 'HEK293_2', 'liver'], 'All possible cut seq on site 0(RNA)')

# intersecT(['breast', 'fibroblast', 'HEK293_2', 'liver'])
# def intersecT(sam):
#     t1 = []
#     for s in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
#         print(locals())
#         # print(s)
#         locals()[s] = list(seqUpDown(s, 'ribo', 15, 15, 50, 0))
#         t1 += s
#     venn4([set(breast), set(fibroblast),set(HEK293_2), set(liver)],
#           ['breast', 'fibroblast', 'HEK293_2', 'liver'], 'Intersect of transcripts(15)')
# intersecT(['breast', 'fibroblast', 'HEK293_2', 'liver'])

# sys.exit()
# save_obj(t, sam + 'CUT' + type)

# seqUpDown('fibroblast', 'ribo', 15, 15, 50, 0)
# seqUpDown('fibroblast', 'rna', 15, 15, 50, 0)
# seqUpDown('breast', 'ribo', 15, 15, 50, 0)
# seqUpDown('breast', 'rna', 15, 15, 50, 0)
# seqUpDown('HEK293_2', 'ribo', 15, 15, 50, 0)
# seqUpDown('HEK293_2', 'rna', 15, 15, 50, 0)
# seqUpDown('liver', 'ribo', 15, 15, 50, 0)
# seqUpDown('liver', 'rna', 15, 15, 50, 0)

# print(len(a))
# print(len(b))
# t = [value for value in list(a) if value  in list(b)]
# tt = [value for value in t if value  in list(c)]
# ttt = [value for value in tt if value  in list(d)]
# print(ttt)
# seqPath = './Try.fa'
# f = open(seqPath, 'w')
# for item in ttt:
#     f.write(item)
#     f.write('\n')
# f.close()

def seq2ngram(cutSeq, k):
    kmer = ''
    for num, line in enumerate(cutSeq):
        line = line.lower()
        l = len(line)
        s = 1
        for i in range(0, l, s):
            if i+k >= l+1:
                break
            kTmp = ''.join(line[i:i+k])
            kmer += kTmp
            kmer += ' '
        kmer += '\n'
    return kmer

def multiSample(*sam):
    pos = []
    neg = []
    ori = ''
    for i in sam:
        ribo = load_obj(i + 'CUT' + 'ribo')
        # print(ribo)
        # sys.exit()
        pos += ribo
        rna = load_obj(i + 'CUT' + 'rna')
        neg += rna
        ori += i
    save_obj(pos, ori + 'CUT' + 'ribo')
    save_obj(neg, ori + 'CUT' + 'rna')
# multiSample('breast', 'fibroblast', 'HEK293_2', 'liver')

def updatecutSeq(cutLen, cutRank):
    sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
    for i, b in enumerate(sam):
        seqUpDown(b, 'ribo', cutLen, cutLen, 50, cutRank)
        seqUpDown(b, 'rna', cutLen, cutLen, 50, cutRank)
# updatecutSeq(100, 0)

def prepareData(sam, k, cutLen, kmer_index=None):
    from keras.preprocessing.text import Tokenizer
    from keras.preprocessing.sequence import pad_sequences
    pos = load_obj(sam + 'CUT' + 'ribo')
    neg = load_obj(sam + 'CUT' + 'rna')
    overlap = [value for value in list(pos) if value in list(neg)]
    overlap = set(overlap)
    if len(overlap) > 0:
        print('overwhelming ' + str(len(overlap)) + ' overlaps')
        print(overlap)
        # print(pos)
        # print(neg)
        for item in overlap:
            # pos.remove(str(item))
            neg.remove(str(item))
    print('need to n-gram %d seqs' % len(pos))
    print('need to n-gram %d seqs' % len(neg))
    posK = seq2ngram(pos, 5)
    negK = seq2ngram(neg, 5)
    pos_seqs = [line[:-1] for line in posK.splitlines()]
    neg_seqs = [line[:-1] for line in negK.splitlines()]
    seqs = pos_seqs + neg_seqs
    if kmer_index is None:
        tokenizer = Tokenizer(num_words=20000)
        tokenizer.fit_on_texts(seqs)
        sequences = tokenizer.texts_to_sequences(seqs)
        kmer_index = tokenizer.word_index
        save_obj(kmer_index, 'kmer_index')
        MAX_LEN = cutLen * 3
        X = pad_sequences(sequences, maxlen=MAX_LEN)
        y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
    else:
        sequences = []
        fail = []
        for t, line in enumerate(seqs):
            l = []
            for i in range(len(line.split(' '))):
                if line.split(' ')[i] not in kmer_index.keys():
                    fail.append(t)
                    break
                l.append(kmer_index[line.split(' ')[i]])
            sequences.append(l)
            # print(line)
            # print(line.split(' '))
            # print(kmer_index[line.split(' ')[0]])
        for ele in sorted(fail, reverse=True):
            del sequences[ele]
        MAX_LEN = cutLen * 3
        X = pad_sequences(sequences, maxlen=MAX_LEN)
        y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
        for ele in sorted(fail, reverse=True):
            y = np.delete(y, ele)
    return X, y
# prepareData('liver', 5, 15)

def multidim_intersect(arr1, arr2):
    arr1_view = arr1.view([('', arr1.dtype)] * arr1.shape[1])
    arr2_view = arr2.view([('', arr2.dtype)] * arr2.shape[1])
    intersected = np.intersect1d(arr1_view, arr2_view)
    return intersected.view(arr1.dtype).reshape(-1, arr1.shape[1])

def create_classifier(c, gamma, verbose=False, shrinking=True, probability=True):
    return svm.SVC(kernel='rbf', C=c, gamma=gamma, decision_function_shape='ovr', max_iter=-1, tol=0.001,
                   verbose=verbose, shrinking=shrinking, probability=probability)


def mlTraining(sam, cutLen = 15):
    from sklearn import metrics, svm
    from sklearn.model_selection import train_test_split, GridSearchCV
    from sklearn.svm import SVC
    from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, f1_score
    from sklearn.metrics.cluster import contingency_matrix
    # grid search k value for SMOTE oversampling for imbalanced classification

    import sys
    sys.path.append("../../../../lncRNAIdentification")
    # from sequence_attributes_orf import SequenceAttributes
    from sequence_attributes_kmer import SequenceAttributes
    pos = SequenceAttributes('./fasta/fibroblastribo15_enst_Uniq.fa', 1)
    neg = SequenceAttributes('./fasta/fibroblastrna15_enst_Uniq.fa', 0)
    # pos = SequenceAttributes('./fasta/breastfibroblastHEK293_2liverribo30_enst_Uniq.fa', 1)
    # neg = SequenceAttributes('./fasta/breastfibroblastHEK293_2liverrna30_enst_Uniq.fa', 0)
    pos_df = pd.DataFrame(pos.process())
    neg_df = pd.DataFrame(neg.process())
    print(pos_df.shape)
    print(neg_df.shape)
    # df.to_csv('trainingfile.txt', index=False)
    # X = pd.concat([pos_df, neg_df.iloc[:pos_df.shape[0]]])
    X = pd.concat([pos_df, neg_df])

    y = X['class']
    # y = np.array([1] * len(pos_df) + [0] * len(neg_df))
    X = X.drop(columns=['id','class'])
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

    X =  np.array(X)
    y = np.array(y)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size = 0.33, random_state = 42)

    # from sklearn.model_selection import StratifiedKFold
    # kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

    # from numpy import mean
    # from sklearn.model_selection import cross_val_score
    # from sklearn.model_selection import RepeatedStratifiedKFold
    # from imblearn.ensemble import BalancedBaggingClassifier
    # from sklearn.ensemble import RandomForestClassifier
    # from imblearn.ensemble import BalancedRandomForestClassifier
    # from imblearn.ensemble import EasyEnsembleClassifier
    # model = EasyEnsembleClassifier(n_estimators=5)
    # model = BalancedBaggingClassifier()
    # model = RandomForestClassifier(n_estimators=5)
    # define evaluation procedure
    # cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=1)
    # # evaluate model
    # scores = cross_val_score(model, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
    # print('Mean ROC AUC: %.3f' % mean(scores))

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
    #
    # model = GridSearchCV(estimator=SVC(),
    #                      param_grid= {'kernel': ['rbf'], 'C': [1, 10, 100, 1000],
    #               'gamma': [0.1, 0.01, 0.001, 0.0001]},
    #                      cv = kfold).fit(X_train, y_train)
    # y_pred = model.predict(X_test)
    # print(accuracy_score(y_test, y_pred))
    # # print(classification_report(y_test, y_pred))
    # print(contingency_matrix(y_test, y_pred))
    # print(roc_auc_score(y_test, y_pred))
    # print(f1_score(y_test, y_pred))
    # from xgb_attributes import XGBAttributes
    # xgbselect = XGBAttributes(X_train, 'svm')
    # f = xgbselect.attributes()
    # print(xgbselect)
    # print(f)


# mlTraining('fibroblast', 15)


# t = 'breast'+'fibroblast'+'HEK293_2'+'liver'
# mlTraining(t, 15)

def crossValidation(sam1, sam2):
    X_train, y_train = prepareData(sam1, 5, 15)
    kmer_index = load_obj('kmer_index')
    X_test, y_test = prepareData(sam2, 5, 15, kmer_index)

    model = GridSearchCV(estimator=SVC(),
                         param_grid= {'kernel': ['rbf'], 'C': [1, 10, 100, 1000],
                  'gamma': [0.1, 0.01, 0.001, 0.0001]},
                         cv = 5).fit(X_train, y_train)
    y_pred = model.predict(X_test)
    print(roc_auc_score(y_test, y_pred))
    print(f1_score(y_test, y_pred))
    print(classification_report(y_test, y_pred))
    print(contingency_matrix(y_test, y_pred))
    return roc_auc_score(y_test, y_pred), f1_score(y_test, y_pred)

# crossValidation('liver', 'breast')

def listToString(s):
    str1 = ''
    return (str1.join(s))


# import copy
# seqpath = 'crossSampleValidate.txt'
# f = open(seqpath, 'a')
# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# for i, b in enumerate(sam):
#     roc, f1 = mlTraining(b)
#     f.write(b)
#     f.write('\t')
#     f.write(b)
#     f.write('\t')
#     f.write(str(f1))
#     f.write('\t')
#     f.write(str(roc))
#     f.write('\n')

#     tt = copy.copy(sam)
#     del tt[i]
#     train = listToString(tt)
#     roc, f1 = crossValidation(train, b)
#     f.write(train)
#     f.write('\t')
#     f.write(b)
#     f.write('\t')
#     f.write(str(f1))
#     f.write('\t')
#     f.write(str(roc))
#     f.write('\n')

#     c = b + '2'
#     roc, f1 = crossValidation(b, c)
#     f.write(b)
#     f.write('\t')
#     f.write(c)
#     f.write('\t')
#     f.write(str(f1))
#     f.write('\t')
#     f.write(str(roc))
#     f.write('\n')
# f.close()


# from venn import venn
# category = {}
# # for sam in ['breast', 'fibroblast', 'HEK293_2', 'liver']:
# sam = 'liver'
# RiboR = ribo[sam]
# RiboR = set([item[:15] for item in RiboR])
# RNA = rna[sam]
# RNA = set([item[:15] for item in RNA])
# seqInfo = threeNTCount(sam)
# a = set([item[:15] for item in seqInfo['transcript_id']])
# category['Ribo'] = RiboR
# category['RNA'] = RNA
# category['3NT'] = a
# venn(category)
# plt.show()
    # venn3(subsets = (len(), 10, 12, 10, 9, 4, 3))


# (sam, type, cutLen, cutStart, orfCutoff, cutRank)
# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# cutVar = 'cutVar.txt'
# f = open(cutVar, 'a')
# f.write('\n')
# ml_result = 'ml_result.txt'
# g = open(ml_result, 'a')
# for b in sam:
#     f.write(b+'CutLength')
#     f.write('\t')
#     f.write('All tID')
#     f.write('\t')
#     f.write('First tID')
#     f.write('\t')
#     f.write('All cutSeq')
#     f.write('\t')
#     f.write('First cutSeq')
#     f.write('\t')
#     f.write('First/All tTD')
#     f.write('\t')
#     f.write('First/All cutSeq')
#     f.write('\t')
#     f.write('First cutSeq/tID')
#     f.write('\n')
#     # for i in list(range(1, 210, 6)):
#     # for i in list(range(1, 50, 1)):
#     # for i in list(range(1, 1200, 12)):
#     for i in list(range(1, 15, 1)):
#         tn, ftn, cn, fcn = seqUpDown(b, 'ribo', i, i, 50, 0)
#         f.write(str(i))
#         f.write('\t')
#         f.write(str(tn))
#         f.write('\t')
#         f.write(str(ftn))
#         f.write('\t')
#         f.write(str(cn))
#         f.write('\t')
#         f.write(str(fcn))
#         f.write('\t')
#         f.write(str(ftn/tn))
#         f.write('\t')
#         f.write(str(fcn/cn))
#         f.write('\t')
#         f.write(str(fcn / ftn))
#         f.write('\n')
# f.close()
##     cutLen = 100
#     g.write(b)
#     g.write('\troc\tf1\n')
#     for i in list(range(20, 41, 1)):
#         seqUpDown(b, 'ribo', i, i, 50, 0)
#         seqUpDown(b, 'rna', i, i, 50, 0)
#         g.write(str(i))
#         g.write('\t')
#         roc, f1 = mlTraining(b, i)
#         g.write(str(roc))
#         g.write('\t')
#         g.write(str(f1))
#         g.write('\n')
        # mlTraining(b+'s')
    # c = b + '1'
    # crossValidation(b, c)
# g.close()



# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# riboOri = []
# rnaOri = []
# for b in sam:
#     ribo = seqUpDown(b, 'ribo', 15, 15, 50, 0)
#     rna = seqUpDown(b, 'rna', 15, 15, 50, 0)
#     riboOri += ribo
#     rnaOri += rna




# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# riboOri = []
# sriboOri = []
# rnaOri = []
# srnaOri = []
# for b in sam:
#     ribo = seqUpDown(b, 'ribo', 15, 15, 50, 0)
#     riboOri.append(len(ribo))
#     sriboOri.append(len(set(ribo)))
#     rna = seqUpDown(b, 'rna', 15, 15, 50, 0)
#     rnaOri.append(len(rna))
#     srnaOri.append(len(set(rna)))
# print(riboOri, sriboOri, rnaOri, srnaOri)
# df = pd.DataFrame({'Ribo-Seq':riboOri, 'RiboUnique': sriboOri,
#                    'RNA-Seq': rnaOri, 'RNAUnique':srnaOri},
#                   index =['breast', 'fibroblast', 'HEK293', 'liver'])
# ax = df.plot.bar(rot=0, title = 'Comparison of repeated cut sequences between Ribo-Seq and RNA-Seq')
# plt.show()


#     df = readGTF(b, 'ribo')
#     # data_wide = df.pivot(columns='day',
#     #                      values='total_bill')
#     # df.pivot(values='TPM')
#     # sns.distplot(df['TPM'])
#     df.TPM.plot.density(figsize=(7, 7),
#                            linewidth=4, xlim=(10000,1e6),)
#     plt.show()
#     df.concat(df)
# print(df1.shape)
# print(df.shape)
#
# # plotting multiple density plot
# data_wide.plot.kde(figsize=(8, 6),
#                    linewidth=4)

# seqUpDown(b, 'ribo', 15, 15, 50, 0)



# save_obj(af['reference_id'], 'liverRiboTmp')
# a = load_obj('liverRiboTmp')
# print(a.shape)


### sample statistics recording
# seqpath = 'samStat.txt'
# f = open(seqpath, 'a')
# f.write('Numbers of unique transcript')
# f.write('\n')
# for s in sam:
#     seqInfo = threeNTCount(s)
#     seqInfo['t_id'] = [item[:15] for item in seqInfo['transcript_id']]
#     f.write(s)
#     f.write('\t')
#     f.write(str(len(set(seqInfo['t_id']))))
#     f.write('\n')
# f.close()

### write fasta for motif discovery
# sam = ['breast', 'fibroblast', 'HEK293_2', 'liver']
# for b in sam:
#     # for i in range(15, 31, 15):
#     i = 30
#     writeFasta(b, 'ribo', i, i, 50, 0, True, True)
#     writeFasta(b, 'ribo', i, i, 50, 0, True, False)
#     writeFasta(b, 'ribo', i, i, 50, 0, False, True)
#     writeFasta(b, 'ribo', i, i, 50, 0, False, False)
#     writeFasta(b, 'rna', i, i, 50, 0, True, True)
#     writeFasta(b, 'rna', i, i, 50, 0, True, False)
#     writeFasta(b, 'rna', i, i, 50, 0, False, True)
#     writeFasta(b, 'rna', i, i, 50, 0, False, False)

def smote(X, y):
    import imblearn
    from imblearn.pipeline import Pipeline
    from imblearn.over_sampling import SMOTE
    from imblearn.under_sampling import RandomUnderSampler
    from numpy import mean
    from sklearn.datasets import make_classification
    from sklearn.model_selection import cross_val_score
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn.tree import DecisionTreeClassifier
    k_values = [1, 2, 3, 4, 5, 6, 7]
    for k in k_values:
        # define pipeline
        model = DecisionTreeClassifier()
        over = SMOTE(sampling_strategy=0.5, k_neighbors=k)
        under = RandomUnderSampler(sampling_strategy=0.5)
        steps = [('over', over), ('under', under), ('model', model)]
        pipeline = Pipeline(steps=steps)
        # evaluate pipeline
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=1)

        scores = cross_val_score(pipeline, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
        score = mean(scores)
        print('> k=%d, Mean ROC AUC: %.3f' % (k, score))
    s = svm.SVC(kernel='linear')
    s.fit(X_train, y_train)
    # print(s.coef_)
    # counter = Counter(y)
    # print(counter)
    # return roc_auc_score(y_test, y_pred), f1_score(y_test, y_pred)