#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import numpy as np
import os
from Bio.Seq import reverse_complement

def windowSearchKmer(seq, k):
    """
    Given a sequence and a length k, returns a set of all kmers within that
        sequence.
    """
    kSet = set()
    for i in range(len(seq) - k + 1):
        kSet.add(seq[i:i+k])
    return(kSet)

def extRight(seq, k, prefSet, genSeqList, seenSet):
    """
    Given:
        1. Sequence
        2. k
        3. Preference Set
        4. Set of kmers already used
        5. List to append generated sequences to
    Recursively extends a sequence to the right if the next nucleotide makes
        a kmer in the preference set and is not in the seen set. Otherwise, 
        appends to the genSeqList.
    """
    queryStart = len(seq) - k + 1
    querySeq = seq[queryStart:]
    canExt = False
    for i in ['A', 'C', 'G', 'T']:
        if querySeq + i in prefSet and querySeq + i not in seenSet:
            canExt = True
            seenSet.add(querySeq + i)
            extRight(seq + i, k, prefSet, genSeqList, seenSet = seenSet)
    if canExt == False:
        genSeqList.append(seq)
        
def extLeft(seq, k, prefSet, genSeqLeftList, seenSet):
    """
    Given:
        1. Sequence
        2. k
        3. Preference Set
        4. Set of kmers already used
        5. List to append generated sequences to
    Recursively extends a sequence to the left if the next nucleotide makes
        a kmer in the preference set and is not in the seen set. Otherwise, 
        appends to the genSeqLeftList.
    """
    querySeq = seq[:k-1]
    canExt = False
    for i in ['A', 'C', 'G', 'T']:
        if i + querySeq in prefSet and i + querySeq not in seenSet:
            canExt = True
            seenSet.add(i + querySeq)
            extLeft(i + seq, k, prefSet, genSeqLeftList, seenSet = seenSet)
    if canExt == False:
        genSeqLeftList.append(seq)

def extendKmers(kmerList, prefSet, k=6):
    """
    Given a list of kmers, generates probes that are extended as far as possible using
        the extRight and extLeft methods. Returns all unique generated probes as a list.
    """
    totalProbes = set()
    # Extend left then right
    for i in kmerList:
        genSeqList = []
        seenSet = set([i])
        # Extend sequence to the left
        extLeft(i, k, prefSet = prefSet, genSeqLeftList = genSeqList, seenSet = seenSet)
        # Extend sequences to the right
        genSeqRightList = []
        for j in genSeqList:
            extRight(j, k, prefSet = prefSet, genSeqList=genSeqRightList, seenSet = seenSet)
        # Update the probe set
        totalProbes.update(genSeqRightList)   
    # Extend right then left
    for i in kmerList:
        genSeqList = []
        seenSet = set([i])
        # Extend sequence to the right
        extRight(i, k, prefSet = prefSet, genSeqList = genSeqList, seenSet = seenSet)
        # Extend sequences to the left
        genSeqLeftList = []
        for j in genSeqList:
            extLeft(j, k, prefSet = prefSet, genSeqLeftList=genSeqLeftList, seenSet = seenSet)
        # Update the probe set
        totalProbes.update(genSeqLeftList)
    totalProbesList = sorted(list(totalProbes))
    return(totalProbesList)       
        
def extCaniR(seq, kSet,kDictNonSpecDiff, k=6):
    """
    Helper function for extend25.
    Given:
        1. Sequence to extend
        2. Set of kmers to use
        3. Dictionary of differences in kmer E-scores
        4. k
    Looks to extend by A,C,G, or T to the right. The one with the greatest difference
        is chosen and the sequence with the additional nucelotide is returned.
    """
    queryStart = len(seq) - k + 1
    querySeq = seq[queryStart:]
    bestCani = ''
    bestCaniPref = -0.5
    for i in ['A', 'C', 'G', 'T']:
        if querySeq + i in kSet:
            if kDictNonSpecDiff[querySeq + i] > bestCaniPref:
                bestCani = i
                bestCaniPref = kDictNonSpecDiff[querySeq + i]
    return(seq + bestCani)

def extCaniL(seq, kSet,kDictNonSpecDiff, k=6):
    """
    Helper function for extend25.
    Given:
        1. Sequence to extend
        2. Set of kmers to use
        3. Dictionary of differences in kmer E-scores
        4. k
    Looks to extend by A,C,G, or T to the left. The one with the greatest difference
        is chosen and the sequence with the additional nucelotide is returned.
    """
    querySeq = seq[:k - 1]
    bestCani = ''
    bestCaniPref = -0.5
    for i in ['A', 'C', 'G', 'T']:
        if i + querySeq in kSet:
            if kDictNonSpecDiff[i + querySeq] > bestCaniPref:
                bestCani = i
                bestCaniPref = kDictNonSpecDiff[i + querySeq]
    return(bestCani + seq)

def extend25(seq, kSet,kDictNonSpecDiff):
    """
    Greedy extension of a sequence. 
    Given a sequence to extend, a set of kmers, and a dictionary with the differences
        between conditions.
    Iterates by trying to extend left, then right. Terminates when the probe length is
        25 or if neither extCaniL or extCaniR can add an additional nucelotide.
    """
    extProbe = seq
    prevLen = len(extProbe)
    while True:
        extProbe = extCaniL(extProbe, kSet,kDictNonSpecDiff)
        if len(extProbe) == 25:
            return(extProbe)
        extProbe = extCaniR(extProbe, kSet,kDictNonSpecDiff)
        if len(extProbe) == 25:
            return(extProbe)
        if len(extProbe) == prevLen:
            return("N")
        prevLen = len(extProbe)

def all25mersLongProbes(totalProbesList):
    """
    Given a list of probes generated by extendKmers
    For probes that are 25bp or more, add all 25mer substrings to a set and 
        return it as a list. 
    """
    all25mer = set()  
    for i in totalProbesList:
        iLen = len(i)
        if iLen >= 25:
            all25mer.update(windowSearchKmer(i, 25))  
    return(sorted(list(all25mer)))

def extendShortProbes(totalProbesList, nonSpecWCSet,kDictNonSpecDiff):
    """
    Given a list of probes generated by extendKmers
    Those that are less than 25bp are extended with the extend25 function. 
    If a probe can be extended to 25bp, it is added to a set. 
    Returns the set of 25bp extended probes.
    """
    shortExtProbeList = []
    for probe in totalProbesList:
        # If the probe has been extended or is less than 25bp
        if len(probe) > 6 and len(probe) < 25:
            #Extend
            extSeq = extend25(probe, nonSpecWCSet,kDictNonSpecDiff)
            # If it was extended to 25bp, append to list
            if extSeq != "N":
                shortExtProbeList.append(extSeq)
    # Remove duplicate probes in the list
    shortExtProbeList = sorted(list(set(shortExtProbeList)))
    return(shortExtProbeList)

def filter_one_orientation(probe_list):
    """
    Given a probe list, ensures there are no duplicates based on 
        orientation. 
    """
    already_added = set()
    probe_rc_list = [reverse_complement(probe) for probe in probe_list]
    return_list = []
    for i in probe_list:
        if i not in already_added:
            return_list.append(i)
            already_added.add(i)
            probe_rc = reverse_complement(i)
            for j in probe_list:
                if probe_rc == j:
                    already_added.add(j)
    return(return_list)
    
if __name__ == '__main__':
    fileTup = [(f"EGR1_2_3/OLS_EGR1_2_3_Escore.txt", "EGR1"),
               (f"CREB1_0_1/OLS_CREB1_0_1_Escore.txt", "CREB1"),
               (f"TBP_28_29/OLS_TBP_28_29_Escore.txt", "TBP"),]
    dataFolder = "../Data"
    os.makedirs(f"{dataFolder}/Diff_UV_Probes/", exist_ok = True)
    for fT in range(len(fileTup)):
        ols = pd.read_csv(f"{dataFolder}/{fileTup[fT][0]}", sep = '\t', skiprows = 7)
        ols["Diff"] = ols["Escore_UV9"] - ols["Escore_WC8"]
        kDictScoreDiff = dict(zip(ols["kmerFwd"], ols["Diff"])) 
        kDictScoreDiff.update(dict(zip(ols["kmerRC"], ols["Diff"])))
        nonSpecWC = ols.query("Escore_WC8 < 0").copy(deep = True)
        nonSpecWC["Diff"] = nonSpecWC["Escore_UV9"] - nonSpecWC["Escore_WC8"]
        nonSpecWCSet = set(list(nonSpecWC["kmerFwd"]) + list(nonSpecWC["kmerRC"]))
        kDictNonSpecDiff = dict(zip(nonSpecWC["kmerFwd"], nonSpecWC["Diff"])) 
        kDictNonSpecDiff.update(dict(zip(nonSpecWC["kmerRC"], nonSpecWC["Diff"])))
        UVpref = ols.query("Preference_Score < 0")
        prefKmers = sorted(list(set(list(UVpref["kmerFwd"]) + list(UVpref["kmerRC"]))))
        prefSet = set(prefKmers)
        totalProbesList = extendKmers(prefKmers, prefSet)

        # Find all 25mers from 25bp+ probes, save to output
        all25merList = all25mersLongProbes(totalProbesList)
        all25merList = filter_one_orientation(all25merList)
        all25DF = pd.DataFrame({"Sequence":all25merList})
        all25DF = all25DF.sort_values(by = "Sequence")
        all25DF.to_csv(f"{dataFolder}/Diff_UV_Probes/{fileTup[fT][1]}_25mers.tsv", sep = '\t', header = None, index = False)
        
        # Extend probes shorter than 25bp to 25bp, save to output
        extShortProbeList = extendShortProbes(totalProbesList, nonSpecWCSet,kDictNonSpecDiff)
        extShortProbeList = filter_one_orientation(extShortProbeList)
        extShortProbesDF = pd.DataFrame({"Sequence":extShortProbeList})
        extShortProbesDF = extShortProbesDF.sort_values(by = "Sequence")
        extShortProbesDF.to_csv(f"{dataFolder}/Diff_UV_Probes/{fileTup[fT][1]}_extShort.tsv", sep = '\t', header = None, index = False)


# In[ ]:




