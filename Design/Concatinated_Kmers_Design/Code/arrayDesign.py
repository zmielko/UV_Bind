#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pandas as pd
from jellyfish import hamming_distance
from Bio.Seq import reverse_complement
from arrayDesignFunctions import *


if __name__ == '__main__':
    # Constants 
    DINUC_LIST = [''.join(p) for p in itertools.product(['A', 'C', 'G', 'T'], repeat=2)]
    TETRANUC_LIST = [''.join(p) for p in itertools.product(['A', 'C', 'G', 'T'], repeat=4)]
    PRIMER = "CGCACACATACACATACACACATACACATACGCAC"
    NAME_LIST = ["TBP", "EGR1", "CREB1"]
    KMER_LOC_LIST = ["../Data/TBP_28_29/TBP_28_29_kmerData.txt",
                   "../Data/EGR1_2_3/EGR1_2_3_kmerData.txt",
                   "../Data/CREB1_0_1/CREB1_0_1_kmerData.txt"]
                         
    KMER_CLUST_LOC_LIST = ["../Data/Diff_UV_Probes/TBP_25mers_cluster.tsv",
                       "../Data/Diff_UV_Probes/EGR1_25mers_cluster.tsv",
                       "../Data/Diff_UV_Probes/CREB1_25mers_cluster.tsv"]
                         
    SHORT_EXT_CLUST_LOC_LIST = ["../Data/Diff_UV_Probes/TBP_extShort_cluster.tsv",
                            "../Data/Diff_UV_Probes/EGR1_extShort_cluster.tsv",
                            "../Data/Diff_UV_Probes/CREB1_extShort_cluster.tsv"]     
    ARRAY_INFO_LIST = ["ArrayInfo_085438.txt", "ArrayInfo_085400.txt", "ArrayInfo_084702.txt"]
    DINUC_MUT_FILE = "../Sequence_Designs/Dinucleotide_Mutations.csv"
    FLANK_FILE = "../Sequence_Designs/Flanks.csv"
    DEBRUIJN_FILE = f"../Sequence_Designs/{ARRAY_INFO_LIST[0]}"
    KD_FILE = "../Sequence_Designs/kD.csv"
    POS_PROBES_FILE = "../Sequence_Designs/Positional_Probes.csv"
    RELA_SEQ = "GGAAATTCCC"
    # Generate Probe DataFrames
    dinucMutDF = createDiProbes(DINUC_MUT_FILE, PRIMER, DINUC_LIST)
    flankMutDF = createFlankProbes(FLANK_FILE, PRIMER, TETRANUC_LIST)
    deBruDF = extractDeBruijn(DEBRUIJN_FILE)
    kDDF = createKDProbes(KD_FILE, PRIMER)
    posProbeDF = createPositionalProbes(POS_PROBES_FILE, PRIMER)
    clustDiffDF = createClustDiffProbesFromLists(NAME_LIST,
                                                 KMER_LOC_LIST,
                                                 KMER_CLUST_LOC_LIST,
                                                 SHORT_EXT_CLUST_LOC_LIST,
                                                 PRIMER,
                                                 3)
    clustCloseDF = createClustCloseProbesFromLists(NAME_LIST,
                                                   KMER_LOC_LIST,
                                                   KMER_CLUST_LOC_LIST,
                                                   SHORT_EXT_CLUST_LOC_LIST,
                                                   RELA_SEQ,
                                                   50,
                                                   PRIMER,
                                                   3)
    # Combine pandas DataFrames of probes
    dfList = [dinucMutDF,flankMutDF,deBruDF,kDDF,posProbeDF,clustDiffDF,clustCloseDF]
    uvArrayDF = pd.concat(dfList)
    # Save output
    saveUniqueProbes(uvArrayDF)
    saveAgilentFormattedArray(uvArrayDF)


# In[ ]:




