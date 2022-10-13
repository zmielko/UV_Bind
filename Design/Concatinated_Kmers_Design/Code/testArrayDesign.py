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
    ARRAY_DESIGN_FILE = "../Array_Design_Full.txt"
    POS_PROBES_FILE = "../Sequence_Designs/Positional_Probes.csv"
    ARRAY_INFO_LIST = ["ArrayInfo_085438.txt", "ArrayInfo_085400.txt", "ArrayInfo_084702.txt"]
    PRIMER = "CGCACACATACACATACACACATACACATACGCAC"
    MAX_LEN = 175906
    HAMMING_DIST_SEQ = "GGAAATTCCC"
    # Tests for the array design
    uvArrayDF = pd.read_csv(ARRAY_DESIGN_FILE, sep = '\t', header = None)
    uvArrayDF = uvArrayDF.rename(columns={1:"Sequence", 2:"Description"})
    for seq, desc in zip(uvArrayDF["Sequence"], uvArrayDF["Description"]):
        testProbeLenMatchDesc(seq, desc)
        testProbePrimer(seq, PRIMER)
        testMultiDNDescription(seq, desc)
        testFlankingProbe(seq, desc)
    testPosProbes(uvArrayDF, POS_PROBES_FILE)
    testDeBruijnProbes(uvArrayDF, ARRAY_INFO_LIST)
    testClusterDescription(uvArrayDF)
    testClusterDistanceDescription(uvArrayDF, HAMMING_DIST_SEQ)
    if len(uvArrayDF) > MAX_LEN:
        raise ValueError(f"Too many probes on the array.\nCurrently: {len(uvArrayDF)}\nMax: {MAX_LEN}")
    print("Array design passed all tests.")


# In[ ]:




