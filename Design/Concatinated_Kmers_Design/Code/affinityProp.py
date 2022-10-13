#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation
from jellyfish import levenshtein_distance
import argparse
from Bio.Seq import reverse_complement

def argumentParse():
    programDescription = 'Affinity Propagation of Probes'
    parser = argparse.ArgumentParser(description=programDescription)
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i','--input',required = True, type=str, help='Input file')
    req.add_argument('-oC','--outputCluster',required = True, type=str, help='Output cluster file as tsv')
    req.add_argument('-oM','--outputMatrix', required = True, type=str,help='Output matrix file as npz')
    req.add_argument('-c','--convIter', required = True, type=int,help='Convergence Iteration parameter')
    req.add_argument('-r', "--randState", type=int, help="Affinity Propagation random_state")
    args = parser.parse_args()
    return(args)

def min_levenshtein_distance(string1, string2):
    """
    Given 2 strings, returns the minimum levenstein distance assuming
        forward and reverse complements are equivalent
    """
    string2_rc = reverse_complement(string2)
    o1_o1 = levenshtein_distance(string1, string2)
    o1_o2 = levenshtein_distance(string1, string2_rc)
    min_distance = min([o1_o1, o1_o2])
    return(min_distance)
    
def levenshteinMatrix(stringArray):
    """
    Given an array of strings, returns the levenshtein distance for 
      each string to every other string as a matrix. 
    """
    listDistList = []
    for string1 in stringArray:
        distList = []
        for string2 in stringArray:
            distList.append(min_levenshtein_distance(string1, string2))
        listDistList.append(distList)
    distMatrix = np.array(listDistList)
    return(distMatrix)

def affinityPropagationParameterOutput(randInput,convIterInput, outputFile):
    """
    Given the input for the random state, convergence iteration,
    and output file, writes a comment indicating the randomstate used.
    If a parameter is None, writes Default. 
    """
    f = open(outputFile, mode = 'w')
    if randInput != None:
        f.write(f"#AP Random State: {randInput}\n")
    else:
        f.write(f"#AP RandomState: Default\n")
    if convIterInput != None:
        f.write(f"#AP Convergence Iteration: {convIterInput}\n")
    else:
        f.write(f"#AP Convergence Iteration: Default\n")
    f.close()

if __name__ == '__main__':
    # Parse arguments
    args = argumentParse()
    # Read input file, convert to array
    print(f"ARGS_INPUT = {args.input}")
    df = pd.read_csv(args.input, sep = '\t', header = None)
    # Construct matrix of levenstein distances, save
    lMatrix = levenshteinMatrix(df[0])
    np.savez_compressed(args.outputMatrix, lMatrix, delimiter='\t')
    # Convert all distances to negative values for a similarity matrix
    lMatrix = lMatrix * -1
    # Run Affinity Propagation
    aProp = AffinityPropagation(affinity="precomputed", random_state = args.randState,
                                convergence_iter = args.convIter, verbose = True)
    aProp.fit(lMatrix)
    labels = aProp.labels_
    # Save randState Information to file
    affinityPropagationParameterOutput(args.randState,args.convIter, args.outputCluster)
    # Add cluster labels to input DataFrame and save
    df["Cluster"] = labels
    df = df.rename(columns={0:"Sequence"})
    df = df.sort_values(by = "Cluster")
    df.to_csv(args.outputCluster, sep = '\t', index = False, mode = 'a')


# In[ ]:




