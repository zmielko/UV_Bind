#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import itertools
import pandas as pd
from jellyfish import hamming_distance
from Bio.Seq import reverse_complement

def testDuplicates(df, colname, processName, dfName):
    """
    Given a pandas dataframe, column name, and name of the process, if there are duplicate
        entires, an error is raised and the duplicated entries are listed.
    """
    dupRows = df[df[colname].duplicated()]
    if len(dupRows) > 0:
        dupEntries = list(set(dupRows[colname]))
        raise ValueError((f"During {processName} in {dfName}, the following entires"
                         f"had duplicates:\n{dupEntries}"))

def testProbeLen(probe, reqLen, funName):
    """
    Tests if a probe is the required length, if not raises a value error. 
    """
    probeLen = len(probe)
    if probeLen != reqLen:
        raise ValueError(f"In function {funName}, probe length should be {reqLen}, len = {probeLen}")
    
    
def dinucMutation(seq, r, description, dinucList, cons = False):
    """
    Given a sequence, seq, generate all possible dinucleotide mutations
    1. seq = sequence to mutate
    2. r = range in sequence to mutate, given as a list or tuple
     ex/ r = [0, 5] or (4,9)
    3. dinucList = a precomputed list of all dinucleotides. 
    4. cons = consolodate. Argument to consolodate matching mutations.
    """
    # Generate all sequences with mutations at given ranges
    mutList, mutDescList = [], []
    for i in range(r[0], r[1]):
        for mut in dinucList:
            mutation = seq[:i] + mut + seq[i + 2:]
            testProbeLen(mutation, 25, "dinucMutation")
            mutDesc = f"DN_{description}_pos{i}->{mut}"
            mutList.append(mutation)
            mutDescList.append(mutDesc)
    mutDF = pd.DataFrame({"Sequence":mutList, "Description":mutDescList})
    # Filter out mutations that match the original
    mutDF = mutDF.query("Sequence != @seq").reset_index(drop=True)
    # Consolodate equivalent mutations, if cons argument is True
    if cons == True:
        mutList, descList  = [], []
        for mutSeq, group in mutDF.groupby(by="Sequence"):
            consDesc = '-'.join(group["Description"])
            mutList.append(mutSeq)
            descList.append(consDesc)
        mutDF = pd.DataFrame({"Sequence":mutList, "Description":descList})
    # Append the original sequence
    mutDF = mutDF.append({"Sequence":seq, "Description":f"DN_{description}_posW->WT"}, ignore_index = True)
    return(mutDF)

def diFlankMutation(seq, coreRange, description, tetranucList):
    """
    Given a sequence, seq, generate all possible dinucleotide flanks
    1. seq = sequence to mutate
    2. coreRange = range in sequence to mutate, given as a list or tuple
     ex/ r = [0, 5] or (4,9)
    3. dinucList = a precomputed list of all dinucleotides. 
    4. cons = consolodate. Argument to consolodate matching mutations.    
    """
    mutList, mutDescList = [], []
    cFlankLeft = seq[0:coreRange[0] - 2]
    cFlankRight = seq[coreRange[1] + 2:]
    core = seq[coreRange[0]:coreRange[1]]
    for tetramer in tetranucList:
        mut = cFlankLeft + tetramer[:2] + core + tetramer[2:] + cFlankRight
        testProbeLen(mut, 25, "diFlankMutation")
        mutList.append(mut)
        mutDescList.append(f"Flank_{description}_{tetramer}")
    return(pd.DataFrame({"Sequence":mutList, "Description":mutDescList}))


def addReps(df, nReps):
    """
    Adds replicates to a dataframe of sequences
    df = Pandas DataFrame of sequences to generate replicates of
     must have a "Sequence" and "Description" column
    nReps = number of replicates, given as an integer
    """
    seqList, desList = [], []
    for seq, des in zip(df["Sequence"], df["Description"]):
        for i in range(nReps):
            seqList.append(seq)
            desList.append(f"{des}_r{i}")
    return(pd.DataFrame({"Sequence":seqList, "Description":desList}))

def sumScore(seq, k, kScoreDict):
    """
    Given a sequence and a dictionary of scores for each kmer of length k,
        returns the sum of those scores for all substrings in the sequence
        of length k.
    """
    score = 0
    for i in range(len(seq) - k + 1):
        score = score + kScoreDict[seq[i:i+k]]
    return(score)

def scoreClusters(df, kDict):
    """
    Given a pandas dataframe of clustered sequences, applies the sumScore
        function and returns a dataframe with only the top scoring sequence
        per cluster.
    """
    df["Score"] = df["Sequence"].apply(lambda x: sumScore(x, 6, kDict))
    df = df.sort_values(by = "Score", ascending = False)
    df = df.drop_duplicates(subset="Cluster")
    return(df)

def closestSubHamming(string, query):
    """
    Given a string and a query string that is the length of the string or
        less, searches all substrings the length of the query. Returns the
        smallest hamming distance across all substrings. 
    """
    qLen = len(query)
    minHam = qLen
    for i in range(len(string) - qLen + 1):
        ham = hamming_distance(query, string[i:i+qLen])
        if ham < minHam:
            minHam = ham
    return(minHam)

def reverseComplementDF(df):
    """
    Given a dataframe with a sequence and description
     ensures the sequence is upper case and applies the reverse_complement 
     function to the sequence. Returns a pandas dataframe with the reverse 
     complement sequences.
    """
    dfOut = df.copy(deep = True)
    seqList, descList  = [], []
    for seq, desc in zip(dfOut["Sequence"], dfOut["Description"]):
        seqList.append(reverse_complement(seq.upper()))
        descList.append(f"{desc}_O2")
    dfOut["Sequence"] = seqList
    dfOut["Description"] = descList
    return(dfOut)

def removeRepString(string):
    """
    Given a string with replicate information at the end seperated by an underscore,
        returns a string without the replicate information.
    """
    strList = string.split("_")[:-1]
    return("_".join(strList))

def createDiProbes(inputFileLoc, primer, dinucList):
    """
    Given a sequence design file for dinucleotide mutations and primer sequence, 
    returns a pandas DataFrame of the probes. 
    """
    # Read file and ensure sequences are upper case
    designDiDF = pd.read_csv(inputFileLoc)
    testDuplicates(designDiDF, "Description", "createDiProbes", "designDiDF")
    designDiDF["Sequence"] = designDiDF["Sequence"].apply(lambda x: x.upper())
    # Generate Dataframes for each sequence, concatinate them
    diMutDFList = []
    for desc, seq, mutSubStr in zip(designDiDF["Description"],
                                    designDiDF["Sequence"],
                                    designDiDF["Mutation_Substring"]):
        rStart = seq.find(mutSubStr)
        rEnd = rStart + len(mutSubStr)
        diMutDFList.append(dinucMutation(seq, [rStart, rEnd],desc,dinucList, cons=True))
    diMutDF = pd.concat(diMutDFList)
    # Add Primer to the sequences
    diMutDF["Sequence"] = diMutDF["Sequence"].apply(lambda x: x + primer)
    # Make 20 replicates for each sequence
    diMutDF = addReps(diMutDF, 20)
    return(diMutDF)

def createFlankProbes(inputFileLoc, primer, tetramerList):
    """
    Given an input file of sequences with a core sequence defined,
        the primer, and a list of all tetramers, returns a pandas dataframe with
        all possible 2-mer flanks on each side of the core sequence.
    """
    designFlankDF = pd.read_csv(inputFileLoc)
    testDuplicates(designFlankDF, "Description", "createFlankProbes", "designFlankDF")
    designFlankDF["Sequence"] = designFlankDF["Sequence"].apply(lambda x: x.upper())
    # Generate Dataframes for each sequence, concatinate them
    flankDFList = []
    for desc, seq, coreSubStr in zip(designFlankDF["Description"],
                                     designFlankDF["Sequence"],
                                     designFlankDF["Core_Sequence"]):
        rStart = seq.find(coreSubStr)
        rEnd = rStart + len(coreSubStr)
        flankDFList.append(diFlankMutation(seq, [rStart, rEnd], desc, tetramerList))
    flankDF = pd.concat(flankDFList)
    # Add Primer to the sequences
    flankDF["Sequence"] = flankDF["Sequence"].apply(lambda x: x + primer)
    # Make 20 replicates for each sequence
    flankDF = addReps(flankDF, 20)
    return(flankDF)

def extractDeBruijn(inputFileLoc):
    """
    Given an array design file, extracts all de Bruijn sequences and returns a pandas
        dataframe with those sequences.
    """
    # Read in a previous array design with the de Bruijn sequences and same primer
    arrayDesignDF = pd.read_csv(inputFileLoc, sep = '\t')
    # Filter the dataframe for the de Bruijn sequences
    arrayDesignDF = arrayDesignDF.dropna(subset=["Sequence"])
    deBruDF = arrayDesignDF[arrayDesignDF["Name"].str.startswith("All_9mer_14bp") |
                           arrayDesignDF["Name"].str.startswith("All_8mer_13bp")].copy(deep = True)
    # Rename columns
    deBruDF = deBruDF[["Sequence", "Name"]]
    deBruDF = deBruDF.rename(columns={"Name":"Description"})
    # Test that the description has no duplicates, all sequences are 60bp
    testDuplicates(deBruDF, "Description", "extractDeBruijn", "deBruDF")
    deBruDF["Sequence"].apply(lambda x: testProbeLen(x, 60, "extractDeBruijn") )
    return(deBruDF) 

def createKDProbes(inputFileLoc, primer):
    """
    Given a file with kD probes, adds 20 replicates to each sequence and returns a 
        pandas dataframe
    """
    kDDF = pd.read_csv(inputFileLoc)
    # Check that there are no duplicates in the description and sequences are 25bp
    testDuplicates(kDDF, "Description", "createKDProbes", "kDDF")
    kDDF["Sequence"].apply(lambda x: testProbeLen(x, 25, "createKDProbes"))
    # Make a pandas dataframe with the reverse complement sequences
    kDDFRC = reverseComplementDF(kDDF)
    # Add primer to sequences
    kDDF["Sequence"] = kDDF["Sequence"].apply(lambda x: x + primer)
    kDDFRC["Sequence"] = kDDFRC["Sequence"].apply(lambda x: x + primer)
    # Concatinate both dataframes, prepend kD_ to the description, add replicates
    kDDF = pd.concat([kDDF, kDDFRC])
    kDDF["Description"] = kDDF["Description"].apply(lambda x: f"kD_{x}")
    kDDF = addReps(kDDF, 20)
    return(kDDF)

def createClustDiffProbes(proteinName, kmerFileLoc, clusterFileLoc, extension):
    """
    Given a protein name, kmer file of Escores for UV and Watson-Crick data, a cluster file,
        and an extension label, returns the top sequence per cluster with a description of the method for
        generating the sequence.
    """
    # Read in kmer file and cluster file
    kmerDF = pd.read_csv(kmerFileLoc, sep = '\t')
    clustDF = pd.read_csv(clusterFileLoc, sep = '\t', skiprows=2)
    # Calculate differences in UV vs Watson-Crick Escores
    kmerDF["Diff"] = kmerDF["Escore_UV9"] - kmerDF["Escore_WC8"]
    # Relate the difference to kmers as a dictionary
    kDiffDict = dict(zip(kmerDF["kmerFwd"], kmerDF["Diff"]))
    kDiffDict.update(dict(zip(kmerDF["kmerRC"], kmerDF["Diff"])))
    # Get the top sequence per cluster, using the kmer difference dictionary
    topSeqPerCluster = scoreClusters(clustDF, kDiffDict)
    topSeqPerCluster["Description"] = topSeqPerCluster["Cluster"].apply(lambda x: f"Clst_{proteinName}_C{x}_{extension}")
    return(topSeqPerCluster)

def createClustCloseProbes(proteinName, kmerFileLoc, clusterFileLoc, extension, seq, nProbes):
    """
    Given:
        1. protein name
        2. Location of the kmer file with UV and Watson-Crick measurements
        3. Location of the cluster file
        4. Extension, indicating the method used to generate the sequence
        5. Sequence to find the closest hamming distance of a substring in the clustered file.
        6. Number of probes to generate
    Returns a pandas dataframe with the sequences that have the smallest hamming distance the length
        of the query sequence
    """
    # Read in kmer file and cluster file
    kmerDF = pd.read_csv(kmerFileLoc, sep = '\t')
    clustDF = pd.read_csv(clusterFileLoc, sep = '\t', skiprows=2)
    # For each sequence, find the closest substring hamming distance
    hDistList = []
    for i in clustDF["Sequence"]:
        # Distance is the minimum of either orientation
        hDistFwd = closestSubHamming(i, seq)
        hDistRC = closestSubHamming(i, reverse_complement(seq))
        hDistList.append(min([hDistFwd, hDistRC]))
    clustDF["hDist"] = hDistList
    clustDF = clustDF.sort_values(by = "hDist").head(n = nProbes).reset_index(drop=True)
    desc = []
    for idx, hDist in enumerate(clustDF["hDist"]):
        desc.append(f"ClstD_{proteinName}_H{hDist}_{extension}_S{idx}")
    clustDF["Description"] = desc
    clustDF = clustDF[["Sequence", "Description"]]
    testDuplicates(clustDF, "Sequence", "createClustCloseProbes", "clustDF")
    testDuplicates(clustDF, "Description", "createClustCloseProbes", "clustDF")
    return(clustDF)    

def createClustDiffProbesFromLists(proteinNameList, 
                                  kmerFileLocList, 
                                  clusterKmersFileLocList,
                                  clusterShortExtFileLocList,
                                  primer,
                                  replicate_number):
    """
    Wrapper function for createClustDiffProbes to generate cluster probes for the array.
    Given a list of protein names, kmer file locations, and cluster file locations
     for both methods, a sequence to compare hamming distances, and the number of
     probes to generate for clost probes to that sequence. 
    Returns a pandas Dataframe with the probes. 
    """
    dfList = []
    # Process kmer clusters
    for protein, kmerFile, clusterKmer in zip(proteinNameList,
                                              kmerFileLocList,
                                              clusterKmersFileLocList):
        clusterDF = createClustDiffProbes(protein, kmerFile, clusterKmer, "TD25")
        dfList.append(clusterDF)
    # Process extended short clusters
    for protein, kmerFile, clusterKmer in zip(proteinNameList,
                                              kmerFileLocList,
                                              clusterShortExtFileLocList):
        clusterDF = createClustDiffProbes(protein, kmerFile, clusterKmer, "EXT")
        dfList.append(clusterDF)
    clusterProbeDF = pd.concat(dfList)
    clusterProbeDF = clusterProbeDF.reset_index(drop = True)
    clusterProbeRCDF = reverseComplementDF(clusterProbeDF)
    clusterProbeDF["Sequence"] = clusterProbeDF["Sequence"].apply(lambda x: x + primer)
    clusterProbeRCDF["Sequence"] = clusterProbeRCDF["Sequence"].apply(lambda x: x + primer)
    clusterProbeDF = addReps(clusterProbeDF, replicate_number)
    clusterProbeRCDF = addReps(clusterProbeRCDF, replicate_number)
    return(pd.concat([clusterProbeDF,clusterProbeRCDF]))

def createClustCloseProbesFromLists(proteinNameList, 
                                  kmerFileLocList, 
                                  clusterKmersFileLocList,
                                  clusterShortExtFileLocList,
                                  sequence,
                                  nProbes,
                                  primer,
                                  replicate_number):
    """
    Wrapper function for createClustCloseProbes. 
    """
    dfList = []
    # Process kmer clusters
    for protein, kmerFile, clusterKmer in zip(proteinNameList,
                                              kmerFileLocList,
                                              clusterKmersFileLocList):
        clusterDF = createClustCloseProbes(protein, kmerFile, clusterKmer, "TD25", sequence, nProbes)
        dfList.append(clusterDF)
    # Process extended short clusters
    for protein, kmerFile, clusterKmer in zip(proteinNameList,
                                              kmerFileLocList,
                                              clusterShortExtFileLocList):
        clusterDF = createClustCloseProbes(protein, kmerFile, clusterKmer, "EXT", sequence, nProbes)
        dfList.append(clusterDF)
    clusterProbeDF = pd.concat(dfList)
    clusterProbeDF = clusterProbeDF.reset_index(drop = True)
    clusterProbeRCDF = reverseComplementDF(clusterProbeDF)
    clusterProbeDF["Sequence"] = clusterProbeDF["Sequence"].apply(lambda x: x + primer)
    clusterProbeRCDF["Sequence"] = clusterProbeRCDF["Sequence"].apply(lambda x: x + primer)
    clusterProbeDF = addReps(clusterProbeDF, replicate_number)
    clusterProbeRCDF = addReps(clusterProbeRCDF, replicate_number)
    return(pd.concat([clusterProbeDF,clusterProbeRCDF]))   

def positionalProbes(desc, insert, baseSeq):
    """
    Given a description of the probe, an insert sequence to move across positions, and a base sequence
        that the insert sequence is inserted into, returns a pandas dataframe with the insert sequences at
        every second position and at the start of the probe that is 59 bases long as a special case. 
    """
    probeList, descList = [], []
    for i in range(int(((len(baseSeq) - len(insert) + 1)))):
        if i % 2 == 0:
            probe = baseSeq[:i] + insert + baseSeq[i + len(insert):]
            probeLen = len(probe)
            testProbeLen(probe, 25, "positionalProbes")
            probeList.append(probe)
    for i in range(len(probeList)):
        descList.append(f"Pos_{desc}_P{i * 2}")
    # Append shortened probe with 0 position to list
    shortProbe = insert[1:] + baseSeq[len(insert):]
    shortProbeLen = len(shortProbe)
    testProbeLen(shortProbe, 24, "positionalProbes")
    probeList.append(shortProbe)
    descList.append(f"Pos_{desc}_P0_Short")
    return(pd.DataFrame({"Sequence":probeList, "Description":descList}))

def createPositionalProbes(positionalProbeFile, primer):
    """
    Wrapper for positional probes
    """
    posProbeInputDF = pd.read_csv(positionalProbeFile)
    positionalProbeDFList = []
    for desc, insert, baseSeq, in zip(posProbeInputDF["Description"],
                                      posProbeInputDF["Sequence_Insert"],
                                      posProbeInputDF["Base_Sequence"]):
        positionalProbeDFList.append(positionalProbes(desc, insert, baseSeq))
    posProbeDF = pd.concat(positionalProbeDFList)
    posProbeDF = posProbeDF.reset_index(drop=True)
    posProbeDF = addReps(posProbeDF, 20)
    posProbeDF["Sequence"] = posProbeDF["Sequence"].apply(lambda x: x + primer)
    return(posProbeDF)

def saveUniqueProbes(probeDF,outputFileLoc = "../Array_Probes.txt"):
    """
    Makes a copy of the array design without replicates and saves it to a file.
    """
    # Generate a copy of the array without replicates
    uvArrayUnique = probeDF.copy(deep = True)
    uvArrayUnique["Description"] = uvArrayUnique["Description"].apply(lambda x: removeRepString(x))
    uvArrayUnique = uvArrayUnique.drop_duplicates()
    uvArrayUnique.to_csv(outputFileLoc, sep = '\t', header = None, index = False)

def saveAgilentFormattedArray(probeDFInput, outputFileLoc = "../Array_Design_Full.txt"):
    """
    Formats a dataframe of probes with a sequence and description in a format that is compatible with
        the format sent to Agilent to order custom microarrays. Saves this as a tab-seperated file.
    """
    probeDF = probeDFInput.copy(deep = True)
    # Generate UV Array file formated for ordering form Agilent
    IDList = []
    for idx, row in enumerate(probeDF['Sequence']):
        ID = f"Ctrl_UVT_{str(idx).zfill(6)}"
        if len(ID) > 15:
            raise ValueError(f"ProbeID {ID} is {len(ID)} characters and needs to be 15 or less")
        IDList.append(ID)
    probeDF['ID'] = IDList
    probeDF = probeDF[['ID', 'Sequence', 'Description']]
    probeDF['NA'] = 'NA|NA'
    probeDF['NA2'] = 'NA'
    probeDF['Desc2'] = probeDF['Description']
    probeDF['chr'] = 'chr1:0-0'
    probeDF.to_csv(outputFileLoc, sep = '\t', header = None, index = False)

def testProbeLenMatchDesc(probe, description):
    """
    Tests a probe to see if it is appropriate for the associated description.
    """
    if len(probe) != 60 and "Short" not in description:
        raise ValueError(f"Probe length should be 60, len = {len(i)}\n{i}")
    if len(probe) != 59 and "Short" in description:
        raise ValueError(f"Probe length should be 59, len = {len(i)}\n{i}")

                         
def testProbePrimer(probe, primer):
    """
    Tests a probe to see if it contains the primer sequence.
    """
    if primer not in probe:
        raise ValueError(f"Primer not found in sequence: {probe}")

    
def testDNDescription(probe, description):
    """
    For dinucleotide tests, parses the description for the probe to see if it
        matches the actual sequence. 
    """
    descList = description.split("_")
    posInfoList = descList[3].split("->")
    dnPos = posInfoList[0][3:]
    if dnPos.isnumeric():
        dnPos = int(dnPos)
        dnMut = posInfoList[1]
        if probe[dnPos:dnPos + 2] != dnMut:
            raise ValueError(f"In probe: {probe}\nthe mutation {dnMut} was not found at position {dnPos}")
    elif dnPos != "W":
        raise ValueError(f"In probe: {probe}\ndescription: {description}\n the postion could not be parsed.")

def testMultiDNDescription(probe, description):
    """
    Checks that there is a maximum of 2 description per probe. If so, 
        applies testDNDescription to each description for a probe.
    """
    if description[:2] == "DN":
        descriptionList = description.split("-D")
        nDesc = len(descriptionList)
        if  nDesc > 2:
            raise ValueError(f"There are {nDesc} description for DN probe {probe}")
        for desc in descriptionList:
            testDNDescription(probe, desc)

def testFlankingProbe(probe, description):
    """
    Given a probe for the all possible 2-mer flanks assay, tests that the description
        matches the probe.
    """
    if description[:4] == "Flank":
        descList = description.split("_")
        querySubstr = descList[3][:2] + descList[2] + descList[3][2:]
        if querySubstr not in probe:
            raise ValueError(f"In flanking probe: {probe}\nThe substring {querySubstr} was not found.")


def testDeBruijnProbes(arrayDF, arrayInfoList):
    """
    Checks that all UV de Bruijn probes are shared across array design files. If so, checks if the shared
        set is contained within the array design
    """
    # Check that all UV de Bruijn probes are shared across array design files
    arrayInfoDF = extractDeBruijn(f"../Sequence_Designs/{arrayInfoList[0]}")
    arrayLen = len(arrayInfoDF)
    for arrayInfoFile in arrayInfoList[1:]:
        newArray = extractDeBruijn(f"../Sequence_Designs/{arrayInfoFile}")
        arrayInfoDF = pd.merge(arrayInfoDF, newArray)
        if len(arrayInfoDF) != arrayLen:
            raise ValueError(f"Array design file {arrayInfoFile} has a different set of de Bruijn sequences")
    # Check that all are present in the array design
    if arrayLen != len(pd.merge(arrayDF, arrayInfoDF)):
        raise ValueError(f"UV array is missing {len(arrayInfoDF) - len(arrayDF)} de bruijn sequences")
        
def testPosProbes(df, posProbeFile):
    """
    Checks that the positional probe descriptions match the sequence. 
    """
    posDF = pd.read_csv(posProbeFile)
    posDict = dict(zip(posDF["Description"], posDF["Sequence_Insert"]))
    for TF in posDict.keys():
        qDF = df[df["Description"].str.startswith(f"Pos_{TF}")].copy(deep = True)
        qDF["Pos"] = qDF["Description"].apply(lambda x: int(x.split('_')[3][1:]))
        # Test that the insert sequence position in the description matches 
        for pos, seq in zip(qDF["Pos"], qDF["Sequence"]):
            if len(seq) == 60:
                qSeq = posDict[TF]
            elif len(seq) == 59:
                qSeq = posDict[TF][1:]
            else:
                raise ValueError(f"Sequence {seq} is not 59bp or 60bp long")
            if seq[pos:pos+len(qSeq)] != qSeq:
                raise ValueError(f"Sequence:{qSeq} for TF {TF}\nshould be at position {pos} in\n{seq}")
        # Test that positions given are in sets of 2
        uniquePos = list(set(qDF["Pos"]))
        for i in uniquePos:
            if i % 2 != 0:
                raise ValueError(f"An odd numbered position, {i}, was used in a positional assay design.")

def testClusterDescription(df):
    clustO1 = df[(df["Description"].str.startswith("Clst_")) &
                     (~df["Description"].str.contains("_O2_"))].copy(deep = True).reset_index(drop=True)
    # Test that all clusters are represented
    clustO1[["TF", "Cluster", "Assay"]] = clustO1["Description"].apply(lambda x: pd.Series(x.split("_")[1:4]))
    clustO1["Cluster"] = clustO1["Cluster"].apply(lambda x: int(x[1:]))
    assayDict = {"TD25":"25mers", "EXT":"extShort"}
    for nameTup, groupDF in clustO1.groupby(by=["TF", "Assay"]):
        cMin = min(groupDF["Cluster"])
        cMax = max(groupDF["Cluster"])
        cSet = set(groupDF["Cluster"])
        if cMin != 0:
            raise ValueError(f"Minimum cluster for {nameTup[0]} with method {nameTup[1]} should be 0, not {cMin}")
        for i in range(cMin, cMax):
            if i not in cSet:
                raise ValueError(f"Cluster {i} is missing from {nameTup[0]}, method {nameTup[1]}")
        # Test that sequence belongs to the cluster in the original file and is the max value
        clusterDF = pd.read_csv(f"../Data/Diff_UV_Probes/{nameTup[0]}_{assayDict[nameTup[1]]}_cluster.tsv",
                    sep = '\t', skiprows=2)
        testDuplicates(clusterDF, "Sequence", "Cluster", f"{nameTup}")
        clusterDict = dict(zip(clusterDF["Sequence"], clusterDF["Cluster"]))
        for seq, cluster in zip(groupDF["Sequence"], groupDF["Cluster"]):
            if cluster != clusterDict[seq[:25]]:
                raise ValueError((f"Sequence {seq} was described as being in cluster {cluster}," 
                                  "found to be in {clusterDict[seq[:25]]}"))
    # Test that for each cluster in O1 there is an O2
    clustO2 = df[(df["Description"].str.startswith("Clst_")) &
                     (df["Description"].str.contains("_O2_"))].copy(deep = True).reset_index(drop=True)
    clustO2[["TF", "Cluster", "Assay"]] = clustO2["Description"].apply(lambda x: pd.Series(x.split("_")[1:4]))
    clustO2["Cluster"] = clustO2["Cluster"].apply(lambda x: int(x[1:]))
    m = pd.merge(clustO1, clustO2, left_on = ["TF", "Cluster", "Assay"], right_on = ["TF", "Cluster", "Assay"])
    for s, rc in zip(m["Sequence_x"], m["Sequence_y"]):
        if s[:25] != reverse_complement(rc[:25]):
            raise ValueError(f"{s} was described as the reverse complement of {rc}")       

def testClusterDistanceDescription(df, testSeq):
    clust = df[(df["Description"].str.startswith("ClstD_"))].copy(deep = True).reset_index(drop=True)
    clust["SubHamming"] = clust["Description"].apply(lambda x: int(x.split('_')[2][1:]))
    for seq, hamDist in zip(clust["Sequence"], clust["SubHamming"]):
        distance1 = closestSubHamming(seq[:25], testSeq)
        distance2 = closestSubHamming(reverse_complement(seq[:25]), testSeq)
        if min([distance1, distance2]) != hamDist:
            raise ValueError(f"Minimum Sub-Hamming Distance does not match description for {seq}")

