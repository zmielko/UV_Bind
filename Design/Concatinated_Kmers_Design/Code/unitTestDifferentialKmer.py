#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import unittest
from differentialKmer_OLS import *



class TestPreferenceExtension(unittest.TestCase):
    def test_extRight_trimer(self):
        """Tests if extRight properly extends a 3mer"""
        test_preference_set = set(["TCA", "CAT", "CAC", "CAG"])
        test_output_list = []
        expected_output_as_set = set(["ATCAT", "ATCAC", "ATCAG"])
        extRight("ATC" , 3, test_preference_set, test_output_list, set(["ATC"]))
        self.assertEqual(set(test_output_list), expected_output_as_set)
    def test_extRight_only_1bp_extension(self):
        """Tests if extRight will extend only 1bp, not repeating a 6mer"""
        test_preference_set = set(["ACATGT", "TGTCAC", "ATGTCA", "CATGTC"])
        test_output_list = []
        expected_output_as_set = set(["TGTCACCATGTCA"])
        extRight("TGTCACCATGTC" , 6, test_preference_set, test_output_list, set(["TGTCACCATGTC","TGTCAC"]))
        self.assertEqual(set(test_output_list), expected_output_as_set)
    def test_extRight_hexamer(self):
        """Tests if extRight properly extends a 6mer"""
        test_preference_set = set(["TTTTTT", "ATTTTT", "TTTTTA", "TTTTAC", "TTTTTC"])
        test_output_list = []
        expected_output_as_set = set(["TTTTTTAC", "TTTTTTC"])
        extRight("TTTTTT" , 6, test_preference_set, test_output_list, set(["TTTTTT"]))
        self.assertEqual(set(test_output_list), expected_output_as_set)
    def test_extLeft_trimer(self):
        """Tests if extLeft will properly extend a trimer"""
        test_preference_set = set(["TCA", "CAT", "CAC", "CAG"])
        test_output_list = []
        expected_output_as_set = set(["TCATC"])
        extLeft("ATC" , 3, test_preference_set, test_output_list, set(["ATC"]))
        self.assertEqual(set(test_output_list), expected_output_as_set)
    def test_extLeft_hexamer(self):
        """Tests if extLeft will properly extend multiple hexamers"""
        test_preference_set = set(["CCCCCC", "ACCCCC", "AAAAAA", "CACACA",
                                   "ACACAC", "TCACAC", "GCACAC"])
        test_output_list = []
        expected_output_as_set = set(["ACACACAG", "TCACACAG", "GCACACAG"])
        extLeft("ACACAG" , 6, test_preference_set, test_output_list, set(["ACACAG"]))
        self.assertEqual(set(test_output_list), expected_output_as_set)
    def test_extendKmers(self):
        test_kmer_list = ["CACGTG", "TATCGC"]
        test_pref_set = set(["ACACGT", "CCACGT", "GCCACG", "ACGTGA", "ATATCG"])
        test_output = extendKmers(test_kmer_list, test_pref_set)
        self.assertEqual(set(test_output), set(["ACACGTGA", "GCCACGTGA", "ATATCGC", "GCCACGTG"]))
        
class TestFilterOneOrientation(unittest.TestCase):
    def test_one_orientation(self):
        test_list = ["CAGCATCGTAGT", "CCCCC", "GGGGG", "ATA", "TAT"]
        test_output = filter_one_orientation(test_list)
        expected_output = ["CAGCATCGTAGT", "CCCCC", "ATA"]
        self.assertEqual(test_output, expected_output)
    def test_filter_orientation(self):
        input_list = ["CCCAAA", "TTTGGG", "CACACA", "CACGTG", "TGTGTG"]
        expected_output = ["CCCAAA", "CACACA", "CACGTG"]
        output = filter_one_orientation(input_list)
        self.assertEqual(output, expected_output)
        
class TestGreedyExtension(unittest.TestCase):
    def test_right_extension(self):
        test_seq = "ACATGCAT"
        test_kSet = set(["TGCATA", "TGCATC", "TGCATG", "CACATG", "ATGCATA"])
        test_dict = {"TGCATA": .2, "TGCATC": .3, "GCATCA": .3, "CACATG": .1, "ATGCATA": 10, "TGCATG":.3}
        test_output = extCaniR(test_seq, test_kSet,test_dict)
        self.assertEqual(test_output, "ACATGCATC")
    def test_left_extension(self):
        test_seq = "CTTTATC"
        test_kSet = set(["ACTTTA", "CCTTTA", "GCTTTA", "TCTTTA", "ACTTT"])
        test_dict = {"ACTTTA":.4, "CCTTTA":.1, "GCTTTA":-5, "TCTTTA":0, "TCTTT":10}
        test_output = extCaniL(test_seq, test_kSet,test_dict)
        self.assertEqual(test_output, "ACTTTATC")
    def test_extend25(self):
        test_seq = "CTAGCTAGCTTAGCTAGCTAGCT"
        test_kSet = set(["TAGCTA", "TAGCTC", "AGCTCG"])
        test_dict = {"TAGCTA":.4, "TAGCTC":.5, "AGCTCG":.4}
        test_output = extend25(test_seq, test_kSet,test_dict)  
        self.assertEqual(test_output, "CTAGCTAGCTTAGCTAGCTAGCTCG")
    def test_extend_short_probes(self):
        test_probe_list = ["C", "CTAGCT", "CTAGCTAGCTTAGCTAGCTAGCT", "CTAGCTAGCTTAGCTAGCTAGCTCGAA"]
        test_kSet = set(["TAGCTA", "TAGCTC", "AGCTCG", "TCGAAT","AC"])
        test_dict = {"TAGCTA":.4, "TAGCTC":.5, "AGCTCG":.4, "TCGAAT":.9 ,"AC":100}
        test_output = extendShortProbes(test_probe_list, test_kSet,test_dict)
        self.assertEqual(test_output, ["CTAGCTAGCTTAGCTAGCTAGCTCG"])

class Test25merGeneration(unittest.TestCase):
    def test_window_search_kmer(self):
        test_probe = "TTTTTTTTAGC"
        expected_output = set(["TTTTTT", "TTTTTA", "TTTTAG", "TTTAGC"])
        self.assertEqual(windowSearchKmer(test_probe, 6), expected_output)
    def test_all_25mers_long_probes(self):
        test_probe_list = ["AAA", "CACGTG", "CTAGCTAGCTAGCTAGCTTAGCGAT",
                           "ATCGATATAGCTATAGATCGATCGTTT", "ACTAGCTAGCTAGCTAGCTTAGCGAT"]
        test_output = all25mersLongProbes(test_probe_list)
        expected_output = set(["CTAGCTAGCTAGCTAGCTTAGCGAT",
                          "ATCGATATAGCTATAGATCGATCGT",
                          "TCGATATAGCTATAGATCGATCGTT",
                          "CGATATAGCTATAGATCGATCGTTT",
                          "ACTAGCTAGCTAGCTAGCTTAGCGA"])
        self.assertEqual(set(test_output), expected_output)

