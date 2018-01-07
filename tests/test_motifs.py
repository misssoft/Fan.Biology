'''
#Test mostprobable.py
#python -m unittest tests.test_motifs
'''
import unittest
from src import motifs

class TestMotifs(unittest.TestCase):

    def test_Count(self):
        input = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        expectedresult = {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}
        result = motifs.Count(input)
        self.assertEqual(expectedresult, result)

    def test_Profile_sample_dataset(self):
        input = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        expectedresult = {'A': [0.2, 0.4, 0.2, 0, 0, 0.4], 'C': [0.4, 0.2, 0.8, 0.4, 0, 0], 'G': [0.2, 0.2, 0, 0.4, 0.2, 0.2], 'T': [0.2, 0.2, 0, 0.2, 0.8, 0.4]}
        result = motifs.Profile(input)
        self.assertEqual(expectedresult, result)

    def test_Consensus(self):
        input = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        expectedresult = 'CACCTA'
        result = motifs.Consensus(input)
        self.assertEqual(expectedresult,result)

    def test_Score(self):
        input = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        expectedresult = 14
        result = motifs.Score(input)
        self.assertEqual(expectedresult,result)

    def test_Pr_1(self):
        input = 'GAGCTA'
        profile = {'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
                   'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
                   'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
                   'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}
        expectedresult = 0.0054
        result = motifs.Pr(input, profile)
        self.assertEqual(expectedresult,result)

    def test_Pr_2(self):
        input = 'ACGGGGATTACC'
        profile = {'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
                   'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
                   'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
                   'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]}
        expectedresult = 0.0008398080000000002
        result = motifs.Pr(input, profile)
        self.assertEqual(expectedresult,result)

    def test_ProfileMostProbablePattern(self):
        testText = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
        k = 5
        A = [0.2, 0.2, 0.3, 0.2, 0.3]
        C = [0.4, 0.3, 0.1, 0.5, 0.1]
        G = [0.3, 0.3, 0.5, 0.2, 0.4]
        T = [0.1, 0.2, 0.1, 0.1, 0.2]
        testProfile = {'A': A, 'C': C, 'G': G, 'T': T}

        expectedresult = 'CCGAG'
        result = motifs.ProfileMostProbablePattern(testText, 5, testProfile)
        self.assertEqual(expectedresult, result)

    def test_GreedyMotifSearch_sample_dataset(self):
        DnaTest = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        expectedresult = ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']

        result = motifs.GreedyMotifSearch(DnaTest, 3, 5)
        self.assertEqual(expectedresult, result)

    def test_GreedyMotifSearch_full_dataset(self):
        DnaTest = ['GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT', 'CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT',
                   'GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT', 'AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC',
                   'ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG']
        expectedresult = ['GTGCGT', 'GTGCGT', 'GCGCCA', 'GTGCCA', 'GCGCCA']
        result = motifs.GreedyMotifSearch(DnaTest, 6, 5)
        self.assertEqual(expectedresult, result)
