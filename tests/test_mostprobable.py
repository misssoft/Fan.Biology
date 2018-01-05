'''
#Test mostprobable.py
#python -m unittest tests.test_mostprobable (under Fan.Biology)
'''
import unittest
from src import mostprobable


class TestMostProbable(unittest.TestCase):

    def test_sample_dataset(self):
        dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        k = 5
        profile = {
            'A':[0.2, 0.2, 0.3, 0.2, 0.3],
            'C':[0.4, 0.3, 0.1, 0.5, 0.1],
            'G':[0.3, 0.3, 0.5, 0.2, 0.4],
            'T':[0.1, 0.2, 0.1, 0.1, 0.2]
        }
        result = mostprobable.ProfileMostProbablePattern(dna,k, profile)
        self.assertEqual(result, 'CCGAG')

    def test_first_kmer(self):
        dna = 'AGCAGCTTTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATCTGAACTGGTTACCTGCCGTGAGTAAAT'
        k = 8
        profile = {
            'A':[0.7,0.2,0.1,0.5,0.4,0.3,0.2,0.1],
            'C':[0.2,0.2,0.5,0.4,0.2,0.3,0.1,0.6],
            'G':[0.1,0.3,0.2,0.1,0.2,0.1,0.4,0.2],
            'T':[0.0,0.3,0.2,0.0,0.2,0.3,0.3,0.1]
        }
        result = mostprobable.ProfileMostProbablePattern(dna,k, profile)
        self.assertEqual(result, 'AGCAGCTT')

    def test_last_kmer(self):
        dna = 'TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA'
        k = 12
        profile = {
            'A': [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.1,0.2,0.3,0.4,0.5],
            'C': [0.3,0.2,0.1,0.1,0.2,0.1,0.1,0.4,0.3,0.2,0.2,0.1],
            'G': [0.2,0.1,0.4,0.3,0.1,0.1,0.1,0.3,0.1,0.1,0.2,0.1],
            'T': [0.3,0.4,0.1,0.1,0.1,0.1,0.0,0.2,0.4,0.4,0.2,0.3]
        }
        result = mostprobable.ProfileMostProbablePattern(dna, k, profile)
        self.assertEqual(result, 'AAGCAGAGTTTA')

    def test_ties(self):
        dna = 'AACCGGTT'
        k = 3
        profile = {
            'A': [1.0,1.0,1.0],
            'C': [0.0,0.0,0.0],
            'G': [0.0,0.0,0.0],
            'T': [0.0,0.0,0.0]
        }
        result = mostprobable.ProfileMostProbablePattern(dna, k, profile)
        self.assertEqual(result, 'AAC')

    def test_full_dataset(self):
        dna = 'TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA'
        k = 5
        profile = {
            'A': [0.2, 0.2, 0.3, 0.2, 0.3],
            'C': [0.4, 0.3, 0.1, 0.5, 0.1],
            'G': [0.3, 0.3, 0.5, 0.2, 0.4],
            'T': [0.1, 0.2, 0.1, 0.1, 0.2]
        }
        result = mostprobable.ProfileMostProbablePattern(dna, k, profile)
        self.assertEqual(result, 'CAGCG')


if __name__ == '__main__':
    unittest.main()
