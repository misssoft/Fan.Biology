'''
#Test randomsearch.py
#python -m unittest tests.test_random_search (under Fan.Biology)
'''
import unittest
from src import random_search


class TestMostProbable(unittest.TestCase):
    def test_sample_dataset(self):
        input = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
        expectedoutput = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        result = random_search.Normalize(input)
        self.assertEqual(expectedoutput, result)

    def test_weighedDie(self):
        dic = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
        result = random_search.WeightedDie(dic)
        self.assertTrue(result=='A' or result=='C' or  result=='G' or result=='T')

    def test_ProfileGeneratedString(self):
        dna = 'AAACCCAAACCC'
        profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
        k = 2

        result = random_search.ProfileGeneratedString(dna,profile,k)
        self.assertTrue(len(result) == 2)