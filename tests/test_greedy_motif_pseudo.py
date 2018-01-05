'''
#Test mostprobable.py
#python -m unittest tests.test_greedy_motif_pseudo (under Fan.Biology)
'''
import unittest
from src import greedy_motif_pseudo


class TestMostProbable(unittest.TestCase):

    def test_CountWithPseudocounts_sample_dataset(self):
        motif = [
            'AACGTA',
            'CCCGTT',
            'CACCTT',
            'GGATTA',
            'TTCCGG'
        ]
        output = {'A': [2, 3, 2, 1, 1, 3], 'C': [3, 2, 5, 3, 1, 1], 'T': [2, 2, 1, 2, 5, 3], 'G': [2, 2, 1, 3, 2, 2]}
        result = greedy_motif_pseudo.CountWithPseudocounts(motif)
        self.assertEqual(result, output)

    def test_CountWithPseudocounts_full_dataset(self):
        motif = [
            'GTACAACTGT',
            'CAACTATGAA',
            'TCCTACAGGA',
            'AAGCAAGGGT',
            'GCGTACGACC',
            'TCGTCAGCGT',
            'AACAAGGTCA',
            'CTCAGGCGTC',
            'GGATCCAGGT',
            'GGCAAGTACC'
            ]
        output = {'A': [3, 4, 4, 4, 7, 5, 3, 3, 2, 4], 'C': [3, 4, 5, 4, 3, 4, 3, 2, 4, 4], 'T': [3, 3, 1, 5, 2, 1, 3, 3, 2, 5], 'G': [5, 3, 4, 1, 2, 4, 5, 6, 6, 1]}
        result = greedy_motif_pseudo.CountWithPseudocounts(motif)
        self.assertEqual(result, output)

    def test_GreedyMotifSearchWithPseudocounts_sample_dataset(self):
        motif = [
            'GGCGTTCAGGCA',
            'AAGAATCAGTCA',
            'CAAGGAGTTCGC',
            'CACGTCAATCAC',
            'CAATAATATTCG'
        ]
        k = 3
        t = 5
        output = [
            'TTC',
            'ATC',
            'TTC',
            'ATC',
            'TTC'
        ]
        result = greedy_motif_pseudo.GreedyMotifSearchWithPseudocounts(motif,k,t)
        self.assertEqual(result, output)

    def test_GreedyMotifSearchWithPseudocounts_first_kmer(self):
        motif = [
                'AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC',
                'ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC',
                'AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT',
                'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
                'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
                'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
                'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
                'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA'
            ]
        k = 5
        t = 8
        output = [
                'AGGCG',
                'ATCCG',
                'AAGCG',
                'AGTCG',
                'AACCG',
                'AGGCG',
                'AGGCG',
                'AGGCG'
            ]
        result = greedy_motif_pseudo.GreedyMotifSearchWithPseudocounts(motif, k, t)
        self.assertEqual(result, output)

    def test_GreedyMotifSearchWithPseudocounts_last_kmer(self):
        motif = [
            'GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG',
            'TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA',
            'AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC',
            'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
            'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
            'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
            'AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG',
            'AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG'
            ]
        k = 5
        t = 8
        output = [
            'AGGCG',
            'TGGCA',
            'AAGCG',
            'AGGCA',
            'CGGCA',
            'AGGCG',
            'AGGCG',
            'AGGCG'
            ]
        result = greedy_motif_pseudo.GreedyMotifSearchWithPseudocounts_Workaround(motif, k, t)
        self.assertEqual(result, output)

    def test_GreedyMotifSearchWithPseudocounts_breaking_ties(self):
        motif = [
            'GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC',
            'TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC',
            'TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG',
            'GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG',
            'GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG',
            'TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG',
            'GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG',
            'AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG',
        ]
        k = 5
        t = 8
        output = [
            'GGCGG',
            'GGCTC',
            'GGCGG',
            'GGCAG',
            'GACGG',
            'GACGG',
            'GGCGC',
            'GGCGC'
        ]
        result = greedy_motif_pseudo.GreedyMotifSearchWithPseudocounts(motif, k, t)
        self.assertEqual(result, output)

    def test_GreedyMotifSearchWithPseudocounts_full_dataset(self):
        motif = [
            'GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC',
            'TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC',
            'TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG',
            'GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG',
            'GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG',
            'TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG',
            'GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG',
            'AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG',
        ]
        k = 3
        t = 8
        output = [
            'GGC',
            'GGC',
            'GGC',
            'GGC',
            'GGC',
            'GGC',
            'GGC',
            'GGC'
        ]
        result = greedy_motif_pseudo.GreedyMotifSearchWithPseudocounts(motif, k, t)
        self.assertEqual(result, output)
