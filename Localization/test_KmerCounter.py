'''
Test the KmerCounter class.
To run the tests...

python -m unittest
'''
import unittest
from KmerCounter import KmerCounter
class Test4mers(unittest.TestCase):
    def setUp(self):
        self.counter = KmerCounter(4)
    def test_4A(self):
        seq='AAAA'
        counts = self.counter.seq_to_kmer_counts(seq)
        self.assertEqual(counts[0],1)
        self.assertEqual(counts.sum(),1)
    def test_8A(self):
        seq='AAAAAAAA'
        counts = self.counter.seq_to_kmer_counts(seq)
        self.assertEqual(counts[0],5)
        self.assertEqual(counts.sum(),5)
    def test_middleN(self):
        seq='AAAANAAAA'
        counts = self.counter.seq_to_kmer_counts(seq)
        self.assertEqual(counts[0],2)
        self.assertEqual(counts.sum(),2)
    def test_middleN(self):
        seq='ACGTACGT'
        counts = self.counter.seq_to_kmer_counts(seq)
        self.assertEqual(counts[0],0)      # no AAAA
        self.assertEqual(counts.max(),2)   # ACTG twice     
        self.assertEqual(counts.sum(),5)   # all kmers counted
    def test_hash_homopolymers(self):
        c = self.counter
        self.assertEqual(c.hash_value('AAAA'),int("00000000",2))
        self.assertEqual(c.hash_value('CCCC'),int("01010101",2))
        self.assertEqual(c.hash_value('GGGG'),int("10101010",2))
        self.assertEqual(c.hash_value('TTTT'),int("11111111",2))  
    def test_hash_combos(self):
        c = self.counter
        self.assertEqual(c.hash_value('ACGT'),int("00"+"01"+"10"+"11",2))
        self.assertEqual(c.hash_value('TCGA'),int("11"+"01"+"10"+"00",2))
        self.assertEqual(c.hash_value('TATG'),int("11"+"00"+"11"+"10",2))
        self.assertEqual(c.hash_value('CTGG'),int("01"+"11"+"10"+"10",2))