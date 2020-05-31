import unittest
from KmerCount import KmerHistogram
from KmerCount import KmerHistNaive
from KmerCount import KmerHistHarvester

class TestKmerHistNaive(unittest.TestCase):
    def test_factory_k2(self):
        k=2
        T="abaababaabaab"
        factory=KmerHistNaive(k)
        factory.add_string(T)
        hist=factory.get_histogram()
        self.assertEqual(k,hist.get_k())
        freq3=hist.get_freq("aa")
        self.assertEqual(3,freq3)
        freq4=hist.get_freq("ba")
        self.assertEqual(4,freq4)
        freq5=hist.get_freq("ab")
        self.assertEqual(5,freq5)
    def test_factory_k3(self):
        k=3
        T="abaababaabaab"
        factory=KmerHistNaive(k)
        factory.add_string(T)
        hist=factory.get_histogram()
        self.assertEqual(k,hist.get_k())
        freq4=hist.get_freq("aba")
        self.assertEqual(4,freq4)

class TestKmerHistogram(unittest.TestCase):
    def setUp(self):
        self.k = 5
        self.hist = KmerHistogram(self.k)
        self.valid_kmer="abcde"
        self.unused_kmer=self.valid_kmer.upper()
        self.invalid_kmer="abc"
        self.count = 20
        for i in range (0,self.count):
            self.hist.increment_count(self.valid_kmer)
    def test_histogram_k(self):
        self.assertEqual(self.k,self.hist.get_k())
    def test_histogram_positive_freq(self):
        self.assertEqual(self.count,self.hist.get_freq(self.valid_kmer))
    def test_histogram_zero_freq(self):
        self.assertEqual(0,self.hist.get_freq(self.unused_kmer))
    def test_histogram_wrong_size(self):
        self.assertRaises(Exception,self.hist.get_freq,self.invalid_kmer)
        self.assertRaises(Exception,self.hist.increment_count,self.invalid_kmer)

class TestKmerHistHarvester(unittest.TestCase):
    def test_smaller_k(self):
        k=3
        factory=KmerHistNaive(k)
        k3=factory.get_histogram().get_k()
        factory=None
        factory=KmerHistHarvester(2)
        k2=factory.get_histogram().get_k()
        self.assertEqual(k3-1, k2)
    def test_combined_freq(self):
        T="abaababaabaab"
        factory=KmerHistNaive(2)
        factory.add_string(T)
        hist2=factory.get_histogram()
        factory=None
        factory=KmerHistNaive(3)
        factory.add_string(T)
        hist3=factory.get_histogram()
        factory=None
        factory=KmerHistHarvester(2)
        factory.extract_from(hist3)
        histE=factory.get_histogram()
        factory=None
        self.assertEqual(hist2.get_items(),histE.get_items())

if __name__ == '__main__':
    '''python3 -m unittest -v TestKmer.py'''
    unittest.main()
