import unittest
import os
from fasta_cleaner import Fasta_Cleaner

class TestFastaCleaner(unittest.TestCase):

    def setUp(self):
        self.FILENAME="test.fa"
    def tearDown(self):
        try:
            os.remove(self.FILENAME)
        except Exception as e:
            pass # This can happen even if the file is removed.

    def test_one_good_sequence(self):
        with open(self.FILENAME,'w') as fa:
            fa.write(">one_good")
            fa.write("ACGTAACCGGTT")
        fc = Fasta_Cleaner(self.FILENAME)
        fc.fix_everything()
        #self.assertEqual(readback,description)


if __name__ == '__main__':
    print ("Running unittest for fasta_cleaner.")
    unittest.main()
