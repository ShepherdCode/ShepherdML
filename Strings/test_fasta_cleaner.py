import unittest
import os
from fasta_cleaner import Fasta_Cleaner

def writeline(fp,text):
    fp.write(text + '\n')

class TestFastaCleaner(unittest.TestCase):

    def setUp(self):
        self.FILE1="test.in.fa"
        self.FILE2="test.out.fa"
        self.FILE3="test.answer.fa"
    def tearDown(self):
        try:
            #os.remove(self.FILE1)
            #os.remove(self.FILE2)
            #os.remove(self.FILE3)
            pass
        except Exception as e:
            pass # This can happen even if the file is removed.

    def test_one_good_sequence(self):
        with open(self.FILE1,'w') as fa:
            writeline(fa,">one_good")
            writeline(fa,"ACGTAACCGGTT")
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE1,self.FILE2))
        self.assertEqual(same,0)

    def test_upper_and_lower_case(self):
        sequence="ACGTAaCcGgTt"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">upper_and_lower")
            writeline(fa,sequence)
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        with open(self.FILE3,'w') as fa:
            writeline(fa,">upper_and_lower")
            writeline(fa,sequence.upper())
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)


if __name__ == '__main__':
    print ("Running unittest for fasta_cleaner.")
    unittest.main()
