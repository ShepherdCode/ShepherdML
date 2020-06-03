import unittest
import os
import shutil #shutil.copy2(self.FILE1,self.FILE3)
from fasta_cleaner import Fasta_Cleaner

def writeline(fp,text):
    fp.write(text)
    fp.write('\n')

class TestFastaCleaner(unittest.TestCase):
    def setUp(self):
        self.FILE1="test.in.fa"
        self.FILE2="test.out.fa"
        self.FILE3="test.answer.fa"
        self.cleanup=True
    def tearDown(self):
        if (self.cleanup):
            try:
                os.remove(self.FILE1)
                os.remove(self.FILE2)
                os.remove(self.FILE3)
            except Exception as e:
                pass # This can happen even if the file is removed.
        else:
            print("Remember to remove the test files.")

    def test_valid_letter(self):
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        self.assertTrue(fc.is_valid('A'))
        self.assertTrue(fc.is_valid('C'))
        self.assertTrue(fc.is_valid('G'))
        self.assertTrue(fc.is_valid('T'))
        self.assertFalse(fc.is_valid('t'))
        self.assertFalse(fc.is_valid('N'))
        self.assertFalse(fc.is_valid('b'))
        self.assertFalse(fc.is_valid('Y'))

    def test_next_valid_position(self):
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        self.assertEqual(fc.next_valid_position('ACGNBG',0),0)
        self.assertEqual(fc.next_valid_position('NBYACG',0),3)
        self.assertEqual(fc.next_valid_position('NBYYYG',0),5)
        self.assertEqual(fc.next_valid_position('NBYnby',0),-1)
        self.assertEqual(fc.next_valid_position('NBYACG',3),3)
        self.assertEqual(fc.next_valid_position('NBYACG',2),3)
        self.assertEqual(fc.next_invalid_position('NBYYYG',0),0)
        self.assertEqual(fc.next_invalid_position('ACGGTT',0),-1)

    def test_one_good_sequence(self):
        with open(self.FILE1,'w') as fa:
            writeline(fa,">one_good")
            writeline(fa,"ACGTAACCGGTT")
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE1,self.FILE2))
        self.assertEqual(same,0)

    def test_invalid_letters(self):
        good_sequence="ACGTAA"
        bad_sequence="BDHKMNRSVWYbdhkmnrsvwy"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">bad_letters")
            writeline(fa,good_sequence+bad_sequence)
            writeline(fa,bad_sequence)
        with open(self.FILE3,'w') as fa:
            writeline(fa,">bad_letters")
            writeline(fa,good_sequence)
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)

    def test_upper_and_lower_case(self):
        sequence="ACGTAaCcGgTt"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">upper_and_lower")
            writeline(fa,sequence)
            writeline(fa,sequence)
        with open(self.FILE3,'w') as fa:
            writeline(fa,">upper_and_lower")
            concat = sequence+sequence
            writeline(fa,concat.upper())
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)

    def test_one_N(self):
        sequence1="ACGTA"
        sequence2="N"
        sequence3="TTCA"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">one_N")
            writeline(fa,sequence1+sequence2+sequence3)
        with open(self.FILE3,'w') as fa:
            writeline(fa,">one_N")
            writeline(fa,sequence1)
            writeline(fa,">one_N-part-2")
            writeline(fa,sequence3)
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)

    def test_multiple_N(self):
        sequence1="ACGTA"
        nsequence="NNNNN"
        sequence2="TTCAA"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">multiple_N")
            writeline(fa,sequence1+nsequence+sequence2+nsequence+sequence1)
        with open(self.FILE3,'w') as fa:
            writeline(fa,">multiple_N")
            writeline(fa,sequence1)
            writeline(fa,">multiple_N-part-2")
            writeline(fa,sequence2)
            writeline(fa,">multiple_N-part-3")
            writeline(fa,sequence1)
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)

    def test_N_at_start_and_end(self):
        sequence1="ACGTA"
        nsequence="NNNNN"
        sequence2="TTCAA"
        with open(self.FILE1,'w') as fa:
            writeline(fa,">start_end_N")
            writeline(fa,nsequence+sequence1+nsequence)
            writeline(fa,nsequence+sequence2+nsequence)
        with open(self.FILE3,'w') as fa:
            writeline(fa,">start_end_N")
            writeline(fa,sequence1)
            writeline(fa,">start_end_N-part-2")
            writeline(fa,sequence2)
        fc = Fasta_Cleaner(self.FILE1,self.FILE2)
        fc.fix_everything()
        same=os.system("diff %s %s"%(self.FILE2,self.FILE3))
        self.assertEqual(same,0)

if __name__ == '__main__':
    print ("Running unittest for fasta_cleaner.")
    unittest.main()
