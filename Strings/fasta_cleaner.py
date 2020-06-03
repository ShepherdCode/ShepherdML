import argparse
import traceback
import os
import sys
from datetime import datetime

class Fasta_Cleaner():
    def __init__(self,infile,outfile,debug=False):
        self.infile = infile
        self.outfile = outfile
        self.DEFLINE_PREFIX='>'
        self.debug=debug

    def fix_everything(self):
        Nchar='N'
        prev_seq=[]
        defline=None
        with open(self.outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                num_seqs = 0
                for line in infa:
                    line=line.rstrip()
                    if line[0]==self.DEFLINE_PREFIX:
                        self.print_prev(outfa,defline,prev_seq)
                        prev_seq=[]
                        defline=line
                        ends_in_n = False  # can we get rid of this?
                        num_seqs += 1
                    else:
                        if ends_in_n:
                            prev_seq.append(0)  # indicate Ns came before
                        ends_in_n = Nchar==line[-1]
                        allcaps = line.upper()
                        allcaps = self.remove_leading_nrun(allcaps)
                        (prefix,suffix)=self.split_first_nrun(allcaps)
                        if len(prefix)>0:
                            prev_seq.append(prefix)
                        while suffix is not None:
                            (prefix,suffix)=self.split_first_nrun(suffix)
                            if len(prefix)>0:
                                prev_seq.append(0)  # indicate Ns came before
                                prev_seq.append(prefix)
                # Last sequence is special case
                self.print_prev(outfa,defline,prev_seq)
        if (self.debug):
            print("Fixed %d sequences."%num_seqs)

    def print_prev(self,outfile,defline,seqs):
        NL='\n'
        part = 1
        if defline is not None and len(seqs)>0:
            outfile.write(defline+NL)
            for seq in seqs:
                if seq==0:  # indicator that there was N here
                    part += 1
                    continuation=defline+"-part-"+str(part)
                    outfile.write(NL+continuation+NL)
                else:
                    outfile.write(seq)
            outfile.write(NL)

    def remove_leading_nrun(self,str):
        Nchar='N'
        nextpos=str.find(Nchar)
        if nextpos==0:
            # yes, it starts with N
            while nextpos<len(str) and str[nextpos]==Nchar:
                # skip to first non-N
                nextpos += 1
            if nextpos>=len(str):
                # string was all N
                str=""
            else:
                # string started with N
                str=str[nextpos:]
        return str

    def split_first_nrun(self,str):
        Nchar='N'
        prefix=str
        suffix = None
        nextpos=str.find(Nchar)
        if nextpos==0:
            raise Exception ("Not expecting a string that starts with N.")
        if nextpos>0:
            prefix=str[:nextpos]
            while nextpos<len(str) and str[nextpos]==Nchar:
                nextpos += 1
            if nextpos < len(str):
                suffix = str[nextpos:]
        return (prefix,suffix)

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''Inspect & fix a FASTA file of DNA sequences.'''
    try:
        args_parse()
        fixer = Fasta_Cleaner(args.infile,args.outfile,args.debug)
        fixer.fix_everything()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
