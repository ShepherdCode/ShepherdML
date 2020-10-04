import argparse
import traceback
import os
import sys
import random

class Reducer ():
    def __init__(self,debug=False):
        self.DEFLINE_PREFIX='>'
        self.debug=debug
        self.test_size = 500
        self.challenge_size = 500
        self.valid_size = 500
        self.train_size = 500
        self.all_genes=[]

    def data_in(self,infile):
        self.all_genes=[]
        with open(infile, 'r') as infa:
            num_seqs = 0
            for line in infa:
                if line[0]==self.DEFLINE_PREFIX:
                    transcript_id=line[1:18]
                    gene_id=line[19:36]
                    self.all_genes.append(gene_id)
                else:
                    seq=line
        print(len(self.all_genes))
        print(self.all_genes[0])
        print(self.all_genes[1])
        print(self.all_genes[len(self.all_genes)-1])

    def reduce():
        pass

    def data_out(self,infile):
        #with open(self.outfile, 'w') as outfa:
        pass

    def reduce(self):
        pass

    def print_prev(self,outfile,defline,seq):
        NL='\n'
        allcaps=seq.upper()
        if self.delete_N and 'N' in allcaps:
            pass
        else:
            outfile.write(defline+NL)
            if self.allcaps:
                outfile.write(allcaps)
            else:
                outfile.write(seq)
            outfile.write(NL)

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''FASTA reduction ala mRNN.'''
    try:
        args_parse()
        fixer = Reducer(args.debug)
        fixer.data_in(args.infile)
        fixer.reduce()
        fixer.data_out(args.infile)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
