import argparse
import traceback
import os
import sys
from datetime import datetime
from math import floor
import random

DEFLINE_PREFIX='>'
NEWLINE='\n'
FIELD_SEPARATOR='|'

class Parser():
    def parse_defline(defline):
        try:
            fields=defline.split(FIELD_SEPARATOR)
            transcript=fields[0]
            gene=fields[1]
            length=int(fields[6])
        except Exception as e:
            print('Cannot parse defline '+defline)
            print(e)
            raise Exception
        return (transcript,gene,length)

class SetAside():
    def __init__(self,infile,debug=False):
        self.debug=debug
        self.infile=infile
        self.test_genes={}

    def process_fasta(self,outfile,include,verbose=True):
        '''Assume FASTA with one line per sequence.'''
        defline=None
        seqline=None
        sn=0  # output sequence number
        with open(outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                for line in infa:
                    if line[0]==DEFLINE_PREFIX:
                        defline = line
                    else:
                        seqline = line
                        (tid,gid,slen)=Parser.parse_defline(defline)
                        if include and gid in self.test_genes or \
                        not include and gid not in self.test_genes:
                            outfa.write(defline)
                            outfa.write(seqline)
                            sn += 1
        if verbose:
            print("Output_file: ",outfile)
            print("Output_sequences: ",sn)

    def choose_test_set(self,tsize,verbose=True):
        all_genes={}
        self.test_genes={}
        with open(self.infile, 'r') as infa:
            for line in infa:
                if line[0]==DEFLINE_PREFIX:
                    (tid,gid,slen)=Parser.parse_defline(line)
                    all_genes[gid] = 1  # signal existence
        shuffled = list(all_genes.keys())
        random.seed(45)
        random.shuffle(shuffled)
        for gid in shuffled[:tsize]:
            self.test_genes[gid] = 1 # signal existence
        if verbose:
            print("Genes_to_test_set: ",len(self.test_genes.keys()))

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Preprocess GenCode FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'testfile', help='test filename (fasta)', type=str)
    parser.add_argument(
        'trainfile', help='train filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        fixer = SetAside(args.infile,args.debug)
        fixer.choose_test_set(1000)
        fixer.process_fasta(args.testfile,True)
        fixer.process_fasta(args.trainfile,False)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
