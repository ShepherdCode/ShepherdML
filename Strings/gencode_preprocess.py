import argparse
import traceback
import os
import sys
#import re
from datetime import datetime

DEFLINE_PREFIX='>'
NEWLINE='\n'

class Filter():
    def processing_callback(self):
        pass # override this
    def max_one_sequence(self):
        pass # TO DO : raise exception
    def process_one_sequence(self,lines):
        self.lines=lines
        self.max_one_sequence()
        self.processing_callback()
        return self.lines

class Remove_Newlines(Filter):
    def processing_callback(self):
        defline=''
        seqline=''
        revised=[]
        for line in self.lines:
            if line[0]==DEFLINE_PREFIX:
                pass

class Gencode_Preprocess():
    def __init__(self,debug=False):
        self.debug=debug
        self.filters=[Filter()]

    def process_fasta(self,infile,outfile):
        lines=[]
        with open(outfile, 'w') as outfa:
            with open(infile, 'r') as infa:
                for line in infa:
                    line=line.rstrip()
                    if line[0]==DEFLINE_PREFIX:
                        if len(lines)>0:
                            self.process_one_seq(lines,outfa)
                            lines=[]
                    lines.append(line)
                if len(lines)>0:
                    self.process_one_seq(lines,outfa)

    def process_one_seq(self,lines,outfile):
        for filter in self.filters:
            lines=filter.process_one_sequence(lines)
        for line in lines:
            outfile.write(line+NEWLINE)

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
    '''Preprocess FASTA from GenCode or Ensembl.'''
    try:
        args_parse()
        fixer = Gencode_Preprocess(args.debug)
        fixer.process_fasta(args.infile,args.outfile)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
