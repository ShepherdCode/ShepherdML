import argparse
import traceback
import os
import sys
#import re
from datetime import datetime

DEFLINE_PREFIX='>'
NEWLINE='\n'

class Sequence():
    def __init__(self):
        self.defline=''
        self.seqline=''

class Filter():
    def processing_callback(self,sequence):
        return sequence # override this
    def process_one_sequence(self,sequence):
        sequence=self.processing_callback(sequence)
        return sequence

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
                            self.write_one_seq(lines,outfa)
                            lines=[]
                    lines.append(line)
                if len(lines)>0:
                    self.write_one_seq(lines,outfa)

    def write_one_seq(self,line_array,outfile):
        seq=Sequence()
        seq.defline=line_array[0]
        for line in line_array[1:]:
            seq.seqline += line
        for filter in self.filters:
            seq=filter.process_one_sequence(seq)
        if len(seq.defline)>0:
            outfile.write(seq.defline+NEWLINE)
            if len(seq.seqline)>0:
                outfile.write(seq.seqline+NEWLINE)

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
