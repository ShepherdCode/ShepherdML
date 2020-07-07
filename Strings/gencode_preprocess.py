import argparse
import traceback
import os
import sys
#import re
from datetime import datetime

class Streamer():
    def __init__(self):
        self.streamers=[]
    def add_streamer(self,streamer):
        self.streamers.append(streamer)
    def input(self,instr):
        self.lines = instr
    def output(self):
        lines=self.lines
        for streamer in self.streamers:
            streamer.input(lines)
            lines=streamer.output()
        return lines

class Gencode_Preprocess(Streamer):
    def __init__(self,debug=False):
        self.DEFLINE_PREFIX='>'
        self.debug=debug
        Streamer.__init__(self)

    def process_fasta(self,infile,outfile):
        lines=[]
        with open(outfile, 'w') as outfa:
            with open(infile, 'r') as infa:
                for line in infa:
                    line=line.rstrip()
                    if line[0]==self.DEFLINE_PREFIX:
                        if len(lines)>0:
                            self.process_seq(lines,outfa)
                            lines=[]
                    lines.append(line)
                if len(lines)>0:
                    self.process_seq(lines,outfa)

    def process_seq(self,lines,outfile):
        NL='\n'
        self.input(lines)
        revised=self.output()
        for line in revised:
            outfile.write(line+NL)

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
