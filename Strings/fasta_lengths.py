import argparse
import traceback
import os
import sys

class Fasta_Lengths():
    def __init__(self,infile,outfile,debug=False):
        self.infile = infile
        self.DEFLINE_PREFIX='>'
        self.debug=debug
        self.noname=False

    def set_noname(self,setting):
        self.noname = setting

    def report(self):
        with open(self.infile, 'r') as infa:
            num_seqs = 0
            for line in infa:
                line=line.rstrip()
                if line[0]==self.DEFLINE_PREFIX:
                    defline=line[1:]
                    num_seqs += 1
                else:
                    slen = len(line)
                    if (self.noname):
                        print("%d"%slen)
                    else:
                        print("%d %s"%(slen,defline))

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        '--noname', help='Print just lengths.',
        action='store_true')
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''FASTA input multiliners output oneliners.'''
    try:
        args_parse()
        s = Fasta_Lengths(args.infile,args.debug)
        s.set_noname(args.noname)
        s.report()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
