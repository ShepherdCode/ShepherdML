import argparse
import traceback
import os
import sys

class Fasta_Oneliner():
    def __init__(self,infile,outfile,debug=False):
        self.infile = infile
        self.outfile = outfile
        self.DEFLINE_PREFIX='>'
        self.debug=debug
        self.allcaps = False
        self.delete_N = False

    def set_allcaps(self,setting):
        self.allcaps=setting

    def set_delete_N(self,setting):
        self.delete_N = setting

    def fix(self):
        build_seq=''
        defline=None
        with open(self.outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                num_seqs = 0
                for line in infa:
                    line=line.rstrip()
                    if line[0]==self.DEFLINE_PREFIX:
                        if num_seqs>0:
                            self.print_prev(outfa,defline,build_seq)
                        build_seq=''
                        defline=line
                        num_seqs += 1
                    else:
                        build_seq += line
                # Last sequence is special case
                self.print_prev(outfa,defline,build_seq)
        if (self.debug):
            print("Processed %d input sequences."%num_seqs)

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
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--delete_N', help='Delete all sequences containing N.',
        action='store_true')
    parser.add_argument(
        '--ignore_case', help='Output DNA all caps.',
        action='store_true')
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''FASTA input multiliners output oneliners.'''
    try:
        args_parse()
        fixer = Fasta_Oneliner(args.infile,args.outfile,args.debug)
        fixer.set_allcaps(args.ignore_case)
        fixer.set_delete_N(args.delete_N)
        fixer.fix()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
