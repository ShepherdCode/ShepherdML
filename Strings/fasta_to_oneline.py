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
        self.minlen = 0
        self.maxlen = None

    def set_min_len(self,minimum):
        self.minlen=minimum

    def set_max_len(self,maximum):
        self.maxlen=maximum

    def set_allcaps(self,setting):
        self.allcaps=setting

    def set_delete_N(self,setting):
        self.delete_N = setting

    def fix(self,verbose=True):
        build_seq=''
        defline=None
        seqs_in=0
        seqs_out=0
        with open(self.outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                for line in infa:
                    line=line.rstrip()
                    if line[0]==self.DEFLINE_PREFIX:
                        if seqs_in>0:
                            seqs_out += self.print_prev(outfa,defline,build_seq)
                        build_seq=''
                        defline=line
                        seqs_in += 1
                    else:
                        build_seq += line
                # Last sequence is special case
                seqs_out += self.print_prev(outfa,defline,build_seq)
        if verbose:
            print("Sequences_input: ",seqs_in)
            print("Sequences_output: ",seqs_out)

    def print_prev(self,outfile,defline,seq):
        NL='\n'
        allcaps=seq.upper()
        if self.delete_N and 'N' in allcaps:
            return 0
        elif self.minlen>0 and len(seq)<self.minlen:
            return 0
        elif self.maxlen is not None and len(seq)>self.maxlen:
            return 0
        outfile.write(defline+NL)
        if self.allcaps:
            outfile.write(allcaps)
        else:
            outfile.write(seq)
        outfile.write(NL)
        return 1

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--minlen', help='minimum sequence length (none)', type=int)
    parser.add_argument(
        '--maxlen', help='maximum sequence length (none)', type=int)
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
        if (args.minlen is not None):
            fixer.set_min_len(args.minlen)
        if (args.maxlen is not None):
            fixer.set_max_len(args.maxlen)
        if (args.ignore_case is not None):
            fixer.set_allcaps(args.ignore_case)
        if (args.delete_N is not None):
            fixer.set_delete_N(args.delete_N)
        fixer.fix()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
