import argparse
import traceback
import os
import sys
from datetime import datetime

def say(say1="",say2=""):
    if args.debug:
        print (datetime.now(),end=' ')
        print (say1,end=' ')
        print (say2)

class Fasta_Cleaner():
    def __init__(self,infile,outfile):
        self.infile = infile
        self.outfile = outfile
        self.DEFLINE_PREFIX='>'
        self.NEWLINE='\n'

    def fix_everything(self):
        NL=self.NEWLINE
        with open(self.outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                for line in infa:
                    if line[0]==self.DEFLINE_PREFIX:
                        defline=line
                        outfa.write(defline)
                    else:
                        partnum = 1
                        allcaps = line.upper()
                        nextpos=allcaps.find('N')
                        if nextpos < 0:
                            outfa.write(allcaps)
                        else:
                            allcaps=allcaps.rstrip()
                            defline=defline.rstrip()
                            prefix=allcaps[0:nextpos]
                            outfa.write(prefix+NL)
                            suffix=allcaps[nextpos+1:]
                            partnum += 1
                            newdefline=defline+"-part-"+str(partnum)
                            outfa.write(newdefline+NL)
                            outfa.write(suffix+NL)

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
        fixer = Fasta_Cleaner(args.infile,args.outfile)
        fixer.fix_everything()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
