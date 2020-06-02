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

    def __init__(self,fastafile):
        self.fastafile = fastafile

    def fix_everything(self):
        with open(self.fastafile, 'r') as fa:
            print("Fixed!")

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'fastafile', help='input filename', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''Inspect & fix a FASTA file of DNA sequences.'''
    try:
        args_parse()
        fixer = Fasta_Cleaner(args.fastafile)
        fixer.fix_everything()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
