import argparse
import traceback
import sys

class SortedSuffixes:
    def __init__(self):
        self.suffixes=[]
        self.T=""
    def set_string(self,T):
        self.T=""
        self.add_string(T)
    def add_string(self,TT):
        # each suffix goes up to the first $ only
        pos=len(self.T) # start counter at end of previous T
        TTT=TT.rstrip("\n\r")+"$" # replace newline with terminator
        self.T = self.T + TTT # add new string to old
        while pos < len(self.T): # add suffixes
            self.suffixes.append(self.T[pos:]+" "+str(pos))
            pos += 1
    def sort(self):
        temp=sorted(self.suffixes)
        self.suffixes=temp
    def print(self,ff=sys.stdout.fileno()):
        try:
            with open (ff,"w") as fout:
                for suffix in self.suffixes:
                    print(suffix,file=fout)
        except IOError:
            print("Cannot write file:",file=sys.stderr)
            print(ff,file=sys.stderr)

class FileParser:
    def __init__(self,filename):
        self.fn = filename
        self.DEFLINE='>' # fasta indicator
        self.NEWLINE="\n\r"  # universal
        self.READONLY="r"
        self.WRITE="w"
    def parse_text(self):
        ss = SortedSuffixes()
        try:
            with open(self.fn,self.READONLY) as f:
                for line in f:
                    if (len(line)>0):
                        #ss.add_string(line)  # debug
                        ss.set_string(line)  # debug
        except IOError:
            print("Cannot read file:",file=sys.stderr)
            print(self.fn,file=sys.stderr)
        return ss
    def parse_fasta(self):
        ss = SortedSuffixes()
        try:
            with open(self.fn,self.READONLY) as f:
                fastT=""
                for line in f:
                    T=line.rstrip(self.NEWLINE)
                    if T[0]==self.DEFLINE:
                        if len(fastT)>0:
                            ss.add_string(fastT)
                        fastT=""
                    else:
                        fastT = fastT + T
                if len(fastT)>0:
                    ss.add_string(fastT)
        except IOError:
            print("Cannot read file:",file=sys.stderr)
            print(self.fn,file=sys.stderr)
        return ss

def demo():
    T="abaababaabaab"
    ss = SortedSuffixes()
    ss.set_string(T)
    ss.sort()
    ss.print()

def args_parse():
    '''Parse command-line arguments.'''
    global args
    description="Output sorted suffixes of input string(s)."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'inFile', help='Input path/file.txt', type=str)
    parser.add_argument(
        'outFile', help='Output file prefix', type=str)
    parser.add_argument(
        '--fasta', help='Assume fasta file.',
        action='store_true')
    parser.add_argument(
        '--nosort', help='No sort. Useful for debug.',
        action='store_true')
    parser.add_argument(
        '--demo', help='Run demo program.',
        action='store_true')
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    '''Print all sorted suffixes.'''
    try:
        args_parse()
        if (args.demo):
            demo()
        else:
            fp=FileParser(args.inFile)
            if (args.fasta):
                ss=fp.parse_fasta()
            else:
                ss=fp.parse_text()
            if not args.nosort:
                ss.sort()
            ss.print(args.outFile)
    except Exception:
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
