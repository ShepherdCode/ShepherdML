import argparse
import traceback
import os
import sys
import numpy

class Stat_Machine():
    def __init__(self,infile,debug=False):
        self.infile = infile
        self.debug=debug

    def readin(self):
        with open(self.infile, 'r') as infa:
            numbers = []
            for line in infa:
                line=line.rstrip()
                value=int(line)
                numbers.append(value)
        return numbers

    def print_stats(self,ary):
        q1,q2,q3 = numpy.quantile(ary,[0.25,0.50,0.75])
        print("Count %d"%len(ary))
        print("Min %d"%numpy.amin(ary))
        print("Quartile1 %d"%q1)
        print("Median %d"%q2)
        print("Quartile3 %d"%q3)
        print("Max %d"%numpy.amax(ary))

    def print_lnrna_counts(self,ary):
        counts=[]
        limits=[]
        limits.append((0,200))
        counts.append(0)
        limits.append((200,1000))
        counts.append(0)
        limits.append((1000,2000))
        counts.append(0)
        limits.append((2000,3000))
        counts.append(0)
        limits.append((3000,30000))
        counts.append(0)
        for num in ary:
            for i in range(len(limits)):
                lim = limits[i]
                if num>=lim[0] and num<lim[1]:
                    counts[i] += 1
        for i in range(len(limits)):
            lim = limits[i]
            cim = counts[i]
            print("From_%d_to_%d %d"%(lim[0],lim[1],cim))

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (txt one num per line)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''FASTA input multiliners output oneliners.'''
    try:
        args_parse()
        s = Stat_Machine(args.infile,args.debug)
        c = s.readin()
        s.print_stats(c)
        s.print_lnrna_counts(c)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
