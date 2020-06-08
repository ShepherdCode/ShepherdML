import argparse
import traceback
import sys

class KmerHistogram:
    def __init__(self,k):
        self.k=k
        self.counts={} # dict
        self.last_kmers=[] # list
    def get_k(self):
        return self.k
    def check_k(self,kmer):
        if len(kmer)!=self.k:
            print("Expected len "+str(self.k)+", received "+kmer)
            raise Exception
    def get_freq(self,kmer):
        self.check_k(kmer)
        if not kmer in self.counts:
            return 0
        return(self.counts[kmer])
    def increment_count(self,kmer,cnt=1):
        self.check_k(kmer)
        if kmer not in self.counts:
            self.counts[kmer]=0
        self.counts[kmer] += cnt
    def get_kmers(self):
        return self.counts.keys()
    def get_items(self):
        return self.counts.items()
    def add_last_kmer(self,kmer):
        self.last_kmers.append(kmer)
    def get_last_kmers(self):
        return(self.last_kmers)
    def print(self,ff=sys.stdout):
        for kmer in sorted(self.counts):
            print("%s = %d"%(kmer,self.counts[kmer]),file=ff)

class FileProcessor:
    def __init__(self,infile,outprefix):
        self.infile = infile
        self.outprefix = outprefix
        self.DEFLINE='>' # fasta indicator
        self.NEWLINE="\n\r"  # universal
        self.READONLY="r"
        self.WRITE="w"
    def make_features(self,wordsize,label):
        '''Make features from a oneline fasta file.'''
        with open(self.infile,self.READONLY) as f:
            is_defline=False
            seqnum=0
            for line in f:
                T=line.rstrip(self.NEWLINE)
                if T[0]==self.DEFLINE:
                    if is_defline:
                        raise Exception('Defline. Is this a oneline fasta? '+str(seqnum))
                    is_defline = True
                    seqname=T[1:]
                else:
                    if not is_defline:
                        raise Exception('Sequence. Is this a oneline fasta? '+str(seqnum))
                    is_defline = False
                    sequence = T
                    features=self.process(sequence,wordsize)
    def process(self,seq,k):
        features = []
        print("Working on "+seq)
        n = len(seq)
        for i in range(n):
            kmer=seq[i:i+k]
            features.append(kmer)  ## needs work
        return features

def args_parse():
    '''Parse command-line arguments.'''
    global args
    description='Convert fasta file to feature file.'
    parser = argparse.ArgumentParser(description)
    parser.add_argument(
        'inFile', help='path/name (oneline fasta)', type=str)
    parser.add_argument(
        'outPrefix', help='output file prefix (csv)', type=str)
    parser.add_argument(
        'label', help='class (int)', type=int)
    parser.add_argument(
        '--k', help='Size of k-mer (4).',
        type=int, default=5)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    '''Translate fasta to feature file.'''
    try:
        args_parse()
        fp=FileProcessor(args.inFile,args.outPrefix)
        fp.make_features(args.k,args.label)
    except Exception:
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
