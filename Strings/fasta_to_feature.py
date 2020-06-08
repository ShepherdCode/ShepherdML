import argparse
import traceback
import sys

class FeatureVector:
    '''Keep count of every k-mer.'''
    ALPHABET_SIZE=4 ## assumed
    def __init__(self,k):
        self.k=k
        self.max=pow(4,self.k)
        self.counts=[0]*self.max
        self.digits={'A':0,'C':1,'G':2,'T':3}
    def kmer_to_int(self,kmer):
        value=0
        power=1
        for i in range(self.k-1,-1,-1):
            letter = kmer[i]
            digit = self.digits[letter]
            value += power * digit
            power = power * 4
        return value
    def increment_count(self,kmer,cnt=1):
        if (len(kmer)!=self.k):
            raise Exception("Invalid k-mer: "+kmer)
        kmerid = self.kmer_to_int(kmer)
        print(kmer)
        print(kmerid)
        self.counts[kmerid] += cnt

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
        features = FeatureVector(k)
        print("Working on "+seq)
        n = len(seq)
        for i in range(n-k+1):
            kmer=seq[i:i+k]
            features.increment_count(kmer)  ## needs work
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
        type=int, default=4)
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
