import argparse
import traceback
import sys
import csv

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
        self.counts[kmerid] += cnt
    def get_array(self):
        return self.counts
    def get_names(self):
        letters=['A','C','G','T']
        names=[]
        for i in range(self.max):
            kmer=""
            number = i
            for pos in range(self.k-1,-1,-1):
                place=pow(4,pos)
                digit=0
                while number>=place:
                    digit += 1
                    number -= place
                kmer += letters[digit]
            names.append(kmer)
        return names

class FileProcessor:
    def __init__(self,infile,outprefix,wordsize):
        self.wordsize = wordsize
        self.infile = infile
        self.outfile = outprefix+".features.csv"
        self.DEFLINE='>' # fasta indicator
        self.NEWLINE="\n\r"  # universal
        self.READONLY="r"
        self.WRITE="w"
        self.features = FeatureVector(wordsize)
        self.uniform_label = None
    def set_uniform_label(self,label):
        self.uniform_label = label
    def make_features(self):
        '''Make features from a oneline fasta file.'''
        with open(self.outfile,self.WRITE,newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            feature_names = self.features.get_names()
            header=['label','seqname']
            header += feature_names
            writer.writerow(header)
            with open(self.infile,self.READONLY) as infile:
                is_defline=False
                seqname=""
                label=""
                for line in infile:
                    T=line.rstrip(self.NEWLINE)
                    if T[0]==self.DEFLINE:
                        if is_defline:
                            raise Exception('Defline. Is this a oneline fasta? '+str(seqnum))
                        is_defline = True
                        fields=T[1:].split(' ')
                        seqname=fields[0]
                        if self.uniform_label is None:
                            label=fields[-1]
                        else:
                            label=self.uniform_label
                    else:
                        if not is_defline:
                            raise Exception('Sequence. Is this a oneline fasta? '+str(seqnum))
                        is_defline = False
                        self.extract_kmer_counts(T)
                        self.process_seq(label,seqname,writer)

    def process_seq(self,label,seqname,writer):
        vec=self.features.get_array()
        row=[label,seqname]
        row += vec
        writer.writerow(row)
    def extract_kmer_counts(self,seq):
        k=self.wordsize
        n = len(seq)
        for i in range(n-k+1):
            kmer=seq[i:i+k]
            self.features.increment_count(kmer)  ## needs work

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
        '--label', help='same label for all (int)', type=int)
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
        fp=FileProcessor(args.inFile,args.outPrefix,args.k)
        if (args.label is not None):
            fp.set_uniform_label(args.label)
        fp.make_features()
    except Exception:
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
