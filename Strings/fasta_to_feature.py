import argparse
import traceback
import sys
import csv

'''
Input: FASTA
    Expect defline starts with '>'
    The rest of the defline is ignored.
    Expect entire sequence on one line.
Output: CSV
    Header line explains the columns.
    Col 1 = sequence number (starts at 1 by default)
    Col 2 = sequence length
    Col 3... count for an individual K-mer, possibly zero.
'''

class FeatureVector:
    '''Keep count of every k-mer.'''
    def __init__(self,k):
        self.ALPHABET_SIZE=4 ## assumed
        self.k=k
        self.max=pow(self.ALPHABET_SIZE,self.k)
        self.counts=[0]*self.max
        self.digits={'A':0,'C':1,'G':2,'T':3}
    def kmer_to_int(self,kmer):
        '''Compute K-mer's position in array or alphabetical list.'''
        value=0
        power=1
        for i in range(self.k-1,-1,-1):
            letter = kmer[i]
            digit = self.digits[letter]
            value += power * digit
            power = power * self.ALPHABET_SIZE
        return value
    def increment_count(self,kmer,cnt=1):
        '''Given K-mer and count, add to existing count.'''
        if (len(kmer)!=self.k):
            raise Exception("Invalid k-mer: "+kmer)
        kmerid = self.kmer_to_int(kmer)
        self.counts[kmerid] += cnt
    def get_array(self):
        '''Return array of counts in alphabetical order.'''
        return self.counts
    def get_names(self):
        '''Return list of all possible K-mers in alphabetical order.'''
        letters=['A','C','G','T']
        names=[]
        for i in range(self.max):
            kmer=""
            number = i
            for pos in range(self.k-1,-1,-1):
                place=pow(self.ALPHABET_SIZE,pos)
                digit=0
                while number>=place:
                    digit += 1
                    number -= place
                kmer += letters[digit]
            names.append(kmer)
        return names
    def next_size_down(self):
        '''Generate new FeatureVector for k-1.'''
        newfv=FeatureVector(self.k-1)
        countup=0
        target=0
        for count in range(len(self.counts)):
            newfv.counts[target]+=count
            countup += 1
            if (countup==self.ALPHABET_SIZE):
                target += 1
                countup=0
        return newfv

class FileProcessor:
    def __init__(self,infile,outprefix,kmax,kmin):
        self.wordsize = kmax
        self.kmin=kmin
        self.infile = infile
        self.outfile = outprefix
        self.DEFLINE='>' # fasta indicator
        self.NEWLINE="\n\r"  # universal
        self.READONLY="r"
        self.WRITE="w"
        self.features = FeatureVector(self.wordsize)
        self.first_number=1
    def set_first_number(self,basenum):
        self.first_number=basenum
    def open_files(self):
        # Waste a little space for clarity.
        # The handle for k=3 is at handles[3].
        self.handles=[None]*(1+self.wordsize)
        self.writers=[None]*(1+self.wordsize)
        for k in range(self.kmin,self.wordsize+1):
            filename = "%s.%dmer.features.csv"%(self.outfile,k)
            csvfile= open(filename,self.WRITE,newline='')
            writer = csv.writer(csvfile, delimiter=',')
            self.handles[k]=csvfile
            self.writers[k]=writer
            fvec = FeatureVector(k)
            feature_names = fvec.get_names()
            header=['seqnum','seqlen']
            header += feature_names
            writer.writerow(header)
    def close_files(self):
        for handle in handles:
            handle.close()
    def make_features(self):
        '''Make csv file from a oneline fasta file.'''
        with open(self.infile,self.READONLY) as infile:
            is_defline=False
            seqnum=self.first_number-1
            for line in infile:
                T=line.rstrip(self.NEWLINE)
                if T[0]==self.DEFLINE:
                    if is_defline:
                        print(T)
                        raise Exception('Unexpected defline. Is this a oneline fasta? '+str(seqnum))
                    is_defline = True
                    seqnum += 1
                else:
                    if not is_defline:
                        print(T)
                        raise Exception('Unexpected sequence. Is this a oneline fasta? '+str(seqnum))
                    is_defline = False
                    seqlen=len(T)
                    self.extract_kmer_counts(T)
                    self.output_features(seqnum,seqlen,self.wordsize,self.features)
                    self.harvester(seqnum,seqlen)

    def output_features(self,seqnum,seqlen,k,features):
        vec=features.get_array()
        row=[seqnum,seqlen]
        row += vec
        writer=self.writers[k]
        writer.writerow(row)

    def harvester(self,seqnum,seqlen):
        '''Extract n-mers for n<K.'''
        k=self.wordsize
        fv=self.features
        while k>self.kmin:
            k -= 1
            fv = fv.next_size_down()
            self.output_features(seqnum,seqlen,k,fv)

    def extract_kmer_counts(self,seq):
        '''Walk the sequence. Increment count per K-mer.'''
        # Assume not cumulative:
        k=self.wordsize
        self.features = FeatureVector(k)
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
        '--base', help='start numbering at (1)', type=int)
    parser.add_argument(
        '--kmax', help='Max k-mer size (4).',
        type=int, default=4)
    parser.add_argument(
        '--kmin', help='Min k-mer size (kmax).',
        type=int, default=4)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    '''Translate fasta to feature file.'''
    try:
        args_parse()
        kmax=args.kmax
        kmin=args.kmin
        if kmin>kmax:
            kmin=kmax
        fp=FileProcessor(args.inFile,args.outPrefix,kmax,kmin)
        fp.open_files()
        if args.base is not None:
            basenum = int(args.base)
            fp.set_first_number(basenum)
        fp.make_features()
    except Exception:
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
