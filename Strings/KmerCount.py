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

class KmerHistFactory:
    def __init__(self,k):
        self.k = k
        self.hist = KmerHistogram(k)
        self.updates=0
    def get_histogram(self):
        return self.hist
    def count_updates(self):
        return self.updates

class KmerHistNaive(KmerHistFactory):
    def __get_range(self,T):
        n=len(T)
        k=self.k
        return range(0,n-k+1)
    def add_string(self,T):
        h=self.hist
        k=self.k
        range_of_i=self.__get_range(T)
        for i in range_of_i:
            kmer=T[i:i+k]
            h.increment_count(kmer)
            self.updates += 1
        last=T[-k:]
        h.add_last_kmer(last)

class KmerHistHarvester(KmerHistFactory):
    def extract_from(self,big_hist):
        self.updates=0
        self.__extract(big_hist)
        self.__endcase(big_hist)
        return self.updates
    def __extract(self,big):
        for (mer,freq) in big.get_items():
            length = len(mer)
            prefix = mer[0:length-1]
            self.hist.increment_count(prefix,freq)
            self.updates += 1
    def __endcase(self,big):
        h=self.hist
        for kmer in big.get_last_kmers():
            suffix = kmer[1:]
            h.increment_count(suffix)
            self.updates += 1
            h.add_last_kmer(suffix)

class FileProcessor:
    def __init__(self,filename,ksize):
        self.fn = filename
        self.k = ksize
        self.hist = None
        self.DEFLINE='>' # fasta indicator
        self.NEWLINE="\n\r"  # universal
        self.READONLY="r"
        self.WRITE="w"
    def parse_text(self):
        factory=KmerHistNaive(self.k)
        with open(self.fn,self.READONLY) as f:
            for line in f:
                T=line.rstrip(self.NEWLINE)
                if (len(T)>0):
                    factory.add_string(T)
        self.hist=factory.get_histogram()
        return factory.count_updates()
    def parse_fasta(self):
        factory=KmerHistNaive(self.k)
        with open(self.fn,self.READONLY) as f:
            fastT=""
            for line in f:
                T=line.rstrip(self.NEWLINE)
                if T[0]==self.DEFLINE:
                    if len(fastT)>0:
                        factory.add_string(fastT)
                    fastT=""
                else:
                    fastT = fastT + T
            if len(fastT)>0:
                factory.add_string(fastT)
        self.hist=factory.get_histogram()
        return factory.count_updates()
    def get_histogram(self):
        return self.hist
    def print_histogram(self,prefix,hist):
        k=str(hist.get_k())
        fn=prefix+"_"+k+".txt"
        with open(fn,self.WRITE) as f:
            hist.print(f)

def demo():
    T="abaababaabaab"
    factory=KmerHistFactory(2)
    factory.add_string(T)
    hist2=factory.get_histogram()
    updates2 = factory.count_updates()
    factory=None
    factory=KmerHistFactory(3)
    factory.add_string(T)
    hist3=factory.get_histogram()
    updates3 = factory.count_updates()
    factory=None
    factory=KmerHistHarvester(hist3)
    histE=factory.get_histogram()
    updatesE = factory.count_updates()
    factory=None
    print("Input string:")
    print(T)
    print("K=3, algorith=naive, updates=%d"%updates3)
    print("K=2, algorith=naive, updates=%d"%updates2)
    updatesC=updates3+updates2
    print("Combined, algorithm=naive, udates=%d"%updatesC)
    print("K=3, algorith=naive, updates=%d"%updates3)
    print("K=2, algorith=harvest, updates=%d"%updatesE)
    updatesC=updates3+updatesE
    print("Combined, algorithm=naive, udates=%d"%updatesC)

def args_parse():
    '''Parse command-line arguments.'''
    global args
    description='Count k-mers.'
    parser = argparse.ArgumentParser(description)
    parser.add_argument(
        'inFile', help='Input path/file.txt', type=str)
    parser.add_argument(
        'outFile', help='Output file prefix', type=str)
    parser.add_argument(
        '--k', help='Size of k-mer (5).',
        type=int, default=5)
    parser.add_argument(
        '--fasta', help='Assume fasta file.',
        action='store_true')
    parser.add_argument(
        '--demo', help='Run demo program.',
        action='store_true')
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    '''Count k-mers in strings.'''
    try:
        args_parse()
        if (args.demo):
            demo()
        else:
            k1=args.k
            fp=FileProcessor(args.inFile,k1)
            cnt1=0
            if args.fasta:
                cnt1=fp.parse_fasta()
            else:
                cnt1=fp.parse_text()
            hist1=fp.get_histogram()
            fp.print_histogram(args.outFile,hist1)
            print("Source=%s, k=%d, updates=%d"%
                (args.inFile,k1,cnt1))
            k2=k1-1
            factory=KmerHistHarvester(k2)
            cnt2=factory.extract_from(hist1)
            hist2=factory.get_histogram()
            fp.print_histogram(args.outFile,hist2)
            print("Source=%s, k=%d, updates=%d"%
                ("histogram",k2,cnt2))
    except Exception:
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
