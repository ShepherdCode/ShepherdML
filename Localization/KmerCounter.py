import numpy as np
class KmerCounter():
    NUCLEOTIDE_BITS={'A':0, 'C':1, 'G':2, 'T':3}
    TOTAL_NUCLEOTIDES=4
    BITS_PER_NUCLEOTIDE=2
    COUNT_TYPE=np.int8
    MAX_COUNT = 2**8-1
    def __init__(self,K=4):
        self.setK(K)
    def setK(self,K):
        self.K=K
        self.VOCABULARY_SIZE = KmerCounter.TOTAL_NUCLEOTIDES**K
    def get_vocabulary_size(self):
        return self.VOCABULARY_SIZE
    def hash_value(self,token):
        '''
        Get the counts array index for the given K-mer.
        Most clients won't call this directly.
        TO DO: For efficiency, cache the previous K-1 letters.
        For example, given AACC then ACCT, no need to recompute ACC.
        '''
        kmer_hash=0
        if len(token)!=self.K:
            raise Exception('Given token not of length K='+str(self.K))
        letters = list(token)
        for letter in letters:
            try:
                additional = KmerCounter.NUCLEOTIDE_BITS[letter]
                kmer_hash = kmer_hash << KmerCounter.BITS_PER_NUCLEOTIDE
                kmer_hash = kmer_hash + additional
            except KeyError:
                # Ignore tokens with N or any non-nucleotide
                return None
        return kmer_hash
    def seq_to_kmer_counts(self,seq):
        '''
        Clients should call setK() once,
        then seq_to_kmer_counts for each DNA sequence.
        Returns a numpy array of one-bit ints,
        so maximum count is 127.
        For example, after setK(4) and counts=seq_to_kmer_counts('AAAAA'),
        then counts(hash_value('AAAA')) will be 2.
        '''
        counts = np.zeros(
            self.VOCABULARY_SIZE, KmerCounter.COUNT_TYPE)
        for p in range(len(seq)-self.K+1):
            token=seq[p:p+self.K]
            hash_value = self.hash_value(token)
            if hash_value is not None:
                if counts[hash_value]<KmerCounter.MAX_COUNT:
                    counts[hash_value] += 1
        return counts