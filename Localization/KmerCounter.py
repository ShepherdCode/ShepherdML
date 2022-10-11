import numpy as np
import regex

class KmerCounter():
    NUCLEOTIDE_BITS={'A':0, 'C':1, 'G':2, 'T':3}
    TOTAL_NUCLEOTIDES=4
    BITS_PER_NUCLEOTIDE=2
    COUNT_TYPE=np.uint8
    MAX_COUNT = 2**8-1
    def __init__(self,K=4):
        self.setK(K)
    def setK(self,K):
        self.K=K
        self.VOCABULARY_SIZE = KmerCounter.TOTAL_NUCLEOTIDES**K
        # for K=4, we want MASK 00 11 11 11
        mask='00'
        for i in range(1,K):
            mask = mask + '11'
        self.MASK=int(mask,2)
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
    def subsequent_value(self,next_letter,prev_value):
        '''
        An optimization used by count_by_letter().
        '''
        bits = KmerCounter.NUCLEOTIDE_BITS[next_letter]  # protected against N
        masked = prev_value & self.MASK   # erase 2 left-most bits of hash_value
        shifted = masked << 2             # make room for 2 bits
        next_value = shifted + bits       # add 2 new bits at right of hash_value
        return next_value
    def seq_to_kmer_counts(self,seq,method=3):
        '''
        Clients should call setK() once,
        then seq_to_kmer_counts for each DNA sequence.
        Returns a numpy array of one-bit ints,
        so maximum count is 127.
        For example, after setK(4) and counts=seq_to_kmer_counts('AAAAA'),
        then counts(hash_value('AAAA')) will be 2.
        The default method is recommended.
        '''
        if method==1:
            return self.count_using_regex(seq)
        elif method==2:
            return self.count_by_letter(seq)
        else:   # 3 default
            return self.count_by_token(seq)
    def count_using_regex(self,seq):
        '''
        This is a reference implementation.
        '''
        counts = np.zeros(
            self.VOCABULARY_SIZE, KmerCounter.COUNT_TYPE)
        token = 'A'*self.K
        i=0
        while token is not None:
            # string.count(substr) does not count overlapping substrings!
            # counts[i] = seq.count(token)
            actual = len(regex.findall(token, seq, overlapped=True))
            capped = min(actual,KmerCounter.MAX_COUNT)
            counts[i] = capped
            i += 1
            token=self.next_kmer(token)
        return counts
    def count_by_token(self,seq):
        counts = np.zeros(
            self.VOCABULARY_SIZE, KmerCounter.COUNT_TYPE)
        for p in range(len(seq)-self.K+1):
            next_token=seq[p:p+self.K]
            hash_value = self.hash_value(next_token)
            if hash_value is not None:
                if counts[hash_value]<KmerCounter.MAX_COUNT:
                    counts[hash_value] += 1
        return counts
    def count_by_letter(self,seq):
        '''
        Optimized way to count K-mers in DNA.
        When hashing a K-mer, re-use the bits from the last K-1 letters.
        Given 'CCCCT', after hashing 'CCCC', just add the bit-value of T.
        On first K-mer, and on K-mers with N, must revert to full count.
        In tests on mRNA, this was no faster than count_by_token.
        '''
        prev_value = None
        next_letter = None
        counts = np.zeros(
            self.VOCABULARY_SIZE, KmerCounter.COUNT_TYPE)
        for p in range(len(seq)-self.K+1):
            next_letter = seq[p+self.K-1]
            if prev_value is not None and \
            next_letter in KmerCounter.NUCLEOTIDE_BITS:
                hash_value = self.subsequent_value(next_letter,prev_value)
                prev_value = hash_value
            else:
                next_token=seq[p:p+self.K]
                hash_value = self.hash_value(next_token)
                prev_value = hash_value
            if hash_value is not None:
                if counts[hash_value]<KmerCounter.MAX_COUNT:
                    counts[hash_value] += 1
        return counts
    def next_kmer(self,token)->str:
        '''
        Designed to be called in a feed-forward loop.
        Example: next_kmer('AAAA') returns 'AAAC'.
        Example: next_kmer('ATTT') returns 'CAAA'.
        Example: next_kmer('TTTT') returns None.
        '''
        p = len(token)-1
        next_token = list(token)
        while p>=0:
            if token[p]=='A':
                next_token[p]='C'
                word=''.join(next_token)
                return word
            if token[p]=='C':
                next_token[p]='G'
                word=''.join(next_token)
                return word
            if token[p]=='G':
                next_token[p]='T'
                word=''.join(next_token)
                return word
            if token[p]=='T':
                next_token[p]='A'
                p = p-1
        return None