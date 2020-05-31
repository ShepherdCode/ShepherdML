from Bucket import Bucket

class BucketList():
    '''Utility class for Manber-Myers algorithm.'''

    def __init__(self,stringT):
        self.buckets=[]
        self.stringT = stringT

    def __str__(self):
        viz = ""
        for b in self.buckets:
            viz = viz+str(b)
            viz = viz+", "
        return viz

    def get_buckets(self):
        return self.buckets

    # Time and space considerations.
    # We need a mechanism to put similar prefixes in the same bucket.
    # We could use a linear-time bucket sort.
    # Here, we use a constant-time dictionary because it is easier to code.
    # The dictionary is not sorted by key but that is ok for now.
    # The dictionary size is not controlled but
    # we could cap dictionary size proportional to alphabet size.
    def string_to_buckets(self):
        '''Get suffixes of T. Bucketize by 1-letter prefixes.'''
        prefix_to_bucket_dictionary = {}
        prev_prefix = ''
        this_prefix = ''
        this_bucket = Bucket(this_prefix,self.stringT) # placeholder, not in dictionary
        self.buckets = []
        L = len(self.stringT)
        for i in range(0,L):
            this_prefix = self.stringT[i:i+1]
            this_suffix = self.stringT[i:L]
            this_index = i+1
            if this_prefix == prev_prefix:
                # Continue adding to same bucket
                pass
            elif this_prefix in prefix_to_bucket_dictionary:
                # Re-use one of the existing buckets
                this_bucket = prefix_to_bucket_dictionary[this_prefix]
            else:
                # Start a brand new bucket
                this_bucket = Bucket(this_prefix,self.stringT)
                prefix_to_bucket_dictionary[this_prefix] = this_bucket
                self.buckets.append(this_bucket)
            this_bucket.add_suffix(this_index)
            prev_prefix = this_prefix

    # Time and space considerations.
    # We need to sort the initial buckets by their one-letter prefixes.
    # We could have sorted them during bucket construction using
    # linear-time bucket sort.
    # We could still sort them using linear-time bucket sort.
    # Here, we use non-linear Python sort because it is easier to code.
    def sort_buckets(self):
        self.buckets.sort(key=Bucket.get_prefix)

    def sort_within_buckets_shortcut(self):
        for b in self.buckets:
            b.sort_suffixes_shortcut()
