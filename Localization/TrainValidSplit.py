from random import Random
import numpy as np
class Splitter():
    def __init__(self):
        self.gid_tid = []
        self.ordered_kmer_counts = []
        self.gene_list=[]
    def set_ids(self,gid_tid_row_label):
        self.gid_tid = gid_tid_row_label
    def set_counts(self,ordered_kmer_counts):
        self.ordered_kmer_counts = ordered_kmer_counts
    def randomize(self):
        gene_set = set()
        for (gid,tid,row,label) in self.gid_tid:
            gene_set.add(gid)
        generator = Random()
        generator.seed(42)
        self.gene_list = list(gene_set)
        generator.shuffle(self.gene_list)  # in-place
    def train_valid_split(self,iteration, partitions):
        # TO DO: make sure we don't leave out the last sequence due to rounding
        partition_size = int(len(self.gene_list) * 1.0/partitions)
        low = int(iteration*partition_size)
        high = int(low+partition_size)    
        valid_genes = set(self.gene_list[low:high])
        # TO DO: Here, we grow lists one at a time.
        # It would be faster with fixed-size numpy arrays?
        X_train = []   
        X_valid = []
        y_train = []
        y_valid = []
        for (gid,tid,row,label) in self.gid_tid:
            if gid in valid_genes:
                X_valid.append(self.ordered_kmer_counts[row])   
                y_valid.append(label)
            else:
                X_train.append(self.ordered_kmer_counts[row])   
                y_train.append(label)
        X_train = np.asarray(X_train)
        y_train = np.asarray(y_train)
        X_valid = np.asarray(X_valid)
        y_valid = np.asarray(y_valid)
        return X_train,y_train,X_valid,y_valid