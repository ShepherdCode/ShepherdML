from random import Random
class Splitter():
    def __init__(self,ordered_gid_tid,ordered_kmer_counts,ordered_labels):
        self.ordered_gid_tid = ordered_gid_tid
        self.ordered_kmer_counts = ordered_kmer_counts
        self.ordered_labels = ordered_labels
        gene_set = set()
        for (gid,tid) in self.ordered_gid_tid:
            gene_set.add(gid)
        generator = Random()
        generator.seed(42)
        self.gene_list = list(gene_set)
        generator.shuffle(self.gene_list)  # in-place
    def train_valid_split(self,iteration, partitions):
        # TO DO: do we sometimes leave out the last sequence?
        partition_size = int(len(self.gene_list) * 1.0/partitions)
        low = int(iteration*partition_size)
        high = int(low+partition_size)    
        valid_genes = set(gene_list[low:high])
        # TO DO: Here, we grow lists one at a time.
        # It would be faster with fixed-size numpy arrays.
        X_train = []   
        X_valid = []
        y_train = []
        y_valid = []
        for i in range(len(ordered_gid_tid)):
            (gid,tid) = ordered_gid_tid[i]     
            if gid in valid_genes:
                X_valid.append(ordered_kmer_counts[i])   
                y_valid.append(ordered_labels[i])
            else:
                X_train.append(ordered_kmer_counts[i])   
                y_train.append(ordered_labels[i])
        X_train = np.asarray(X_train)
        y_train = np.asarray(y_train)
        X_valid = np.asarray(X_valid)
        y_valid = np.asarray(y_valid)
        return X_train,y_train,X_valid,y_valid