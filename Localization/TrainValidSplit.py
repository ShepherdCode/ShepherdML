from random import Random
from csv import reader
import numpy as np
import math

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
    def load_sequences(self,gencode_file,labeled_genes):
        gid_tid = []
        with open(gencode_file,'r') as gencode:
            header = None
            csv = reader(gencode)
            rownum = 0
            for row in csv:
                if header is None:
                    header = row
                else:
                    tran_id = row[0]
                    gene_id = row[1]
                    if gene_id in labeled_genes:
                        label = labeled_genes[gene_id]
                        gid_tid.append((gene_id,tran_id,rownum,label))
                rownum += 1
        return gid_tid
    def load_labels(self,atlas_file,cells,RCI_THRESHOLDS):
        gene_labels = {}
        with open(atlas_file,'r') as atlas:
            header = None
            genes_considered = 0
            positives = 0
            csv = reader(atlas)
            for row in csv:
                if header is None:
                    header = row
                else:
                    gene = row[0]
                    genes_considered += 1
                    rci = float(row[1+cells])
                    if not math.isnan(rci):
                        # GENERATE BINARY LABELS
                        if rci < RCI_THRESHOLDS[0]:
                            label = 0   # RCI<0
                            gene_labels[gene]=label
                        elif rci >= RCI_THRESHOLDS[1]:
                            label = 1   #RCI>=0
                            positives += 1
                            gene_labels[gene]=label
        # TO DO: make this optional or not necessary
        print('ATLAS Genes:',genes_considered)
        print('Filter for genes with RCI labels from cell line:',cells)
        print('Labeled genes:',len(gene_labels))
        print('Positive labels:',positives)
        return gene_labels
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