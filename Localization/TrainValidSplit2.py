from random import Random
from csv import reader
import numpy as np
import statistics
from math import isnan
'''
Sample usage for training:
splitter = Splitter2()
splitter.load_counts_universe(Kmer_counts_filepath)
universe = splitter.get_gene_universe(atlas_train_file, cell_line_index)
partitions = splitter.get_train_valid_partitions(universe)
for p in partitions:
    train_genes,valid_genes = p[0],p[1]
    thresholds = splitter.choose_thresholds(train_genes:
    X_train,y_train = splitter.get_X_y(train_genes, thresholds)
    X_valid,y_valid = splitter.get_X_y(valid_genes, thresholds)
    
Sample usage for testing:
splitter = Splitter2()
splitter.load_counts_universe(Kmer_counts_filepath)
train_genes = splitter.get_gene_universe(atlas_train_file, cell_line_index)
thresholds = splitter.choose_thresholds(train_genes)
test_genes = splitter.get_gene_universe(atlas_test_file, cell_line_index)
X_test,y_test = splitter.get_X_y(test_genes, thresholds)
'''
class Splitter2():
    def __init__(self):
        self.ordered_gid_tid = []
        self.ordered_kmer_counts = []
        self.SEED = 42
        self.PARTITIONS=5
    def set_num_partitions(self,num):
        self.PARTITIONS=num
    def get_num_partitions(self):
        return self.PARTITIONS
    def get_train_valid_partitions(self, gene_to_rci:dict):  
        '''
        Returns list of pairs.
        Length of list = number of partitions.
        Each pair contains two dicts of gene_id to rci,
        one for training set, other for validation set.
        '''
        gene_list = list()
        for (gid,rci) in gene_to_rci.items():
            gene_list.append(gid)
        generator = Random()
        generator.seed(self.SEED)  # for self-consistency
        generator.shuffle(gene_list)  # in-place
        #
        partition_size = int(len(gene_list) * 1.0/self.PARTITIONS)
        partitions=[]
        for part in range(self.PARTITIONS):
            low = int(part*partition_size)
            high = int(low+partition_size)   
            valid_gene_to_rci={}
            train_gene_to_rci={}
            for pos in range(len(gene_list)):
                gene = gene_list[pos]
                rci=gene_to_rci[gene]
                if pos>=low and pos<high:
                    valid_gene_to_rci[gene]=rci
                else:
                    train_gene_to_rci[gene]=rci
            partitions.append((train_gene_to_rci,valid_gene_to_rci))
        return partitions
    def choose_thresholds(self, gene_to_rci:dict, method='mean'):
        '''
        The gene_to_rci should represent one training set.
        Returns a pair of float to be used as upper and lower threshold.
        When method=='mean', upper=lower=mean.
        When method=='one_z', upper=mean+stddev, lower=mean=stddev.
        '''
        vals=gene_to_rci.values()
        avg = statistics.mean(vals)
        if method=='one_z':
            std = statistics.stdev(vals)
            thresholds=(avg-std,avg+std)
        else:  # method=='mean'
            thresholds=(avg,avg)
        return thresholds
    def get_X_y(self, gene_to_rci:dict, thresholds:tuple):
        '''
        Construct two numpy arrays.
        Use previously-loaded counts universe for X values.
        Filter by genes in given dict whose values exceed thresholds.
        Use thresholds to make y values that are 0 or 1.
        '''
        X_list = []
        y_list = []
        for i in range(len(self.ordered_gid_tid)):
            (gene_id,tran_id)=self.ordered_gid_tid[i]
            if gene_id in gene_to_rci:
                rci = gene_to_rci[gene_id]
                label = None
                if rci<thresholds[0]:
                    label = 0
                elif rci>=thresholds[1]:
                    label = 1
                if label is not None:
                    kmer_counts=self.ordered_kmer_counts[i]
                    X_list.append(kmer_counts)
                    y_list.append(label)
        X_array = np.asarray(X_list)
        y_array = np.asarray(y_list)
        return X_array,y_array
    def get_gene_universe(self, atlas_file, cells):
        '''
        Returns dict of gene_id to RCI value,
        filtered for only genes with a numeric RCI in given cell line.
        '''
        rci_values = {}
        with open(atlas_file,'r') as atlas:
            header = None
            genes_considered = 0
            csv = reader(atlas)
            for row in csv:
                if header is None:
                    header = row
                else:
                    gene = row[0]
                    genes_considered += 1
                    rci = float(row[1+cells]) # e.g. choose value from 5th cell line
                    if not isnan(rci):
                        rci_values[gene]=rci
        print('Loaded values for cell line',cells)
        print('Selected',len(rci_values),'values out of',genes_considered)
        return rci_values
    def load_counts_universe(self, counts_file)->None:
        '''
        Load K-mer counts from a csv file. 
        Example file: CNRCI_coding_train_counts.K4.csv
        Save the data in list instance variables.
        Use the data later for get_X().
        '''
        self.ordered_gid_tid = []
        self.ordered_kmer_counts = []
        with open(counts_file,'r') as kc:
            header = None
            csv = reader(kc)
            for row in csv:
                if header is None:
                    header = row
                else:
                    gene_id = row.pop(0) 
                    tran_id = row.pop(0)
                    self.ordered_gid_tid.append( (gene_id,tran_id) )
                    kmercounts_one = np.asarray(row,dtype=np.int8)
                    self.ordered_kmer_counts.append(kmercounts_one)
