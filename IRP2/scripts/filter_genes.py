'''
Prepare to decide which genes pass filter.
    Count read pairs such that the following could be applied.
The filter is applied to reads from homozygous tissue
    mapped to the heterozygous reference (true and false parent).
The filter involves 3 thresholds:
    min 50 reads per replicate
    min 200 reads per sample
    min 5-fold mapping preference for true parent.
    
This code was called filter_genes.py, then filter_homozygous.py, then back to filter_genes.py
'''
import argparse
import traceback
import csv
import sys

# PSEUDOCOUNT=1

class Replicate ():
    """
    Data struct.
    Holds input counts per gene for one biological replicate of one sample.
    """
    def __init__(self,name,sample,true_parent):
        self.my_name=name
        self.my_sample=sample
        self.true_parent=true_parent
        self.good_count_per_gene={}
        self.bad_count_per_gene={}
    def add(self,gene,allele,count):
        if gene not in self.good_count_per_gene.keys():
            self.good_count_per_gene[gene]=0
            self.bad_count_per_gene[gene]=0
        if (allele==self.true_parent):
            self.good_count_per_gene[gene] += count
        else:
            self.bad_count_per_gene[gene] += count
    def report(self):
        recs=str(len(self.good_count_per_gene))
        say("Replicate: "+self.my_name+" true parent: "+self.true_parent+" genes: "+recs)
    def get_name(self):
        return self.my_name
    def get_sample(self):
        return self.my_sample
    def get_genes(self):
        return self.good_count_per_gene.keys()
    def get_good_count(self,gene):
        return self.good_count_per_gene[gene]
    def get_bad_count(self,gene):
        return self.bad_count_per_gene[gene]
    def get_total_count(self,gene):
        total = self.good_count_per_gene[gene]
        total += self.bad_count_per_gene[gene]
        return total

class AllInput ():
    """
    Data struct.
    Holds input collection of Replicates
    """
    def __init__(self):
        self.replicates={}
    def set_basic_facts(self,keystr,sample,parent):
        if keystr not in self.replicates.keys():
            self.replicates[keystr]=Replicate(keystr,sample,parent)
    def add_count(self,keystr,gene,allele,count):
        rep = self.replicates[keystr]
        rep.add(gene,allele,count)
    def report(self):
        for keyval in self.replicates.values():
            keyval.report()
    def get_keys(self):
        return self.replicates.keys()
    def get_value(self,keystr):
        return self.replicates[keystr]

class Gene_Statistics ():
    """
    Data struct.
    Holds computed statistics for one gene.
    """
    def __init__(self,gene):
        self.gene=gene
        self.lowest_one_replicate=-1
        self.lowest_one_sample={}
        self.total_correct=0
        self.total_incorrect=0
    def get_fold (self):
        no= 1.0*max(1,self.total_incorrect)
        yes=1.0*self.total_correct
        fold = (yes-no)/no
        return fold
    def get_lowest_rep (self):
        return self.lowest_one_replicate
    def include (self, sample, count, good, bad):
        prev = self.lowest_one_replicate
        if prev < 0 or prev > count:
            self.lowest_one_replicate = count
        if sample not in self.lowest_one_sample.keys():
            self.lowest_one_sample[sample]=0
        self.lowest_one_sample[sample] += count
        self.total_correct += good
        self.total_incorrect += bad
    def get_lowest_samp (self):
        return min(self.lowest_one_sample.values())

class Calculator ():
    """
    Data analysis.
    Compute statistics for every gene across every replicate.
    """
    def __init__(self,dataset):
        self.dataset = dataset
        self.stats_per_gene = {}
    def analyze (self):
        for keystr in self.dataset.get_keys():
            replicate = self.dataset.get_value(keystr)
            sample = replicate.get_sample()
            for gene in replicate.get_genes():
                rep_count = replicate.get_total_count(gene)
                good = replicate.get_good_count(gene)
                bad = replicate.get_bad_count(gene)
                if gene not in self.stats_per_gene.keys():
                    self.stats_per_gene[gene]=Gene_Statistics(gene)
                stat = self.stats_per_gene[gene]
                stat.include(sample,rep_count,good,bad)
    def report (self):
        print("gene\tMinOneRep\tMinOneSamp\tFold")
        for gene in self.stats_per_gene.keys():
            gene_stat = self.stats_per_gene[gene]
            print("%s\t%d\t%d\t%f" % ( gene,
                gene_stat.get_lowest_rep(),
                gene_stat.get_lowest_samp(),
                gene_stat.get_fold()))

class Data_Loader (object):
    """
    Populate data structures from files.
    """
    def __init__(self,model):
        self.model = model
        self.accum = AllInput()

    def validate_model(self,filename):
        valid = True
        try:
            with open(filename, 'r') as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter="\t")
                rownum=0
                for row in reader: # row is a dict
                    if valid:
                        rownum += 1
                        if 'FilterSample' not in row.keys():
                            say("Error: model is missing FilterSample in row "+rownum)
                            valid = False
                        if 'FilterRep' not in row.keys():
                            say("Error: model is missing FilterRep in row "+rownum)
                            valid = False
                        if 'FilterParent' not in row.keys():
                            say("Error: model is missing FilterParent in row "+rownum)
                            valid = False
                        if 'CountsFilename' not in row.keys():
                            say("Error: model is missing CountsFilename in row "+rownum)
                            valid = False
        except IOError:
            say("Could not read file: "+filename)
        if rownum==0:
                say("Error: model file does not contain expected headers")
                valid = False
        return valid

    def input_counts(self,filename):
        try:
            with open(filename, 'r') as tsvfile:
                reader = csv.DictReader(tsvfile, delimiter="\t")
                for row in reader: # row is a dict
                    sample = row['FilterSample']
                    rep    = row['FilterRep']
                    parent = row['FilterParent']
                    filename = row['CountsFilename']
                    self.load_more_counts(sample,rep,parent,filename)
        except IOError:
            say("Could not read file: "+filename)
            return False
        return True

    def load_more_counts(self,sam,rep,par,filename):
        keystr = str(sam)+'-'+str(rep)
        self.accum.set_basic_facts(keystr,sam,par)
        with open(filename, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader: # row is a dict
                gene = str(row['gene'])
                allele = str(row['allele'])
                count = int(row['count'])
                self.accum.add_count(keystr,gene,allele,count)

    def summarize_inputs(self):
        self.accum.report()

    def load_all(self):
        say("Validate model from "+self.model)
        if not self.validate_model(self.model):
            raise Exception("Invalid model.")
        say("Input counts from files")
        if not self.input_counts(self.model):
            raise Exception("Invalid inputs.")
        say("Input data summary:")
        self.summarize_inputs()
        return self.accum

class Tests ():
    def __init__(self,filename):
        self.filename=filename
    def write_test_files(self):
        with open(self.filename, 'wt') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            tsv_writer.writerow(['FilterSample','FilterRep','FilterParent','Library','SeqDate','Parent','Stage','RepName','SeqName','CountsFilename'])
            tsv_writer.writerow(['1','1','Col0','Lib1','2019','Col-0','Early','BR1','NameOfSequencingRun1','testfile1.counts.tsv'])
            tsv_writer.writerow(['1','1','Col0','Lib1','2020','Col-0','Early','BR1','NameOfSequencingRun2','testfile2.counts.tsv'])
            tsv_writer.writerow(['1','2','Col0','Lib2','2019','Col-0','Early','BR2','NameOfSequencingRun3','testfile3.counts.tsv'])
        with open('testfile1.counts.tsv', 'wt') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            tsv_writer.writerow(['gene','allele','count','indel-free','spliced'])
            tsv_writer.writerow(['GENE1','Col0','25','0','0'])
            tsv_writer.writerow(['GENE1','Tsu1','25','0','0'])
        with open('testfile2.counts.tsv', 'wt') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            tsv_writer.writerow(['gene','allele','count','indel-free','spliced'])
            tsv_writer.writerow(['GENE1','Col0','25','0','0'])
            tsv_writer.writerow(['GENE1','Tsu1','25','0','0'])
        with open('testfile3.counts.tsv', 'wt') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            tsv_writer.writerow(['gene','allele','count','indel-free','spliced'])
            tsv_writer.writerow(['GENE1','Col0','25','0','0'])
            tsv_writer.writerow(['GENE1','Tsu1','25','0','0'])

def say(statement):
    if (args.debug):
        print(statement, file=sys.stderr, end="\n")

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='Map reads and revise the consensus.')
    parser.add_argument('model', help = 'Defines the inputs (tsv)', type=str)
    parser.add_argument('--test', help = 'Generate sample input files', action = 'store_true')
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 filter_homozygous.py --help
    """
    try:
        args_parse()
        if args.test:
            args.model='testmodel.tsv'
            tests = Tests(args.model)
            tests.write_test_files()
        loader = Data_Loader(args.model)
        loaded_data = loader.load_all()
        calculator = Calculator(loaded_data)
        calculator.analyze()
        calculator.report()
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
