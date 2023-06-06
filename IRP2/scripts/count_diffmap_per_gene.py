"""
Count the number of differentially mapped read pairs per gene+allele.
Input: *.out files written by differential_mapping.py
Input format: one per line, no headers, space separated values
    K00253:118:HCLGFBBXY:4:1101:14072:1156 AT1G43580.1 Tsu1 SNP 2 5 1 4
Output format: one per line, with headers, tab separated values
    gene allele pairs
    AT1G01020.1 Tsu1 7
This program replaces make_tsv.sh
"""
import traceback
import argparse
import sys
import csv

class Input_Line:
    def __init__(self,readname,gene,allele):
        self.readname=readname
        self.gene=gene
        self.allele=allele
    def set_evidence(self,evidence_type,ev1,ev2,ev3,ev4):
        self.evidence_type=evidence_type
        self.ev1=ev1
        self.ev2=ev2
        self.ev3=ev3
        self.ev4=ev4
    def get_gene(self):
        return self.gene
    def get_allele(self):
        return self.allele

class Gene_Count:
    def __init__(self,gene):
        self.gene=gene
        self.alleles={}
    def increment(self,allele):
        if allele not in self.alleles.keys():
            self.alleles[allele]=0
        self.alleles[allele] = self.alleles[allele]+1
    def get_gene(self):
        return self.gene
    def get_allele_counts(self):
        return self.alleles.items()

class All_Counts:
    def __init__(self):
        self.gene_counts={}
    def increment(self,input_line):
        gene=input_line.get_gene()
        allele=input_line.get_allele()
        if gene not in self.gene_counts.keys():
            self.gene_counts[gene]=Gene_Count(gene)
        gc=self.gene_counts[gene]
        gc.increment(allele)
    def get_genes(self):
        return self.gene_counts.keys()
    def get_gene_counts(self,gene):
        return self.gene_counts[gene]

class Streamer:
    def __init__(self, filename):
        self.infile=filename
        self.counts=All_Counts()
    def print(self):
        print("gene\tallele\tcount")
        all_genes = self.counts.get_genes()
        for one_gene in all_genes:
            one_count = self.counts.get_gene_counts(one_gene)
            for pair in one_count.get_allele_counts():
                allele = pair[0]
                cnt = pair[1]
                print("%s\t%s\t%d" % (one_gene, allele, cnt))
    def stream(self):
        with open(self.infile) as csvfile:
            fn=['read','gene','allele','evtype','ev1','ev2','ev3','ev4']
            reader = csv.DictReader(csvfile,delimiter=' ',fieldnames=fn)
            for row in reader:
                readname=row['read']
                gene=row['gene']
                allele=row['allele']
                inline = Input_Line(readname,gene,allele)
                inline.set_evidence(row['evtype'],row['ev1'],row['ev2'],row['ev3'],row['ev4'])
                self.counts.increment(inline)

def say(statement):
    if (args.debug):
        print(statement, file=sys.stderr, end="\n")

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='Count diff maps per gene.')
    parser.add_argument('infile', help = 'Text file from differential_mapping.py', type=str)
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 count_diffmap_per_gene.py --help
    """
    try:
        args_parse()
        analyzer = Streamer(args.infile)
        analyzer.stream()
        analyzer.print()
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Consider running with --debug')
