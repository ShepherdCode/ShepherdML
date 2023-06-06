"""
Apply IRP rules to each pair of read pair alignments.
Read pairs were mapped to artificial heterozygous genome.
Retained best 2 maps per pair.
To qualify as informative, pair must:
    Have one concordant mapping to each haplotype.
    Map to same locus in both haplotypes.
    Map significantly better to one haploytpe.
"""
import pysam
import traceback
#import subprocess
#import logging
#import pathlib
import argparse
import sys
#import pprint

MIN_MAPQ=5
MIN_READS=4
MAX_TEMPLATE=400
NUM_TRANSCRIPTS=1

class Report:
    total_alignments = 0
    wrong_num_mapping = 0
    genome_imbalance = 0
    low_mapq = 0
    discordant_mappings = 0
    different_transcripts = 0
    indels_both_genomes = 0
    total_indels = 0
    total_snps = 0

    def print():
        say("Report")
        say("%10d = Total alignments parsed"
            % Report.total_alignments)
        say("%10d = Wrong number of mappings. Require = %d"
            % (Report.wrong_num_mapping, MIN_READS))
        say("%10d = Imbalance of genomes. Require = 2+2"
            % Report.genome_imbalance)
        say("%10d = Low MapQ/uniqueness. Require >= %d"
            % (Report.low_mapq, MIN_MAPQ))
        say("%10d = Discordant mappings. Require < %d"
            % (Report.discordant_mappings, MAX_TEMPLATE))
        say("%10d = Different transcripts. Require = %d"
            % (Report.different_transcripts, NUM_TRANSCRIPTS))
        say("%10d = Indel to both genomes. Require < 2"
            % Report.indels_both_genomes)
        say("%10d = Informative pairs decided by indel."
            % Report.total_indels)
        say("%10d = Informative pairs decided by SNP."
            % Report.total_snps)
        totalir=Report.total_indels+Report.total_snps
        say("%10d = Total informative pairs."
            % totalir)

class AlignedRead:
    """
    Data structure holds 1 aligned read (without the mate pair).
    """
    def __init__(self):
        self.clear()
    def clear(self):
        self.gene=""
        self.genome=""
        self.cigar=""
        self.num_mismatch=0
        self.mapq=0
        self.has_indels=False
        self.template=0
    def reset(self,gene,cigar,nm,mapq,tlen):
        position=gene.index("_") # Expect AT3G61260.1_Tsu1
        # TO DO: raise exception otherwise
        self.gene=gene[0:position]
        self.genome=gene[position+1:]
        self.cigar=cigar
        self.num_mismatch=nm
        self.mapq=mapq
        self.template=tlen
        self.has_indels=("I" in cigar or "D" in cigar)
    def print(self):  # just for debugging
        print(" read %s %s %s %d %d %d" %
            (self.gene,self.genome,self.cigar,self.num_mismatch,
            self.mapq,self.template))

class ReadGroup:
    """
    Data structure holds 4 read alignments for one mate pair.
    Expect all four have same name.
    UiO mates have same ID.
    TO DO: For systems with /1 and /2 read suffixes, chop the sufix.
    Assume 2 reads that map to genomeA, and 2 that map to genomeB, are mates.
    TO DO: Use sequence IDs and mate ID to identify mates for sure.
    """
    SNP_RATIO=2 # believe indel-free pair map if #SNPs * 2 still better than 2nd best
    SNP_MAX=1 # believe indel-free pair map if #SNPs <= 1 and 2nd best has indels
    def __init__(self):
        self.reads=[AlignedRead(),AlignedRead(),AlignedRead(),AlignedRead()]
        self.clear()
    def clear(self):
        for read in self.reads:
            read.clear()
        self.num_reads=0
        self.read_name=""
    def set_name(self,name):
        self.read_name=name
    def add_read(self,gene,cigar,nmtag,mapq,tlen):
        Report.total_alignments += 1
        self.num_reads += 1
        if self.num_reads <= 4:
            this_read = self.reads[self.num_reads-1]
            this_read.reset(gene,cigar,nmtag,mapq,tlen)
            #this_read.print()  # just for debug
    def get_min_mapq(self):
        return min(self.reads[0].mapq,self.reads[1].mapq,
            self.reads[2].mapq,self.reads[3].mapq)
    def is_discordant(self):
        discord=False
        for r in self.reads:
            if abs(r.template)>MAX_TEMPLATE:
                discord=True
        return discord
    def is_balanced_genomes(self):
        genome_count={}
        for r in self.reads:
            if not r.genome in genome_count.keys():
                genome_count[r.genome]=1
            else:
                genome_count[r.genome] += 1
        list1=genome_count.values()
        balance = (len(list1)==2 and max(list1)==2 and min(list1)==2)
        return balance
    def get_genomes(self):
        genomes=[]
        for r in self.reads:
            if not r.genome in genomes:
                genomes.append(r.genome)
        return genomes
    def get_transcripts(self):
        transcripts=[]
        for r in self.reads:
            if not r.gene in transcripts:
                transcripts.append(r.gene)
        return transcripts
    def get_mates(self,genome):
        list=[]
        for r in self.reads:
            if r.genome == genome:
                list.append(r)
        return list
    def indels_both_genomes(self):
            (genomeA,genomeB)=self.get_genomes()
            pairA=self.get_mates(genomeA)
            pairB=self.get_mates(genomeB)
            indelA=(pairA[0].has_indels or pairA[1].has_indels)
            indelB=(pairB[0].has_indels or pairB[1].has_indels)
            return (indelA and indelB)
    def choose_genome_by_indel(self):
            (genomeA,genomeB)=self.get_genomes()
            pairA=self.get_mates(genomeA)
            pairB=self.get_mates(genomeB)
            indelA=(pairA[0].has_indels or pairA[1].has_indels)
            indelB=(pairB[0].has_indels or pairB[1].has_indels)
            snpA= pairA[0].num_mismatch + pairA[1].num_mismatch
            snpB= pairB[0].num_mismatch + pairB[1].num_mismatch
            if indelA and not indelB and snpB <= ReadGroup.SNP_MAX:
                return genomeB
            elif indelB and not indelA and snpA <= ReadGroup.SNP_MAX:
                return genomeA
            return None
    def choose_genome_by_snp(self):  # assuming no indels
            (genomeA,genomeB)=self.get_genomes()
            pairA=self.get_mates(genomeA)
            pairB=self.get_mates(genomeB)
            snpA= pairA[0].num_mismatch + pairA[1].num_mismatch
            snpB= pairB[0].num_mismatch + pairB[1].num_mismatch
            if ReadGroup.SNP_RATIO*snpA < snpB: # examples: 0 vs 1, 1 vs 3, 2 vs 5
                return genomeA
            elif snpA > ReadGroup.SNP_RATIO*snpB:
                return genomeB
            return None
    def output(self):
        #print(self.get_genomes())
        if (self.num_reads != MIN_READS):
            Report.wrong_num_mapping += 1
        elif self.get_min_mapq() < MIN_MAPQ:
            Report.low_mapq += 1
        elif self.is_discordant():
            Report.discordant_mappings += 1
        elif not self.is_balanced_genomes():
            Report.genome_imbalance += 1
        elif len(self.get_transcripts())!=NUM_TRANSCRIPTS:
            Report.different_transcripts += 1
        elif self.indels_both_genomes():
            Report.indels_both_genomes += 1
        else:
            genome=self.choose_genome_by_indel()
            #print('genome by indel',genome)
            if genome != None:
                Report.total_indels += 1
                r = self.reads
                print("%s %s %s INDEL %s %s %s %s snp %d %d %d %d" %
                    (self.read_name,r[0].gene,r[0].genome,
                    r[2].cigar,r[3].cigar,
                    r[0].cigar,r[1].cigar,
                    r[2].num_mismatch,r[3].num_mismatch,
                    r[0].num_mismatch,r[1].num_mismatch))
            else:
                genome=self.choose_genome_by_snp()
                #print('genome by snp',genome)
                if genome != None:
                    Report.total_snps += 1
                    r = self.reads
                    print("%s %s %s SNP %d %d %d %d indel %s %s %s %s" %
                        (self.read_name,
                        r[0].gene,genome,
                        r[2].num_mismatch,r[3].num_mismatch,
                        r[0].num_mismatch,r[1].num_mismatch,
                        r[2].cigar,r[3].cigar,
                        r[0].cigar,r[1].cigar))

class Streamer:
    """
    Stream the BAM file, process each alignment as it arrives.
    Require that alignment arrive sorted by read ID.
    Expect four alignments per read ID: 2 mates * 2 parental haplotypes.
    """
    def __init__(self, filename):
        self.infile=filename
        self.bamfile = pysam.AlignmentFile(self.infile)

    def stream(self):
        prev_read_name=""
        curr_read_name=""
        read_group=ReadGroup()
        for read in self.bamfile.fetch(until_eof=True):
            #print(read.qname,read.reference_name,read.cigarstring)
            #print(type(read))
            #print(dir(read))
            #print(read.cigarstring)
            curr_read_name=read.qname
            if (curr_read_name != prev_read_name):
                read_group.output()
                read_group.clear()
                read_group.set_name(curr_read_name)
            prev_read_name=curr_read_name
            # Start next read
            read_group.add_read(
                read.reference_name,
                read.cigarstring,
                read.get_tag("NM"),
                read.mapping_quality,
                read.template_length)
        read_group.output()  # last input line is special case

def say(statement):
    #if args.debug:
    print(statement, file=sys.stderr, end="\n")

def args_parse():
    parser = argparse.ArgumentParser(description='Apply IRP rules to mapped pairs.')
    parser.add_argument('bam', help='Alignments file (bam)', type=str)
    #parser.add_argument('genomeA', help = 'Suffix of gene names e.g. Col0', type=str)
    #parser.add_argument('genomeB', help = 'Suffix of gene names e.g. Tsu1', type=str)
    parser.add_argument('--debug', action='store_true')
    return parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    """
    Command line invocation:
    $ python3 Differential_Mapping.py --help
    """
    args = None
    try:
        args = args_parse()
        analyzer = Streamer(args.bam)
        analyzer.stream()
        if args.debug is not None:
            Report.print()
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print('Consider running with --debug')
