import argparse
import traceback
import os
import sys
import random

class Reducer ():
    def __init__(self,debug=False):
        self.DEFLINE_PREFIX='>'
        self.TOKEN_SEPARATOR='|'
        self.debug=debug
        self.test_size = 500
        self.challenge_size = 500
        self.valid_size = 500
        self.train_size = 16000
        self.all_genes={} # every gene in input stream
        self.transcript_to_gene={} # map every transcript to gene
        self.set_aside={} # genes removed as train candidates
        self.test_set={} # genes set aside for testing
        self.challenge_set={} # genes set aside for challenge
        self.valid_set={} # genes set aside for validation
        self.train_set={} # transcripts selected for training
        self.infile=None

    def data_in(self,infile):
        self.infile=infile
        with open(infile, 'r') as infa:
            num_seqs = 0
            for line in infa:
                if line[0]==self.DEFLINE_PREFIX:
                    words=line[1:].split(self.TOKEN_SEPARATOR)
                    transcript_id=words[0]
                    gene_id=words[1]
                    self.all_genes[gene_id]=1
                    self.transcript_to_gene[transcript_id]=gene_id

    def reduce(self):
        test_set=[]
        for i in range(self.test_size):
            gene_id=self.get_random_not_taken()
            test_set.append(gene_id)
        self.test_set=test_set
        print("Size of test set is %d"%len(self.test_set))
        challenge_set=[]
        for i in range(self.challenge_size):
            gene_id=self.get_random_not_taken()
            challenge_set.append(gene_id)
        self.challenge_set=challenge_set
        print("Size of challenge set is %d"%len(self.challenge_set))
        valid_set=[]
        for i in range(self.valid_size):
            gene_id=self.get_random_not_taken()
            valid_set.append(gene_id)
        self.valid_set=valid_set
        print("Size of valid set is %d"%len(self.valid_set))

    def make_train(self):
        train_set=[]
        print("Selecting train set...")
        while (len(train_set)<self.train_size):
            tr_list=list(self.transcript_to_gene.keys())
            limit=len(tr_list)
            rnd=random.randrange(0,limit)
            transcript_id=tr_list[rnd]
            gene_id=self.transcript_to_gene[transcript_id]
            self.transcript_to_gene.pop(transcript_id)
            if gene_id not in self.set_aside:
                train_set.append(transcript_id)
        self.train_set=train_set
        print("Size of train set is %d"%len(self.train_set))

    def get_random_not_taken(self):
        list_genes=list(self.all_genes.keys())
        num_genes = len(list_genes)
        while (True):
            rnd=random.randrange(0,num_genes)
            gene_id=list_genes[rnd]
            if gene_id not in self.set_aside:
                self.set_aside[gene_id] = 1
                return gene_id

    def data_out(self,infile):
        self.genes_out(infile,'test',self.test_set)
        self.genes_out(infile,'challenge',self.challenge_set)
        self.genes_out(infile,'valid',self.valid_set)
        self.transcripts_out(infile,'train',self.train_set)

    def transcripts_out(self,infile,name,set):
        FN=name+'.'+self.infile
        with open(FN, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                good_seq=False
                for line in infa:
                    if line[0]==self.DEFLINE_PREFIX:
                        words=line[1:].split(self.TOKEN_SEPARATOR)
                        transcript_id=words[0]
                        if transcript_id in set:
                            good_seq=True
                        else:
                            good_seq=False
                    if good_seq:
                        outfa.write(line)

    def genes_out(self,infile,name,set):
        # Output first transcript for each gene in set.
        # TO DO: output a random transcript for each gene in set.
        FN=name+'.'+self.infile
        with open(FN, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                good_seq=False
                for line in infa:
                    if line[0]==self.DEFLINE_PREFIX:
                        words=line[1:].split(self.TOKEN_SEPARATOR)
                        transcript_id=words[0]
                        gene_id=words[1]
                        if gene_id in set:
                            set.remove(gene_id)
                            good_seq=True
                        else:
                            good_seq=False
                    if good_seq:
                        outfa.write(line)

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''FASTA reduction ala mRNN.'''
    try:
        args_parse()
        fixer = Reducer(args.debug)
        fixer.data_in(args.infile)
        fixer.reduce()
        fixer.make_train()
        fixer.data_out(args.infile)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
