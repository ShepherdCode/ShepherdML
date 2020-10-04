import argparse
import traceback
import os
import sys
import random

class Reducer ():
    def __init__(self,debug=False):
        self.DEFLINE_PREFIX='>'
        self.debug=debug
        self.test_size = 500
        self.challenge_size = 500
        self.valid_size = 500
        self.train_size = 16000
        self.all_genes={}
        self.transcript_to_gene={}
        self.set_aside={}
        self.train_set_transcripts={}

    def data_in(self,infile):
        TOKEN_SEPARATOR='|'
        with open(infile, 'r') as infa:
            num_seqs = 0
            for line in infa:
                if line[0]==self.DEFLINE_PREFIX:
                    words=line.split(TOKEN_SEPARATOR)
                    transcript_id=words[0]
                    gene_id=words[1]
                    self.all_genes[gene_id]=1
                    self.transcript_to_gene[transcript_id]=gene_id

    def reduce(self):
        test_set=[]
        for i in range(self.test_size):
            gene_id=self.get_random_not_taken()
            test_set.append(gene_id)
        print("Size of test set is %d"%len(test_set))
        challenge_set=[]
        for i in range(self.challenge_size):
            gene_id=self.get_random_not_taken()
            challenge_set.append(gene_id)
        print("Size of challenge set is %d"%len(challenge_set))
        valid_set=[]
        for i in range(self.valid_size):
            gene_id=self.get_random_not_taken()
            valid_set.append(gene_id)
        print("Size of valid set is %d"%len(valid_set))

    def make_train(self):
        train_set=[]
        while (len(train_set)<self.train_size):
            tr_list=list(self.transcript_to_gene.keys())
            limit=len(tr_list)
            rnd=random.randrange(0,limit)
            transcript_id=tr_list[rnd]
            gene_id=self.transcript_to_gene[transcript_id]
            self.transcript_to_gene.pop(transcript_id)
            if gene_id not in self.set_aside:
                train_set.append(transcript_id)
        self.train_set_transcripts=train_set

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
        #with open(self.outfile, 'w') as outfa:
        pass

    def print_prev(self,outfile,defline,seq):
        NL='\n'
        allcaps=seq.upper()
        if self.delete_N and 'N' in allcaps:
            pass
        else:
            outfile.write(defline+NL)
            if self.allcaps:
                outfile.write(allcaps)
            else:
                outfile.write(seq)
            outfile.write(NL)

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
