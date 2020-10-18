import argparse
import traceback
import os
import sys
from datetime import datetime
from math import floor
import random

DEFLINE_PREFIX='>'
NEWLINE='\n'
FIELD_SEPARATOR=' '

class Parser():
    def parse_defline(defline):
        try:
            fields=defline.split(FIELD_SEPARATOR)
            transcript=fields[0]
            gene=fields[1]
            length=int(fields[2])
        except Exception as e:
            print('Cannot parse defline '+defline)
            print(e)
            raise Exception
        return (transcript,gene,length)

class Gencode_Splitter():
    def __init__(self,debug,infile):
        self.debug=debug
        self.infile=infile
        self.id_to_len={}
        self.outfile=None

    def set_output_file(self,outfile):
        self.outfile=outfile

    def process_fasta(self,good_tid,verbose=True):
        '''Assume FASTA with one line per sequence.'''
        defline=None
        seqline=None
        sn=0  # output sequence number
        with open(self.outfile, 'w') as outfa:
            with open(self.infile, 'r') as infa:
                for line in infa:
                    if line[0]==DEFLINE_PREFIX:
                        defline = line
                    else:
                        seqline = line
                        (tid,gid,slen)=Parser.parse_defline(defline)
                        if tid in good_tid:
                            slen=int(slen)
                            # The defline '>' is part of the tid
                            defline="%s %s len=%d seq=%d\n"%(tid,gid,slen,sn)
                            outfa.write(defline)
                            outfa.write(seqline)
                            sn += 1
        if verbose:
            print("Output_filename: ",self.outfile)
            print("Output_sequences: ",sn)

    def remove_ids(self,bad_ids,verbose=True):
        plen=len(self.id_to_len.keys())
        dlen=len(bad_ids.keys())
        for id in bad_ids:
            del self.id_to_len[id]
        nlen=len(self.id_to_len.keys())
        if verbose:
            print("Num_ids_before_delete: ",plen)
            print("Num_ids_to_delete: ",dlen)
            print("Num_ids_after_delete: ",nlen)

    def load_ids(self,verbose=True):
        self.id_to_len={}
        num_deflines = 0
        num_genes = 0
        num_transcripts = 0
        num_too_short = 0
        with open(self.infile, 'r') as infa:
            for line in infa:
                if line[0]==DEFLINE_PREFIX:
                    num_deflines += 1
                    (tid,gid,slen)=Parser.parse_defline(line)
                    self.id_to_len[tid]=slen
                    num_transcripts += 1
        if verbose:
            print("Infile: ",self.infile)
            print("Input_Deflines: ",num_deflines)
            print("Retained_Transcripts: ",num_transcripts)

    def make_test_set(self,tsize,verbose=True):
        REPEATABLE = 45
        shuffled=list(self.id_to_len.keys())
        random.seed(REPEATABLE)
        random.shuffle(shuffled)
        test_list = shuffled[:tsize]
        test_set={}
        for id in test_list:
            test_set[id]=1 # signal existence
        if verbose:
            print("Test set size: ",len(test_set))
        return test_set

    def make_train_set(self,maxlen,tsize,verbose=True):
        REPEATABLE = 55
        shuffled=list(self.id_to_len.keys())
        random.seed(REPEATABLE)
        random.shuffle(shuffled)
        train_list = shuffled[:tsize]
        train_set={}
        selected=0
        for id in train_list:
            if selected < tsize:
                if self.id_to_len[id] <= maxlen:
                    selected += 1
                    train_set[id]=1 # signal existence
        if verbose:
            print("Train set size: ",len(train_set))
        return train_set

    def remove_duplicates(self,good_tid,infile,verbose=True):
        '''GenCode includes ~60 identical sequences from different genes.
        Sort by len, then test neighbors for seq identity.'''
        tuples=[]
        with open(infile, 'r') as infa:
            (tid,gid,slen) = ('','','')
            seq=''
            # Iterate through multi-line FASTA file.
            # TO DO: use a sequence iterator.
            for line in infa:
                if line[0]==DEFLINE_PREFIX:
                    if len(seq)>0:
                        tuple=(tid,seq)
                        tuples.append(tuple)
                    (tid,gid,slen)=Parser.parse_defline(line)
                    seq=''
                elif tid in good_tid:
                    seq=seq+line.rstrip()
            if len(seq)>0:
                tuple=(tid,seq) # tuple = id, sequence
                tuples.append(tuple)
        sorted_tuples=sorted(tuples,key=Gencode_Preprocess.getSeq)
        reduced_set = {}
        for i in range(1,len(sorted_tuples)):
            this_tid=sorted_tuples[i][0]
            prev_seq=sorted_tuples[i-1][1]
            this_seq=sorted_tuples[i][1]
            if prev_seq!=this_seq:
                reduced_set[this_tid]=1
        if verbose:
            print("Transcripts_After_Remove_Dupes: ",len(reduced_set.keys()))
        return reduced_set

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Split FASTA into test and train sets.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output tag (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        splitter = Gencode_Splitter(args.debug,args.infile)
        splitter.load_ids()
        test_ids = splitter.make_test_set(1000)
        testfile = 'test.'+args.outfile
        splitter.set_output_file(testfile)
        splitter.process_fasta(test_ids)
        trainfile = 'train.'+args.outfile
        splitter.remove_ids(test_ids)
        train_ids=splitter.make_train_set(1000,16000)
        splitter.set_output_file(trainfile)
        splitter.process_fasta(train_ids)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
