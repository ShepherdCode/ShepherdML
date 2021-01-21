import argparse
import traceback
import os
import sys
from datetime import datetime
from math import floor
import random

DEFLINE_PREFIX='>'
NEWLINE='\n'
FIELD_SEPARATOR='|'

class Parser():
    def parse_defline(defline):
        try:
            fields=defline.split(FIELD_SEPARATOR)
            transcript=fields[0]
            gene=fields[1]
            length=int(fields[6])
        except Exception as e:
            print('Cannot parse defline '+defline)
            print(e)
            raise Exception
        return (transcript,gene,length)

class Gencode_Preprocess():
    def __init__(self,debug=False):
        self.debug=debug

    def getLen(tuple):
        '''Expect (tid,gid,slen). Return slen.'''
        return tuple[2]

    def getSeq(tuple):
        '''Expect (tid,seq). Return seq.'''
        return tuple[1]



    def process_fasta(self,infile,outfile,good_tid,verbose=True):
        '''Assume FASTA with one line per sequence.'''
        defline=None
        seqline=None
        sn=0  # output sequence number
        with open(outfile, 'w') as outfa:
            with open(infile, 'r') as infa:
                for line in infa:
                    if line[0]==DEFLINE_PREFIX:
                        defline = line
                    else:
                        seqline = line
                        (tid,gid,slen)=Parser.parse_defline(defline)
                        if tid in good_tid:
                            slen=int(slen)
                            # The defline '>' is part of the tid
                            defline="%s %s %d %d\n"%(tid,gid,slen,sn)
                            outfa.write(defline)
                            outfa.write(seqline)
                            sn += 1
        if verbose:
            print("Output_sequences: ",sn)

    def load_transcript_metadata(self,infile):
        '''Load metadata from FASTA file.
        Load data structure of transcript ID & length
        organized by gene ID.'''
        transcripts_per_gene={}
        with open(infile, 'r') as infa:
            for line in infa:
                if line[0]==DEFLINE_PREFIX:
                    (tid,gid,slen)=Parser.parse_defline(line)
                    if not gid in transcripts_per_gene:
                        transcripts_per_gene[gid]=[]
                    tuple=(tid,gid,slen)
                    transcripts_per_gene[gid].append(tuple)
        return transcripts_per_gene

    def subset_transcripts(self,metadata,criteria,verbose=True):
        '''Choose which transcripts IDs to keep.
        Criteria must be one of short|long|median|all.
        Rely on length provided in GenCode defline.'''
        good_tid_dict={}
        num_genes = 0
        num_transcripts = 0
        for transcripts_one_gene in metadata.values():
            num_genes += 1
            num_transcripts += len(transcripts_one_gene)
            if criteria=='all':
                for good_tuple in transcripts_one_gene:
                    good_tid=good_tuple[0] # 0=tid, 1=gid, 2=slen
                    good_tid_dict[good_tid]=1 # 1=exists
            else:
                sorted_tuples=sorted(transcripts_one_gene,key=Gencode_Preprocess.getLen)
                if criteria=='long':
                    good_tuple = sorted_tuples[-1]
                elif criteria=='median':
                    # given 2 median transcripts, we'll take the shorter
                    middle = floor((len(sorted_tuples)-1)/2)
                    good_tuple = sorted_tuples[middle]
                elif criteria=='short':
                    good_tuple = sorted_tuples[0]
                else:
                    sys.exit("Invalid criteria: "+criteria)
                good_tid=good_tuple[0] # 0=tid, 1=gid, 2=slen
                good_tid_dict[good_tid]=1 # 1=exists
        num_good = len(good_tid_dict.keys())
        if verbose:
            print("Infile: ",infile)
            print("Input_Genes: ",num_genes)
            print("Input_Transcripts: ",num_transcripts)
            print("Retained_Transcripts: ",num_good)
        return good_tid_dict

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

    def random_subset(self,ids,max_transcripts):
        random_set={}
        shuffled = list(ids.keys())
        if max_transcripts is None: # subset = entire set
            for tid in shuffled:
                random_set[tid] = 1 # signal existence
        else:
            random.seed(45)
            random.shuffle(shuffled)
            for tid in shuffled[:max_transcripts]:
                random_set[tid] = 1 # signal existence
        return random_set

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Preprocess GenCode FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--keep', help='max transcripts (no limit)', type=int)
    parser.add_argument(
        '--subset', help='short|long|median (all)', type=str, default="all")
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        infile=args.infile
        outfile=args.outfile
        subset=args.subset
        keep=args.keep
        fixer = Gencode_Preprocess(args.debug)
        metadata = fixer.load_transcript_metadata(infile)
        keep_transcripts = fixer.subset_transcripts(metadata,subset)
        keep_transcripts = fixer.remove_duplicates(keep_transcripts,infile)
        keep_transcripts = fixer.random_subset(keep_transcripts,keep)
        fixer.process_fasta(infile,outfile,keep_transcripts)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
