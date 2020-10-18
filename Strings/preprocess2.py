import argparse
import traceback
import os
import sys
#import re
from datetime import datetime
#import statistics
from math import floor

DEFLINE_PREFIX='>'
NEWLINE='\n'
FIELD_SEPARATOR='|'

class Sequence():
    def __init__(self):
        self.defline=''
        self.seqline=''

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

class Filter():
    def processing_callback(self,sequence):
        return sequence # override this
    def process_one_sequence(self,sequence):
        sequence=self.processing_callback(sequence)
        return sequence

class Filter_N(Filter):
    def processing_callback(self,sequence):
        if 'N' in sequence.seqline:
            return Sequence() #empty
        return sequence

class All_Caps(Filter):
    def processing_callback(self,sequence):
        sequence.seqline=sequence.seqline.upper()
        return sequence

class Size_Limit(Filter):
    def __init__(self,min,max):
        self.min=min
        self.max=max
        Filter.__init__(self)
    def processing_callback(self,sequence):
        slen = len(sequence.seqline)
        if (slen<self.min or slen>self.max):
            return Sequence() # empty
        return sequence

class Pretty_Defline(Filter):
    def __init__(self):
        self.sequence_counter=0
        Filter.__init__(self)
    def processing_callback(self,sequence):
        defline=sequence.defline
        sn=self.sequence_counter
        if len(defline)>0:
            (transcript,gene,length)=Parser.parse_defline(defline)
            slen=int(length)
            sequence.defline="%s %s len%d seq%d"%(transcript,gene,slen,sn)
            self.sequence_counter += 1
        return sequence

class Filter_By_ID(Filter):
    def __init__(self,keepers):
        self.keepers=keepers
        Filter.__init__(self)
    def processing_callback(self,sequence):
        if len(sequence.defline)>0:
            (tid,gid,slen)=Parser.parse_defline(sequence.defline)
            if tid in self.keepers:
                return sequence
        return Sequence()

class Gencode_Preprocess():
    def __init__(self,debug=False):
        self.debug=debug
        self.filters=[Filter()]

    def add_filter(self,filter):
        self.filters.append(filter)

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
                            defline=">%s %s len%d seq%d"%(tid,gid,slen,sn)
                            outfa.write(defline)
                            outfa.write(seqline)
                            sn += 1
        if verbose:
            print("Output_sequences: ",sn)

    def getLen(tuple):
        '''Expect (tid,gid,slen). Return slen.'''
        return tuple[2]

    def getSeq(tuple):
        '''Expect (tid,seq). Return seq.'''
        return tuple[1]

    def choose_transcript_per_gene(self,infile,criteria='max',verbose=True):
        '''Criteria can be one of: min, max, median.
        Rely on length provided in GenCode defline.'''
        transcripts_per_gene={}
        num_deflines = 0
        num_genes = 0
        num_transcripts = 0
        with open(infile, 'r') as infa:
            for line in infa:
                if line[0]==DEFLINE_PREFIX:
                    num_deflines += 1
                    (tid,gid,slen)=Parser.parse_defline(line)
                    if not gid in transcripts_per_gene:
                        transcripts_per_gene[gid]=[]
                        num_genes += 1
                    tuple=(tid,gid,slen)
                    transcripts_per_gene[gid].append(tuple)
                    num_transcripts += 1
        good_tid_dict={}
        for transcripts_one_gene in transcripts_per_gene.values():
            sorted_tuples=sorted(transcripts_one_gene,key=Gencode_Preprocess.getLen)
            if (criteria=='max'):
                good_tuple = sorted_tuples[-1]
            elif (criteria=='median'):
                # given 2 median transcripts, we'll take the shorter
                middle = floor((len(sorted_tuples)-1)/2)
                good_tuple = sorted_tuples[middle]
            else: # (criteria=='min'):
                good_tuple = sorted_tuples[0]
            good_tid=good_tuple[0] # 0=tid, 1=gid, 2=slen
            good_tid_dict[good_tid]=1 # 1=exists
        num_good = len(good_tid_dict.keys())
        if verbose:
            print("Infile: ",infile)
            print("Input_Deflines: ",num_deflines)
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

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Fix a FASTA file.')
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''Preprocess FASTA from GenCode.'''
    try:
        args_parse()
        infile=args.infile
        outfile=args.outfile
        fixer = Gencode_Preprocess(args.debug)
        keep_transcripts = fixer.choose_transcript_per_gene(infile,'median')
        keep_transcripts = fixer.remove_duplicates(keep_transcripts,infile)
        #fixer.add_filter(All_Caps())
        #fixer.add_filter(Filter_N())
        #fixer.add_filter(Filter_By_ID(keep_transcripts))
        #fixer.add_filter(Size_Limit(200,30000))
        #fixer.add_filter(Pretty_Defline()) # must be last
        fixer.process_fasta(infile,outfile,keep_transcripts)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
