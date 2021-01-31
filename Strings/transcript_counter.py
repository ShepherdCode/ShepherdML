import argparse
import traceback
import os
import sys
from datetime import datetime
from math import floor
import random
import numpy as np

class AnnotationSet():
    def __init__(self,debug=False):
        self.debug=debug
        self.gene_to_transcripts={}
        self.transcript_ids={}
        self.lines_ignored=0
        self.comments_ignored=0
    def count_ignored(self):
        return self.lines_ignored
    def add_ignored(self):
        self.lines_ignored += 1
    def count_comments(self):
        return self.comments_ignored
    def add_comment(self):
        self.comments_ignored += 1
    def count_genes(self):
        return len(self.gene_to_transcripts)
    def count_transcripts(self):
        return len(self.transcript_ids)
    def add_transcript(self,gid,tid):
        if tid in self.transcript_ids:
            print("WARNING: duplicate transcript id ",tid)
        else:
            self.transcript_ids[tid]=1
        if gid not in self.gene_to_transcripts:
            self.gene_to_transcripts[gid]=[]
        self.gene_to_transcripts[gid].append(tid)
    def histogram_transcripts_per_gene(self):
        MAX=50
        tpg=np.zeros(MAX)
        for gid in self.gene_to_transcripts:
            cnt = len(self.gene_to_transcripts[gid])
            if cnt>=MAX:
                cnt=MAX-1   # stuff overflow in last bin
            tpg[cnt] += 1
        for i in range(MAX):
            print('%d\t%f'%(i,tpg[i]))

class GFF_Parser():
    '''With GenCode GFF in mind.'''
    def __init__(self,filename,debug=False):
        self.filename=filename
        self.debug=debug
    def load_file(self,type):
        '''Type should be 'nc' or 'pc' or None'''
        COMMENT_PREFIX='#'
        PSEUDOAUTOSOMALREGION='PAR_Y'
        infile=self.filename
        annot = AnnotationSet()
        lno = 0
        fields=None
        line=None
        try:
            with open(infile, 'r') as infa:
                for line in infa:
                    lno += 1
                    if line.startswith(COMMENT_PREFIX):
                        annot.add_comment()
                    else:
                        gff_line=self.parse_line(line)
                        if gff_line['ID'].endswith(PSEUDOAUTOSOMALREGION):
                            annot.add_ignored()
                        elif (type is None or
                        type=='nc' and 'gene_type' in gff_line and gff_line['gene_type']=='lncRNA' or
                        type=='pc' and 'gene_type' in gff_line and gff_line['gene_type']=='protein_coding'):
                            if gff_line['entity']=='transcript':
                                gid=gff_line['gene_id']
                                tid=gff_line['transcript_id']
                                annot.add_transcript(gid,tid)
        except Exception as e:
            print('Problem reading file %s'%infile)
            print('Encountered at line %d'%lno)
            print(line)
            raise(e)
        return annot
    def parse_line(self,oneline):
        FIELD_SEPARATOR='\t'
        gff_line={}
        try:
            fields=oneline.split(FIELD_SEPARATOR)
            gff_line['chr']=fields[0]
            gff_line['entity']=fields[2]
            gff_line['start_pos']=fields[3]
            gff_line['stop_pos']=fields[4]
            gff_line['strand']=fields[6]
            gff_line['extra']=fields[8]
            self.parse_extra(gff_line)
        except Exception as e:
            print('Cannot parse this line: '+oneline)
            print(e)
            raise Exception
        return (gff_line)
    def parse_extra(self,gff_line):
        FIELD_SEPARATOR=';'
        extra=gff_line['extra']
        if extra is not None:
            fields=extra.split(FIELD_SEPARATOR)
            for field in fields:
                if field.startswith('ID='):
                    value=field[3:]
                    gff_line['ID']=value
                elif field.startswith('gene_id='):
                    value=field[8:]
                    gff_line['gene_id']=value
                elif field.startswith('transcript_id='):
                    value=field[14:]
                    gff_line['transcript_id']=value
                elif field.startswith('gene_type='):
                    value=field[10:]
                    gff_line['gene_type']=value

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Preprocess GenCode GFF3 file.')
    parser.add_argument(
        'infile', help='input filename (GFF3)', type=str)
    parser.add_argument(
        '--type', help='pc or nc (all)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        infile=args.infile
        debug=args.debug
        if args.type is None:
            type=None
        else:
            type=args.type.lower()
        parser = GFF_Parser(infile,debug)
        ant = parser.load_file(type)
        print("%d comments ignored"%ant.count_comments())
        print("%d PAR annotations ignored"%ant.count_ignored())
        print("%d genes"%ant.count_genes())
        print("%d transcripts"%ant.count_transcripts())
        ant.histogram_transcripts_per_gene()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
