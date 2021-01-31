import argparse
import traceback
import os
import sys
from datetime import datetime
from math import floor
import random

class AnnotationSet():
    def __init__(self,debug=False):
        self.debug=debug
        self.transcripts_per_gene={}
        self.exons_per_transcript={}
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
        return len(self.transcripts_per_gene)
    def add_gene(self,id):
        if id in self.transcripts_per_gene:
            print("Warning: id exists ",id)
        self.transcripts_per_gene[id] = 0
    def count_transcripts(self):
        return len(self.exons_per_transcript)
    def add_transcript(self,id):
        self.exons_per_transcript[id] = 0

class GFF_Parser():
    '''With GenCode GFF in mind.'''
    def __init__(self,filename,debug=False):
        self.filename=filename
        self.debug=debug
    def load_file(self):
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
                        elif gff_line['entity']=='gene':
                            id=gff_line['gene_id']
                            annot.add_gene(id)
                        elif gff_line['entity']=='transcript':
                            id=gff_line['transcript_id']
                            annot.add_transcript(id)
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

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Preprocess GenCode GFF3 file.')
    parser.add_argument(
        'infile', help='input filename (GFF3)', type=str)
    #parser.add_argument(
    #    '--keep', help='num transcripts (no limit)', type=int)
    #parser.add_argument(
    #    '--subset', help='short|long|median (all)', type=str, default="all")
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        infile=args.infile
        debug=args.debug
        parser = GFF_Parser(infile,debug)
        ant = parser.load_file()
        print("%d comments ignored"%ant.count_comments())
        print("%d PAR annotations ignored"%ant.count_ignored())
        print("%d genes"%ant.count_genes())
        print("%d transcripts"%ant.count_transcripts())
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
