import argparse
import traceback
import os
import sys
import csv

class FastaLabeler():
    def __init__(self,debug=False):
        self.debug_mode=debug
        self.noncoding='nc' #LncAtlas indicator
        self.measure='CNRCI' # LncAtlas indicator
        self.database={}
    def set_fasta(self,infile,outfile):
        self.say("infile="+infile)
        self.say("outfile="+outfile)
        self.infile=infile
        self.outfile=outfile
    def set_atlas(self,filename):
        self.say("atlas="+filename)
        self.atlas=filename
    def set_celltype(self,celltype):
        self.say("celltype="+celltype)
        self.celltype=celltype
    def parse_atlas(self):
        DB=self.database
        with open(self.atlas,'r') as csvfile:
            reader=csv.DictReader(csvfile)
            for row in reader:
                gene=row['ENSEMBL ID']
                celltype=row['Data Source']
                measure=row['Data Type']
                score=row['Value']
                genetype=row['Biotype']
                if self.is_numeric(score) and \
                genetype==self.noncoding and \
                measure==self.measure and \
                celltype==self.celltype: \
                    DB[gene]=score

    def is_numeric(self,str):
        try:
            val=float(str)
            return True
        except ValueError:
            return False

    def say(self,message):
        if (self.debug_mode):
            print(message)

    def just_do_it(self):
        self.parse_atlas()
        defline=""
        sequence=""
        gene=""
        DB = self.database
        with open(self.infile,'r') as innie:
            with open(self.outfile,'w') as outie:
                for line in innie:
                    if (line[0]=='>'):
                        defline=""
                        rawline=line.rstrip()
                        fields=line.split(' ')
                        geneversion=fields[0][1:]
                        fields=geneversion.split('.')
                        gene=fields[0]
                        if gene in DB.keys():
                            score=DB[gene]
                            defline = rawline + ' ' + score + '\n'
                    else:
                        if len(defline)>0:
                            sequence=line
                            outie.write(defline)
                            outie.write(sequence)

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Label FASTA with LncAtlas score.')
    parser.add_argument(
        'atlasfile', help='LncAtlas filename (csv)', type=str)
    parser.add_argument(
        'celltype', help='LncAtlas ID like HeLa.S3', type=str)
    parser.add_argument(
        'infile', help='input filename (fasta)', type=str)
    parser.add_argument(
        'outfile', help='output filename (fasta)', type=str)
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''Label FASTA with LncAtlas score.'''
    try:
        args_parse()
        labeler = FastaLabeler(args.debug)
        labeler.set_fasta(args.infile,args.outfile)
        labeler.set_atlas(args.atlasfile)
        labeler.set_celltype(args.celltype)
        labeler.just_do_it()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
