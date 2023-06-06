'''
Given two (homozygote) fasta files,
make one (diploid) fasta file.

Usage:
python make_diploid_fasta.py MxM.fasta SxS.fasta diploid.fasta --name1='MxM' --name2='SxS'
'''

import os
import argparse

class Merger():
    def __init__(self):
        self.name1='MOM'
        self.name2='DAD'
        self.oneline=True
    def set_name1_for_defline(self,name):
        self.name1=name
    def set_name2_for_defline(self,name):
        self.name2=name
    def set_oneline(self,oneline=True):
        '''Whether to remove newlines within a sequence'''
        self.oneline=oneline
    def _fix_defline(self,original,name):
        revised = original
        if original[0]=='>':
            pilons = original.find('_pilon')
            # chop any suffix like _pilon_pilon_pilon
            if pilons > 0:
                revised = original[:pilons]
            # append the name
            revised = revised + '_' + name + '\n'
        return revised
    def _process(self,outstream,instream,name):
        num_deflines = 0
        growing_sequence = ''
        for line in instream:
            if line.startswith('>'):
                num_deflines += 1
                if self.oneline and len(growing_sequence)>0:
                    growing_sequence += '\n'
                    outstream.write(growing_sequence)
                    growing_sequence = ''
                line = self._fix_defline(line,name)
                outstream.write(line)
            elif self.oneline:
                line = line.strip()
                growing_sequence += line
            else:
                outstream.write(line)
        if len(growing_sequence)>0:
            growing_sequence += '\n'
            outstream.write(growing_sequence)
        return num_deflines

    def merge(self,file1,file2,output):
        print('Input file1:',file1)
        print('Input file2:',file2)
        print('Output file:',output)
        with open (output,'w') as fout:
            with open (file1,'r') as fin:
                count1 = self._process(fout,fin,self.name1)
            with open (file2,'r') as fin:
                count2 = self._process(fout,fin,self.name2)
        print(file1,'has',count1,'deflines.')
        print(file2,'has',count2,'deflines.')
        if count1 != count2:
            print('WARNING: the input files differ by number of deflines.')

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='merge fasta files')
    parser.add_argument('input_fasta_1', help = 'filename, read-only', type = str)
    parser.add_argument('input_fasta_2', help = 'filename, read_only', type = str)
    parser.add_argument('--output_fasta', help = 'filename, overwrite (diploid.fasta)',
        type = str, default = 'diploid.fasta')
    parser.add_argument('--name1', help = 'for defline (MOM)',
        type = str, default = 'MOM')
    parser.add_argument('--name2', help = 'for defline (DAD)',
        type = str, default = 'DAD')
    parser.add_argument('--debug', action = 'store_true')

    args = parser.parse_args()

if __name__ == '__main__':
    try:
        args_parse()
        f1 = args.input_fasta_1
        f2 = args.input_fasta_2
        fout = args.output_fasta
        n1 = args.name1
        n2 = args.name2
        merger = Merger()
        if n1 is not None:
            merger.set_name1_for_defline(n1)
        if n2 is not None:
            merger.set_name2_for_defline(n2)
        merger.merge(f1,f2,fout)
    except Exception as e:
        print('ERROR! {}'.format(e))
        if args.debug:
            print(traceback.format_exc())
        else:
            print('Run with --debug for traceback.')
