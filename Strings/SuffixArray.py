import csv
import argparse
import traceback
import os
import sys
import csv
from datetime import datetime
from BucketList import BucketList

'''
Suffix Array algorithms. Code by Jason Miller.
This one program operates in either of two modes: create or search.
The mode (create or search) is determined by whether the file exists.
Usage to create a suffix array given a string:
$ rm file.txt; python3 SuffixArray.py no-such-file.txt target-string
Usage to search a file for a given string.
$ python3 SuffixArray.py existing-file.txt query-string
'''

def debug(say1="",say2=""):
    '''Utility: log a message with the time.'''
    if args.debug:
        print (datetime.now(),end=' ')
        print (say1,end=' ')
        print (say2)

class SuffixArray():

    def __init__(self,saFile,inString):
        '''Constructor.'''
        self.saFile = saFile
        self.inString = inString
        self.TERMINATOR = '$'
        self.my_suffix_array = []
        self.stringT = None  # target to search
        self.stringP = None  # pattern or query to search for

    def __str__(self):
        viz = ""
        for s in self.my_suffix_array:
            offset = s - 1
            viz = viz + self.stringT[offset:]
            viz = viz + ", "
        return viz

    def index_at(self,rank):
        '''Utility: use rank to lookup an index stored in the SA.'''
        index = None
        if (rank<len(self.my_suffix_array)):
            index = int(self.my_suffix_array[rank])
        return index

    def suffix_at(self,rank):
        '''Utility: use rank and the SA to get a suffix of T.'''
        suffix = None
        index = self.index_at(rank)
        if (index):
            offset = index - 1  # position 1 -> offset 0
            suffix = self.stringT[offset:]
        return suffix

    def load_suffix_array(self):
        '''Utility: load SA into memory from file.'''
        with open(self.saFile, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            self.stringT = next(reader)[0]
            self.my_suffix_array = next(reader)

    def save_suffix_array(self):
        '''Utility: save in-memory SA to file.'''
        with open(self.saFile, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow([self.stringT])
            writer.writerow(self.my_suffix_array)

    def search_naive(self):
        '''Algorithm: simple Python string search of T looking for P.'''
        # find() returns position of first hit, starting at 0.
        # find() returns -1 for not found.
        # We add one to present one-based counting.
        debug("Start search_naive()")
        position = 1+self.stringT.find(self.stringP)
        debug("Done search_naive()")
        return position

    def search_binary(self):
        '''Algorithm: simple binary search for P in SA.'''
        debug("Start search_binary()")
        NOT_FOUND = 0
        start = 0
        end = len(self.my_suffix_array) # past the end
        rank = 0
        answer = NOT_FOUND
        while (end>start):
            rank = start+(end-start)//2  # depend on decimal truncation
            suffix = self.suffix_at(rank)
            if (suffix.startswith(self.stringP)):
                index = self.index_at(rank)
                answer = index
                end = 0
            elif (end==start+1):
                answer = NOT_FOUND
                end = 0
            elif (self.stringP < suffix):
                end = rank
            else:
                start = rank+1
        debug("Done search_binary()")
        return answer

    def create_manber_myers(self):
        '''Algorithm: build SA by Manber-Myers algorithm.'''
        bucket_list = BucketList(self.stringT)
        bucket_list.string_to_buckets()
        bucket_list.sort_buckets()
        debug("Initial bucket list",bucket_list)
        bucket_list.sort_within_buckets_shortcut()
        debug("Sorted buckets",bucket_list)
        my_suffix_array = []
        for b in bucket_list.get_buckets():
            for s in b.get_suffixes():
                self.my_suffix_array.append(s)
        debug("Suffix array",self.my_suffix_array)
        debug(self)

    def create_naive(self):
        '''Algorithm: build SA just using Python sort and substring.'''
        debug("Start create_naive()")
        L = len(self.stringT)
        one_tuple = ('end',0)  # Just a placeholder.
        self.my_suffix_array = []
        for i in range(0,L):
            one_tuple=(self.stringT[i:L],i+1)
            self.my_suffix_array.append(one_tuple)
        self.my_suffix_array.sort(key=get_first_of_tuple)
        debug("Done create_naive()")

    def main(self,build):
        '''
        In build mode, build SA from string in input file and save SA.
        In search mode, load SA from file and search for input string.
        '''
        if not build:
            debug("Mode = Search Suffix Array.")
            self.stringP = self.inString
            self.load_suffix_array()
            if args.naive:
                index = self.search_naive()
            else:
                index = self.search_binary()
            print(index)
        else:
            debug("Mode = Build Suffix Array.")
            self.stringT = self.string_from_file()
            if (args.naive):
                self.create_naive()
            else:
                self.create_manber_myers()
            self.save_suffix_array()

    def string_from_file(self):
        filename=self.inString
        concat=""
        try:
            with open (filename,"r") as fin:
                for line in fin:
                    trim=line.rstrip("\n\r")
                    concat= concat+trim
                    # This allows many lines per string.
                    # Assume the TERMINATOR was already present.
        except IOError:
            print("Cannot read file "+filename)
        return concat        

def get_first_of_tuple(tup):
    '''Utility.'''
    return tup[0]

def get_second_of_tuple(tup):
    '''Utility.'''
    return tup[1]

def args_parse():
    '''Parse command-line arguments.'''
    global args
    parser = argparse.ArgumentParser(
        description='Suffix Array tools.')
    parser.add_argument(
        'saFile', help='suffix array file to read or create', type=str)
    parser.add_argument(
        'mode', help='build or search', type=str)
    parser.add_argument(
        'inString', help='pattern for search, or file for build', type=str)
    parser.add_argument(
        '--naive', help='Use a simple algorithm.',
        action='store_true')
    parser.add_argument(
        '--debug', help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''Respond to command-line invocation.'''
    try:
        args_parse()
        sa = SuffixArray(args.saFile, args.inString)
        sa.main(args.mode=="build")
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
