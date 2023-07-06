import sys
import re
import argparse
import pysam
'''
Align reads to the M parent and again to the S parent.
$ samtools sort -n -T . -@ 4 -o Sorted.M.bam Primary.M.bam
$ samtools sort -n -T . -@ 4 -o Sorted.S.bam Primary.S.bam
$ python bam_two_targets.py > ml_stats.csv

Although pysam can return the group of alignments for a read,
that requires a bam index and random access.
We prefer to sort both bam files and stream them in tandem.

For testing, use a reproducible 0.1% subset (seed 123).
$ samtools view -bo subset.M.bam -s 123.001 Sorted.M.bam
'''
class read_alignment():
    '''Capture one read and its alignment to one transcript.'''

    def __init__(self,rid,tid,read_num,
    score,ed,mm,hqmm,go,ge,hqins,hqdel,ins,dels,mat,prim,rlen,span):
        self.rid = rid # long string
        self.tid = tid # string containing gene or transcript name
        # BAM files contain a flag for primary or secondary alignment
        self.primary = prim  # True or False
        # these fields are integers
        self.read_num = read_num # 1 or 2
        self.align_score = score
        self.edit_distance = ed
        self.mismatches = mm
        self.high_qual_mismatch = hqmm
        self.gap_opens = go
        self.gap_extends = ge
        self.high_qual_ins = hqins
        self.high_qual_del = hqdel
        self.insertions = ins
        self.deletions = dels
        self.matches = mat
        self.rlen = rlen # length of one read in bases
        self.span = span # length of the pair projection on the transcript

    def __str__(self):
        # This method supports the main reporting feature of this program.
        s  = str(self.align_score)+','
        s += str(self.edit_distance)+','
        s += str(self.matches)+','
        s += str(self.mismatches)+','
        s += str(self.high_qual_mismatch)+','
        s += str(self.gap_opens)+','
        s += str(self.gap_extends)+','
        s += str(self.insertions)+','
        s += str(self.deletions)+','
        s += str(self.high_qual_ins)+','
        s += str(self.high_qual_del)
        return s

    def __repr__(self):
        return str(self)

    def dump_IDs(self,filename='read_ids.csv'):
        with open (filename,'w') as fout:
            for k,v in self.hash.items():
                fout.write('%s,%d\n'%(k,v))

class four_alignments():
    '''
    Capture 4 best alignments of 1 read pair to 1 target in 2 parents.
    Input groups are lists of instances of read_alignment.
    Both groups should concern the same readname ie read pair.
    Both groups should concern the same target in different parents.
    Together, the two groups should be exhaustive for the readname.
    '''
    def __init__(self,group1,group2):
        self.alignments=[None,None,None,None] # order P1R1,P1R2,P2R1,P2R2
        for record in group1:
            self.add_alignment(1, record) # keeps the 4 best
        for record in group2:
            self.add_alignment(2, record) # keep the 4 best
        self.different_targets = False

    def allow_different_targets(self):
        self.different_targets = True

    def add_alignment(self,parent,alignment):
        '''
        Retain the best by alignment score.
        Retain at most 4 alignments: P1R1,P1R2,P2R1,P2R2.
        Insist that all alignments concern same read.
        '''
        if not alignment.primary:
            # No thanks. We don't want secondary alignments.
            return False
        # Choose list position.
        # 0=P1R1, 1=P1R2, 2=P2R1, 3=P2R2
        read_num = alignment.read_num
        index = (parent-1)*2 + (read_num-1)
        if self.alignments[index] is not None:
            this_AS = alignment.align_score
            prev_AS = self.alignments[index].align_score
            if prev_AS >= this_AS:
                # No thanks. We already have a better one.
                # This should never happen.
                # Configure each aligner to chose ONE primary target,
                # so we get consistency within a read pair.
                # Force the aligner to break ties.
                # Do not use the STAR option to mark ties all primary.
                return False
        # Ok, store this alignment as one-of-four.
        self.alignments[index] = alignment
        return True

    def is_complete(self):
        if None in self.alignments:
            # Insist on 4 alignments.
            # Implicitly, this requires 2 aligns to 2 parents.
            return False
        if not self.different_targets:
            # Insist on 4 alignments to same gene (in different genomes)
            same_tid=self.alignments[0].tid
            if self.alignments[1].tid != same_tid or \
                self.alignments[2].tid != same_tid or \
                self.alignments[3].tid != same_tid:
                return False
        return True

    def get_preferred_parent(self):
        # Identify the parent with the higher combined alignment score.
        if not self.is_complete():
            return 0
        sum_score_1 = self.alignments[0].align_score+self.alignments[1].align_score
        sum_score_2 = self.alignments[2].align_score+self.alignments[3].align_score
        if sum_score_1==sum_score_2:
            return 0
        elif sum_score_1>sum_score_2:
            return 1
        return 2  # sum_score_1<sum_score_2

    def show(self):
        # This method implements the main reporting feature of this program.
        if self.is_complete():
            # Stats for each alignment (2 reads times 2 parents).
            msg =  str(self.alignments[0])+','
            msg += str(self.alignments[1])+','
            msg += str(self.alignments[2])+','
            msg += str(self.alignments[3])+','
            # Assume same for first read pair and second read pair.
            msg += str(self.alignments[0].rlen)+','
            msg += str(self.alignments[1].rlen)+','
            # Assume same span given for both reads of a pair.
            msg += str(self.alignments[0].span)+','
            msg += str(self.alignments[2].span)+','
            # Preferred is parent 1 or 2 (or 0 undecided).
            msg += str(self.get_preferred_parent())+','
            # This feature is only meaningful if
            # self.require_same_target==True
            msg += self.alignments[0].tid

            print(msg) # TO DO: redirect to file
            return True  # enable caller to count outputs
        else:
            # For debugging, explain why this is incomplete
            #print(self.alignments[0])
            #print(self.alignments[1])
            #print(self.alignments[2])
            #print(self.alignments[3])
            #raise Exception
            pass
        return False

class bam_line_parser():
    def __init__(self,irp,maxq):
        self.irp = irp  # boolean
        self.MAXQ=maxq  # SAM ASCII encoding of highest quality score ('K').

    def parse_cigar(self,cigar,mdstr,quals):
        '''
        Assume a high-quality mismatch or indel is real,
        while a low-quality mismatch or indel is due to sequencing error.
        Use the clue provided by the base call quality score.
        Although scores have 256 values, we'll use a binary representation.
        Most scores by far are the maximum value, encoded as 'K'.
        So we consider that high quality and everything else low quality.
        The SAM/BAM format makes quality extraction very difficult.
        Only the cigar captures where the read indels are.
        Example cigar: 1M2I3M4D5M means
        1 align, 2 extra read bases, 3 align, 4 extra ref bases, 5 align).
        Only the MD str captures where the mismatches are.
        Example MD str: 5^GA2C5 means
        5 match, GA in ref only, 2 match, C in ref is mismatch, 5 match.
        '''
        hqmm = 0
        hqins = 0
        hqdel = 0
        insertions = 0
        deletions = 0
        matches = 0
        DIGITS = '0123456789'
        # ACGT indicates reference letter, different from read.
        # M aligns, could be match or mismatch
        # I read has extra bases not in reference
        # D reference has extra bases not in read
        # N reference has intron not in read
        # S soft clip a read end - skip these
        # H hard clip - read as shown was trimmed before alignment
        # P padding - the reference bases are just padding
        BASES = 'ACGTABCDEFGHIJKLMNOPQRSTUVWXYZ'
        MAXQ=self.MAXQ  # SAM ASCII encoding of highest quality score.
        # Start at position 0 and parse cigar left to right.
        pos = 0
        # Modify the quality string if cigar specifies indels.
        mquals=''  # accumulate and save this for the MD string parse
        while len(cigar)>0:
            match = re.search('^[0-9]+',cigar)
            numstr = match.group(0)
            number = int(numstr)
            cigar = cigar[len(numstr):]
            letter = cigar[0]
            if letter == 'M' or letter == 'S':  # align but may or may not match
                # Update our qual position and the mquals
                mquals += quals[pos:pos+number]
                pos += number
            elif letter == 'D':  # letters in ref only
                # Don't update our qual position or the mquals
                deletions += number
                if quals[pos]==MAXQ and quals[pos+1]==MAXQ:
                    # Require high quality calls on both sides of deletion
                    hqdel += number
            elif letter == 'I':   # letters in read only
                # Do update our qual position, but don't update the mquals
                insertions += number
                for i in range(number):
                    if quals[pos]==MAXQ:
                        # Require high quality at inserted base
                        hqins += 1
                    pos += 1
            elif letter == 'N': # intron
                pass
            cigar = cigar[1:]
        # parse md string left to right
        pos = 0
        hold = mdstr
        while len(mdstr)>0:
            if mdstr[0]=='^':  # redundant with cigar D, so ignore it
                mdstr = mdstr[1:]
                while mdstr[0] in BASES:
                    mdstr = mdstr[1:]
            elif mdstr[0] in DIGITS: # matched bases
                match = re.search('^[0-9]+',mdstr)
                numstr = match.group(0)
                number = int(numstr)
                pos += number
                matches += number
                mdstr = mdstr[len(numstr):]
            elif mdstr[0] in BASES: # mismatched bases
                # Redundant with MM field, so mostly ignore this
                # This counts even soft clipped bases, which bowtie AS ignores
                if mquals[pos] == MAXQ:
                    # Require high quality score at mismatched base
                    hqmm += 1
                pos += 1
                mdstr = mdstr[1:]
            else:
                raise Exception ('Cannot process mdstr: '+hold)
        return hqmm,hqins,hqdel,insertions,deletions,matches

    def parse_line(self,fields):
        '''
        Process one record of a BAM file.
        Assume record has been converted to string and split on tabs,
        similar to samtools view | split('\t').
        Return an instance of read_alignment.
        '''
        # bit=1 means this is a secondary alignment
        RID = fields[0]
        FLAGS = int(fields[1])  # Bitfield of basic info like whether mate aligned.
        # We don't need any secondary alignments.
        PRIMARY = not (FLAGS & 0x100)
        READ_NUM = 2
        if FLAGS & 0x40:
            READ_NUM = 1
        TID = fields[2]
        if self.irp:
            # Fasta files from the IRP1 pipeline have parent suffixes like tid100_MxM.
            # We strip the suffix to use same tid in both parents.
            TID = TID.split('_')[0]
        REF_POS_THIS = int(fields[3]) # Where read alignment starts.
        CIGAR = fields[5] # Cryptic summary of matches, mismatches, indels.
        REF_POS_MATE = int(fields[7])  # Where mate alignment starts.
        # Insert len = aligned pair's end-to-end span along the transcript.
        # Bowtie makes it negative if mate is upstream of this read.
        SPAN = abs(int(fields[8]))
        SEQ = fields[9]
        RLEN = len(SEQ)
        QUALS = fields[10]
        # Loop through SAM's optional fields.
        # (We require certain ones even though SAM does not.)
        # samtools view tab-delimited optional fields appear in any order
        ALIGN_SCORE = -1000 # In case none provided, assume a very negative score.
        MISMATCHES = 0
        GAP_OPENS = 0
        GAP_EXTENDS = 0
        EDIT_DIST = 0
        MDSTRING = ''
        for j in range(11,len(fields)):
            OPTIONAL_FIELD=fields[j][:5]
            OPTIONAL_VALUE=fields[j][5:]
            if OPTIONAL_FIELD=='AS:i:':
                ALIGN_SCORE = int(OPTIONAL_VALUE)
            elif OPTIONAL_FIELD=='XM:i:' or OPTIONAL_FIELD=='nM:i:' :
                # Bowtie outputs XM, STAR outputs nM
                MISMATCHES = int(OPTIONAL_VALUE)
            elif OPTIONAL_FIELD=='XO:i:':
                GAP_OPENS = int(OPTIONAL_VALUE)
            elif OPTIONAL_FIELD=='XG:i:':
                GAP_EXTENDS = int(OPTIONAL_VALUE)
            elif OPTIONAL_FIELD=='NM:i:':
                EDIT_DIST = int(OPTIONAL_VALUE)
            elif OPTIONAL_FIELD=='MD:Z:':
                MDSTRING = OPTIONAL_VALUE
        # compute hiqh quality mismatches
        HQMM,HQINS,HQDEL,INS,DELS,MAT=self.parse_cigar(CIGAR,MDSTRING,QUALS)
        # save one alignment
        alignment = read_alignment(
            RID,TID,READ_NUM,ALIGN_SCORE,EDIT_DIST,
            MISMATCHES,HQMM,GAP_OPENS,GAP_EXTENDS,HQINS,HQDEL,
            INS,DELS,MAT,
            PRIMARY,RLEN,SPAN)
        return alignment

class bam_file_parser():
    '''
    This is an iterator that streams one bam file.
    Each next() returns a list of alignments for one readname.
    Each alignment is an instance of read_alignment.
    '''
    def __init__(self,filename,irp,maxq):
        handle = pysam.AlignmentFile(filename, "rb")
        # Fetch puts error message on console, like
        # [E::idx_find_and_load] Could not retrieve index file
        # Pysam gives no way to suppress this.
        print('Error ignored.',file=sys.stderr)
        self.iterator = handle.fetch(until_eof=True)
        self.prev_aln = None
        self.done = False
        self.line_parser = bam_line_parser(irp,maxq)

    def __iter__(self):
        return self

    def __next__(self):
        if self.done:
            raise StopIteration()
        group = []  # list of alignments
        readname = None
        if self.prev_aln is not None:
            # We previously streamed first record of next group.
            # Now, that record starts this new group.
            readname = self.prev_aln.rid
            group.append(self.prev_aln)
        try:
            while True:
                pysam_record = next(self.iterator)
                # Pysam provides a python-style data structure.
                # It differs from 'samtools view' in several ways,
                # including 0-based positions and quality scores.
                # For equivalence with 'samtools view',
                # we convert the pysam record to list of string
                # and parse those strings.
                fields = pysam_record.to_string().split('\t')
                # Convert strings to instance of read_alignment.
                self.prev_aln = self.line_parser.parse_line(fields)
                rn = self.prev_aln.rid
                if readname is None:
                    # Special case: first record in the file
                    readname = rn
                if rn==readname:
                    # Add this alignment to our group and continue streaming.
                    group.append(self.prev_aln)
                else:
                    # We streamed too far. Return the group that we have.
                    # Hold onto prev_rec for next time.
                    break
        except StopIteration:
            self.prev_aln = None
            self.done = True  # return a group for the last time
        return group

class tandem_file_walker():
    '''
    Stream two BAM files in tandem.
    Apply a bam_file_parser to each file.
    Collect all the BAM records from both files pertaining to one readname.
    That is, all the alignments for one read pair.
    Pass each collection to four_alignments for processing.
    '''
    def __init__(self,bamfile1,bamfile2,irp,diff,maxq):
        bf1 = bam_file_parser(bamfile1,irp,maxq)
        bf2 = bam_file_parser(bamfile2,irp,maxq)
        self.parser1 = iter(bf1)
        self.parser2 = iter(bf2)
        self.irp = irp
        self.diff = diff

    def go(self):
        '''
        Extract records from two BAM files, one for each parent.
        Assume BAM records are sorted by readname.
        Maintain a curcor in each BAM file.
        Advance both cursors to the next readname common to both files.
        Load and process all the lines containing that readname.
        '''
        total = 0
        try:
            while True:  # loop till EOF in either file
                group1 = next(self.parser1)
                rn1 = group1[0].rid
                group2 = next(self.parser2)
                rn2 = group2[0].rid
                while rn1 != rn2: # advance cursors in tandem
                    if rn1 < rn2:
                        group1 = next(self.parser1)
                        rn1 = group1[0].rid
                    if rn2 < rn1:
                        group2 = next(self.parser2)
                        rn2 = group2[0].rid
                pair = four_alignments(group1,group2)
                if self.diff:
                    pair.allow_different_targets()
                if pair.show():
                    total += 1
        except StopIteration:
            # No need to special case the last read because
            # our parsers raise the exception AFTER last read.
            print('Total %d'%total, file=sys.stderr)

def args_parse():
    parser = argparse.ArgumentParser(description='Parse two BAM files in tandem.')
    parser.add_argument('bam1', help='Reads aligned to parent 1 (BAM)', type=str)
    parser.add_argument('bam2', help='Reads aligned to parent 2 (BAM)', type=str)
    parser.add_argument('--irp', help='Expect allele in transcript ID', action='store_true')
    parser.add_argument('--diff', help='Allow different target per parent', action='store_true')
    parser.add_argument('--maxq', help='Max quality encoding (K)', default='K')
    parser.add_argument('--debug', action='store_true')
    return parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    args = None
    try:
        args = args_parse()
    except Exception as e:
        print("\nError parsing command line.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print('Consider running with --debug')
    tfw=tandem_file_walker(args.bam1,args.bam2,
        args.irp,args.diff,args.maxq)
    tfw.go()
