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

    def __init__(self,read,score,ed,mm,hqmm,go,ge,hqins,hqdel,prim,rlen,span):
        # Most fields are integers
        self.read = read # 1 or 2
        self.align_score = score
        self.edit_distance = ed
        self.mismatches = mm
        self.high_qual_mismatch = hqmm
        self.gap_opens = go
        self.gap_extends = ge
        self.high_qual_ins = hqins
        self.high_qual_del = hqdel
        self.primary = prim  # True or False
        self.rlen = rlen # length of one read in bases
        self.span = span # length of the pair projection on the transcript

    def is_primary(self):
        # Aligners mark some alignments as secondary.
        return prim

    def __str__(self):
        s  = str(self.align_score)+','
        s += str(self.edit_distance)+','
        s += str(self.mismatches)+','
        s += str(self.high_qual_mismatch)+','
        s += str(self.gap_opens)+','
        s += str(self.gap_extends)+','
        s += str(self.high_qual_ins)+','
        s += str(self.high_qual_del)
        return s

    def __repr__(self):
        return str(self)

    def dump_IDs(self,filename='read_ids.csv'):
        with open (filename,'w') as fout:
            for k,v in self.hash.items():
                fout.write('%s,%d\n'%(k,v))

class pair_alignments():
    '''Capture 4 alignments of a read pair.'''
    def __init__(self,rid):
        self.rid=rid  # invariant for 4 reads in group
        self.read_lengths=[None,None] # invariant for 2 parents
        self.parent_spans=[None,None] # invariant for 2 reads
        self.targets=[None,None,None,None]
        self.alignments=[None,None,None,None]

    def is_complete(self):
        '''Insist on 4 alignments covering 1 transcript, 2 parents.'''
        if None in self.read_lengths: return False
        if None in self.parent_spans: return False
        if None in self.targets: return False
        if None in self.alignments: return False
        if not (
            self.targets[0]==self.targets[1] and
            self.targets[0]==self.targets[2] and
            self.targets[0]==self.targets[3] ):
            return False
        return True

    def get_rid(self):
        return self.rid

    def get_preferred_parent(self):
        # Identify the target with the higher combined alignment score.
        if not self.is_complete():
            return 0
        as1 = self.alignments[0].align_score+self.alignments[1].align_score
        as2 = self.alignments[2].align_score+self.alignments[3].align_score
        if as1==as2:
            return 0
        elif as1>as2:
            return 1
        return 2

    def _index_(self,parent,read):
        # Choose list position.
        # 0=P1R1, 1=P1R2, 2=P2R1, 3=P2R2
        return (parent-1)*2 + (read-1)

    def add_alignment(self,parent,read,alignment,tid,allele):
        index = self._index_(parent,read)
        if self.alignments[index] is not None:
            this_AS = alignment.align_score
            prev_align = self.alignments[index]
            prev_AS = prev_align.align_score
            if prev_AS >= this_AS:
                # No thanks. We already have a better one.
                # Bowtie best-2 would never do this.
                # But STAR gives up to 10 aligns per read.
                return False
        self.alignments[index] = alignment
        self.targets[index] = tid
        self.alleles[index] = allele
        return True

    def __str__(self):
        return self.rid

    def show(self):
        if self.is_complete():
            # Show the stats for each alignment (2 reads times 2 parents).
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
            # The primary alignment is to parent 1 or 2 (or 0 undecided).
            msg += str(self.primary)+','
            # By design, the whole read group aligns to one transcript.
            msg += self.targets[0] # assume 4 of the same
            print(msg)
            return True
        return False

class pair_maker():
    '''
    Convert
    from group of alignments (sam read pair, many targets)
    to a pair_alignments (two comparable alignments).
    '''
    def __init__(self,irp=False):
        self.irp = irp

    def parse_cigar(self,cigar,mdstr,quals):
        '''
        Assume a high-quality mismatch or indel is real,
        while a low-quality mismatch or indel is due to sequencing error.
        Use the clue provided by the base call quality score.
        Although scores have 256 values, we'll use a binary representation.
        Most scores are the maximum value, encoded as 'F'.
        So we consider 'F' high quality and everything else low quality.
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
        DIGITS = '0123456789'
        BASES = 'ACGTNABCDEFGHIJKLMNOPQRSTUVWXYZ'
        MAXQ='F'  # SAM encoding of highest quality score
        # parse cigar left to right
        pos = 0
        mquals=''
        while len(cigar)>0:
            match = re.search('^[0-9]+',cigar)
            numstr = match.group(0)
            number = int(numstr)
            cigar = cigar[len(numstr):]
            letter = cigar[0]
            if letter == 'M':  # align but may or may not match
                # Update our qual position and the mquals
                mquals += quals[pos:pos+number]
                pos += number
            elif letter == 'D':  # letters in ref only
                # Don't update our qual position or the mquals
                if quals[pos]==MAXQ and quals[pos+1]==MAXQ:
                    hqdel += number
            elif letter == 'I':
                # Do update our qual position, but don't update the mquals
                for i in range(number):
                    if quals[pos]==MAXQ:
                        hqins += 1
                    pos += 1
            cigar = cigar[1:]
        # parse md string left to right
        pos = 0
        hold = mdstr
        while len(mdstr)>0:
            if mdstr[0]=='^':  # redundant with cigar D, so ignore it
                mdstr = mdstr[1:]
                while mdstr[0] in BASES:
                    mdstr = mdstr[1:]
            elif mdstr[0] in DIGITS: # matched bases, so ignore them
                match = re.search('^[0-9]+',mdstr)
                numstr = match.group(0)
                pos += int(numstr)
                mdstr = mdstr[len(numstr):]
            elif mdstr[0] in BASES: # mismatched bases
                if mquals[pos] == MAXQ:
                    hqmm += 1
                pos += 1
                mdstr = mdstr[1:]
            else:
                raise Exception ('Cannot process mdstr: '+hold)
        return hqmm,hqins,hqdel

    def make_pair(self,group1,group2):
        pa = pair_alignments()
        for record in group1:
            aln = parse_line(record.to_string())
            if aln.is_primary():
                pa.add_alignment(aln)
        for record in group2:
            aln = parse_line(record.to_string())
            if aln.is_primary():
                pa.add_alignment(aln)
        return pa

    def parse_line(self,line):
        '''Process one line of `samtools view alignments.bam` '''
        # samtools view tab-delimited required fields appear in standard order
        fields=line.strip().split('\t')
        # bit=1 means this is a secondary alignment
        RID = fields[0]
        FLAGS = int(fields[1])  # Bitfield of basic info like whether mate aligned.
        PRIMARY = not (FLAGS & 0x100):
        READ = 2
        if FLAGS & 0x40:
            READ = 1
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
        HQMM,HQINS,HQDEL=parse_cigar(CIGAR,MDSTRING,QUALS)
        # save one alignment
        alignment = read_alignment(
            READ,ALIGN_SCORE,EDIT_DIST,
            MISMATCHES,HQMM,
            GAP_OPENS,GAP_EXTENDS,HQINS,HQDEL,
            RLEN,SPAN)
        return alignment

class bam_file():
    '''
    This is an iterator that streams one bam file.
    Each next() returns the alignment group for one read pair.
    '''
    def __init__(self,filename):
        handle = pysam.AlignmentFile(filename, "rb")
        # Fetch puts error message on console, like
        # [E::idx_find_and_load] Could not retrieve index file
        # Pysam gives no way to suppress this.
        print('Error ignored.',file=sys.stderr)
        self.iterator = handle.fetch(until_eof=True)
        self.prev_rec = None
        self.done = False

    def __iter__(self):
        return self

    def __next__(self):
        '''Returns a list of pysam alignment records.'''
        if self.done:
            raise StopIteration()
        group = []
        readname = None
        if self.prev_rec is not None:
            readname = self.prev_rec.to_string().split('\t')[0]
            group.append(self.prev_rec)
        try:
            while True:
                self.prev_rec = next(self.iterator)
                rn = self.prev_rec.to_string().split('\t')[0]
                if readname is None:
                    readname = rn
                if rn==readname:
                    group.append(self.prev_rec)
                else:
                    # We read past this group.
                    # Return this group
                    # but hold onto prev rec.
                    break
        except StopIteration:
            self.prev_rec = None
            self.done = True  # return a group for the last time
        return group

class tandem_file_walker():
    def __init__(self,bamfile1,bamfile2):
        bf1 = bam_file(bamfile1)
        bf2 = bam_file(bamfile2)
        self.iter1 = iter(bf1)
        self.iter2 = iter(bf2)
        self.irp = False

    def set_IRP_mode(self):
        self.irp = True

    def go(self):
        '''Assume records are sorted by read ID.'''
        maker = pair_maker(self.irp)
        try:
            while True:  # till EOF in either file
                group1 = next(self.iter1)
                rn1 = group1[0].to_string().split('\t')[0]
                group2 = next(self.iter2)
                rn2 = group2[0].to_string().split('\t')[0]
                while rn1 != rn2: # advance cursors in tandem
                    if rn1 < rn2:
                        group1 = next(self.iter1)
                        rn1 = group1[0].to_string().split('\t')[0]
                    if rn2 < rn1:
                        group2 = next(self.iter2)
                        rn2 = group1[0].to_string().split('\t')[0]
                #print(rn1,len(group1),rn2,len(group2))
                maker = pm.make_pair(group1,group2)
                if pa.is_complete:
                    pa.show()
        except:
            print('done')

def args_parse():
    parser = argparse.ArgumentParser(description='Parse two BAM files in tandem.')
    parser.add_argument('bam1', help='Reads aligned to parent 1 (BAM)', type=str)
    parser.add_argument('bam2', help='Reads aligned to parent 2 (BAM)', type=str)
    parser.add_argument('--irp', help='Expect allele in transcript ID', action='store_true')
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
    tfw=tandem_file_walker(args.bam1,args.bam2)
    if args.irp:
        tfw.set_IRP_mode()
    tfw.go()
