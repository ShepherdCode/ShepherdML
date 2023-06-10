import sys
import re
'''
$ samtools view Primary.bam | python bam_to_ml_csv.py > Primary.csv

The bam file must be sorted by read ID,
which is how aligners generate it by default.

File format documentation:
https://samtools.github.io/hts-specs/SAMv1.pdf
https://samtools.github.io/hts-specs/SAMtags.pdf
https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Recent updates:
Show transcript ID.
Show insert lenths.
Show read lengths.
Show high-quality mismatch count.

TO DO:
Encode qualities. Maybe use a bit vec with 256 quals, 0=F, 1=other.
Save cigar or convert it to bit vec of match, mismatch, missing.
Add feature for z-score deviation from mean insert length.

MxM is parent 1. SxS is parent 2
'''
class read_alignment():
    '''Capture one read and its alignment to one transcript.'''

    def __init__(self,score,ed,mm,hqmm,go,ge,hqins,hqdel):
        self.align_score = score
        self.edit_distance = ed
        self.mismatches = mm
        self.high_qual_mismatch = hqmm
        self.gap_opens = go
        self.gap_extends = ge
        self.high_qual_ins = hqins
        self.high_qual_del = hqdel

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
        self.alleles=[None,None,None,None]
        self.alignments=[None,None,None,None]

    def _complete_(self):
        '''Insist on 4 alignments covering 1 transcript, 2 parents.'''
        if None in self.read_lengths: return False
        if None in self.parent_spans: return False
        if None in self.targets: return False
        if None in self.alleles: return False
        if None in self.alignments: return False
        if not (
            self.targets[0]==self.targets[1] and
            self.targets[0]==self.targets[2] and
            self.targets[0]==self.targets[3] and
            self.alleles[0]==self.alleles[1] and
            self.alleles[2]==self.alleles[3] and
            self.alleles[0]!=self.alleles[2] ):
            return False
        return True

    def get_rid(self):
        return self.rid

    def _index_(self,parent,read):
        # 0=P1R1, 1=P1R2, 2=P2R1, 3=P2R2
        return (parent-1)*2 + (read-1)

    def add_alignment(self,parent,read,alignment,tid,allele):
        index = self._index_(parent,read)
        self.alignments[index] = alignment
        self.targets[index] = tid
        self.alleles[index] = allele

    def add_parent_span(self,parent,span):
        self.parent_spans[parent-1] = span

    def add_read_length(self,read,length):
        self.read_lengths[read-1] = length

    def __str__(self):
        return self.rid

    def show_if_complete(self):
        if not self._complete_():
            return False
        # Show the stats for each alignment (2 reads times 2 parents).
        msg =  str(self.alignments[0])+','
        msg += str(self.alignments[1])+','
        msg += str(self.alignments[2])+','
        msg += str(self.alignments[3])+','
        # Show the stats for each read pair (2 pairs).
        msg += str(self.read_lengths[0])+','
        msg += str(self.read_lengths[1])+','
        # Show the stats for each parent (2 parents).
        msg += str(self.parent_spans[0])+','
        msg += str(self.parent_spans[1])+','
        # Show the stats for the whole read group.
        msg += self.targets[0] # assume 4 of the same
        print(msg)
        return True

def parse_cigar(cigar,mdstr,quals):
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

def process_stdin(parent1,parent2):
    '''Online processing of the output of "samtools view alignments.bam"'''
    read_group = None
    printed=0
    incomplete=0
    for line in sys.stdin:
        # samtools view tab-delimited required fields appear in standard order
        fields=line.strip().split('\t')
        RID = fields[0]
        FLAGS = int(fields[1])  # Bitfield of basic info like whether mate aligned.
        if (FLAGS & 0x100):
            PRIMARY = 2
        else:
            PRIMARY = 1
        if FLAGS & 0x40:
            READ = 1
        else:
            READ = 2
        # Parse transcript ID like jg27855.t1_SxS i.e. transcript.suffix_allele
        TID,SUFFIX = fields[2].split('.')
        ALLELE = SUFFIX.split('_')[1]
        if ALLELE==parent1:
            PARENT = 1   # e.g. parent MxM
        else:
            PARENT = 2   # e.g. parent SxS
            if ALLELE != parent2:
                raise Exception('Unrecognized allele: '+ALLELE)
        REF_POS_THIS = int(fields[3]) # Where read alignment starts.
        CIGAR = fields[5] # Cryptic summary of matches, mismatches, indels.
        REF_POS_MATE = int(fields[7])  # Where mate alignment starts.
        # Insert len = aligned pair's end-to-end span along the transcript.
        # Bowtie makes it negative if mate is upstream of this read.
        SPAN = abs(int(fields[8]))
        SEQ = fields[9]
        RLEN = len(SEQ)
        QUALS = fields[10]

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
            elif OPTIONAL_FIELD=='XM:i:':
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
            ALIGN_SCORE,EDIT_DIST,
            MISMATCHES,HQMM,
            GAP_OPENS,GAP_EXTENDS,HQINS,HQDEL)

        # Use read ID to determine whether to start a new read group.
        # Here, rely on assumption that inputs are sorted by read ID.
        if read_group is None:
            # start the first read group
            read_group = pair_alignments(RID)
        elif RID != read_group.get_rid():
            # finish the previous read group
            if read_group.show_if_complete():
                printed += 1
            else:
                incomplete += 1
            # and start a new read group
            read_group = pair_alignments(RID)

        # Accumulate alignments grouped by read ID.
        read_group.add_alignment(PARENT,READ,alignment,TID,ALLELE)
        if READ==1:  # same for both reads, arbitrarily use read 1
            read_group.add_parent_span(PARENT,SPAN)
        if PARENT==1:  # same for both parents, arbitrarily use parent 1
            read_group.add_read_length(READ,RLEN)

        # Feedback while running
        if printed%100000 == 0 or incomplete%100000 == 0 and incomplete>0:
            print('Progress: Printed %d, Incomplete %d'%(printed,incomplete),file=sys.stderr)
    # finish the last read group
    if read_group is not None:
        if read_group.show_if_complete():
            printed += 1
        else:
            incomplete += 1
        read_group = None

if __name__ == '__main__':
    process_stdin('MxM','SxS')
