import sys
'''
$ samtools view Primary.bam | python bam_to_ml_csv.py > Primary.csv

See https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

Updated Jun 6 to show transcript ID and insert lenths.
TO DO: Encode qualities. Maybe use a bit vec with 256 quals, 0=F, 1=other.
TO DO: Save cigar or convert it to bit vec of match, mismatch, missing.
TO DO: add feature for deviation from mean insert length

MxM is parent 1. SxS is parent 2
'''
class read_alignment():
    '''Capture one read and its alignment to one transcript.'''
    def __init__(self,ip,span,score,ed,mm,go,ge,tid):
        self.is_primary = ip
        self.span = span
        self.align_score = score
        self.edit_distance = ed
        self.mismatches = mm
        self.gap_opens = go
        self.gap_extends = ge
        self.transcript = tid
        
    def __str__(self):
        s = 'S,'
        if self.is_primary:
            s = 'P,'
        s += str(self.align_score)+','
        s += str(self.edit_distance)+','
        s += str(self.mismatches)+','
        s += str(self.gap_opens)+','
        s += str(self.gap_extends)
        return s
    
    def __repr__(self):
        return str(self)
    
class read_IDs():
    '''Convert read ID from very long string to int.
       So far, no use for this class.'''
    def __init__(self):
        self.hash =  dict()
        next_iid = 1
     
    def convert_id(self,rid):
        if rid in self.hash:
            iid = self.hash[rid]
        else:
            iid = next_iid
            self.hash[rid] = iid
            next_iid += 1
        return iid
    
    def dump_IDs(self,filename='read_ids.csv'):
        with open (filename,'w') as fout:
            for k,v in self.hash.items():
                fout.write('%s,%d\n'%(k,v))
    
class pair_alignments():
    '''Capture 4 alignments: P1 R1, P1 R2, P2 R1, P2 R2.
    Two reads per pair, each aligned to same transcript in two parents.''' 
    def __init__(self,rid):
        self.rid=rid
        self.targets=set()
        self.alleles=set()
        self.alignments=[None,None,None,None]
        
    def get_rid(self):
        return self.rid
    
    def has_four(self):
        missing = (None in self.alignments)
        return not missing
    
    def _add_(self,alignment,index,tid,allele):
        self.alignments[index] = alignment
        self.targets.add(tid)
        self.alleles.add(allele)
        
    def add_P1_R1(self,alignment,tid,allele):
        self._add_(alignment,0,tid,allele)
    def add_P1_R2(self,alignment,tid,allele):
        self._add_(alignment,1,tid,allele)
    def add_P2_R1(self,alignment,tid,allele):
        self._add_(alignment,2,tid,allele)
    def add_P2_R2(self,alignment,tid,allele):
        self._add_(alignment,3,tid,allele)
        
    def __str__(self):
        return self.rid
    
    def show_if_good(self):
        good = self.has_four() and \
            len(self.targets)==1 and \
            len(self.alleles)==2
        if good:
            # Show the basic alignments stats for 2 reads times 2 parents.
            msg =  str(self.alignments[0])+','
            msg += str(self.alignments[1])+','
            msg += str(self.alignments[2])+','
            msg += str(self.alignments[3])+','
            # Get the span of the pair alignments to either parent.
            # In bowtie output, reads 0 and 1 have same span, opposite sign.
            # In bowtie output, reads 2 and 3 have same span, opposite sign.
            span_p1 = self.alignments[0].span 
            span_p2 = self.alignments[2].span 
            msg += str(span_p1)+','+str(span_p2)+','
            # Because we filtered for four-reads-to-one-gene,
            # we only need report the first gene (actually transcript) ID.
            msg += self.alignments[0].transcript
            print(msg)

def process_stdin(parent1,parent2):
    '''Online processing of the output of "samtools view alignments.bam"'''
    read_group = None
    for line in sys.stdin:

        # samtools view tab-delimited required fields appear in standard order
        fields=line.strip().split('\t')
        RID = fields[0]
        FLAGS = int(fields[1])  # Bitfield of basic info like whether mate aligned.
        IS_PRIMARY = not (FLAGS & 0x100) # This bit is 1 on secondary alignments.
        is_read1 = FLAGS & 0x40  # This bit is 1 on first read of a pair.
        # Parse transcript ID like jg27855.t1_SxS i.e. transcript.suffix_allele
        TID,SUFFIX = fields[2].split('.') 
        ALLELE = SUFFIX.split('_')[1]
        if ALLELE==parent1: 
            is_parent1 = True   # e.g. parent MxM
        else:
            is_parent1 = False   # e.g. parent SxS
            if ALLELE != parent2:
                raise Exception('Unrecognized allele: '+ALLELE)
        REF_POS_THIS = int(fields[3]) # Where read alignment starts.
        CIGAR = fields[5] # Cryptic summary of matches, mismatches, indels.
        REF_POS_MATE = int(fields[7])  # Where mate alignment starts.
        # Insert len = aligned pair's end-to-end span along the transcript.
        # Bowtie makes it negative if mate is upstream of this read.
        SPAN = abs(int(fields[8]))

        # samtools view tab-delimited optional fields appear in any order
        ALIGN_SCORE = -1000 # In case none provided, assume a very negative score.
        MISMATCHES = 0
        GAP_OPENS = 0
        GAP_EXTENDS = 0
        EDIT_DIST = 0
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
        # accumulate
        alignment = read_alignment(IS_PRIMARY,SPAN,ALIGN_SCORE,
            EDIT_DIST,MISMATCHES,GAP_OPENS,GAP_EXTENDS,TID)
        if read_group is None:
            # start the first read group
            read_group = pair_alignments(RID)
        elif RID != read_group.get_rid():
            # finish the previous read group
            read_group.show_if_good()
            # and start a new read group
            read_group = pair_alignments(RID)
        # accumulate four alignments per read pair
        if is_parent1:
            if is_read1:
                read_group.add_P1_R1(alignment,TID,ALLELE)
            else:
                read_group.add_P1_R2(alignment,TID,ALLELE)
        else:
            if is_read1:
                read_group.add_P2_R1(alignment,TID,ALLELE)
            else:
                read_group.add_P2_R2(alignment,TID,ALLELE)
    # finish the last read group
    if read_group is not None:
        read_group.show_if_good()
        read_group = None

if __name__ == '__main__':
    process_stdin('MxM','SxS')
