import sys
'''
$ samtools view Primary.bam | python bam_to_ml_csv.py > Primary.csv

See https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

MxM is parent 1. SxS is parent 2
'''
class read_alignment():
    def __init__(self,ip,il,score,ed,mm,go,ge):
        self.is_primary = ip
        self.insert_length = il
        self.align_score = score
        self.edit_distance = ed
        self.mismatches = mm
        self.gap_opens = go
        self.gap_extends = ge
    def __str__(self):
        # TO DO: add feature for deviation from mean insert length
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
class pair_alignments():
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
    def show(self):
        good = self.has_four() and \
            len(self.targets)==1 and \
            len(self.alleles)==2
        if not good:
            msg = 'WARN: '+self.rid+\
                ' hit '+str(self.targets)+\
                ' in '+str(self.alleles)
            print(msg, file=sys.stderr)
            return
        msg =  str(self.alignments[0])+','
        msg += str(self.alignments[1])+','
        msg += str(self.alignments[2])+','
        msg += str(self.alignments[3])
        print(msg)

def process_stdin(parent1,parent2):
    read_group = None
    for line in sys.stdin:
        # samtools view tab-delimited required fields appear in standard order
        fields=line.strip().split('\t')
        RID = fields[0]
        FLAGS = int(fields[1])
        IS_PRIMARY = not (FLAGS & 0x100)
        is_read1 = FLAGS & 0x40
        # assume transcript ID like jg27855.t1_SxS
        TID,SUFFIX = fields[2].split('.')
        ALLELE = SUFFIX.split('_')[1]
        if ALLELE==parent1:
            is_parent1 = True
        else:
            is_parent1 = False
            if ALLELE != parent2:
                raise Exception('Unrecognized allele: '+ALLELE)
        TPOS = int(fields[3]) # transcript position
        CIGAR = fields[5]
        MATE_POS = int(fields[7])
        INSERT_LEN = abs(int(fields[8]))
        # samtools view tab-delimited optional fields appear in any order
        ALIGN_SCORE = -1000 # if non provided, assume a very negative score
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
        alignment = read_alignment(IS_PRIMARY,INSERT_LEN,ALIGN_SCORE,
            EDIT_DIST,MISMATCHES,GAP_OPENS,GAP_EXTENDS)
        if read_group is None:
            read_group = pair_alignments(RID)
        elif RID != read_group.get_rid():
            read_group.show()
            read_group = pair_alignments(RID) # move on
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
        if read_group.has_four():
            read_group.show()
            read_group = None
    # handle any stragglers
    if read_group is not None:
        read_group.show()
        read_group = None

if __name__ == '__main__':
    process_stdin('MxM','SxS')
