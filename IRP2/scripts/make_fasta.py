'''
Write two fasta files.
Both have the same number of sequences.
Each sequence in one has exactly one match in the other.
Both IDs are recorded in each defline.
For ease of processing, remove internal newlines from sequences.
'''

def load_best(filename):
    best_forward = dict()
    best_reverse = dict()
    with open (filename,'r') as fin:
        for line in fin:
            line=line.strip()
            fields=line.split('\t')
            acov1=fields[7]
            acov2=fields[8]
            tid1=fields[9]
            tid2=fields[10]
            alignment1=(acov1,tid2)
            alignment2=(acov2,tid1)
            if tid1 not in best_forward.keys():
                best_forward[tid1]=alignment1
            else:
                prev = best_forward[tid1]
                if prev is not None:
                    if prev[0]==alignment1[0]:
                        #print('Exclude tie',prev,alignment1)
                        best_forward[tid1] = None # exclude genes with ties
                    elif prev[0] < alignment1[0]:
                        best_forward[tid1] = alignment1
            if tid2 not in best_reverse.keys():
                best_reverse[tid2]=alignment2
            else:
                prev = best_reverse[tid2]
                if prev is not None:
                    if prev[0] == alignment2[0]:
                        #print('Exclude tie',prev,alignment2)
                        best_reverse[tid2] = None # exclude genes with ties
                    elif prev[0] < alignment2[0]:
                        best_reverse[tid2] = alignment2
    return best_forward,best_reverse

def require_reflection(dic1,dic2):
    kills = list()
    for	orig in dic1.keys():
        if dic1[orig] is None:
            kills.append(orig)
        else:
            fwd_len,fwd_tid = dic1[orig]
            if fwd_tid not in dic2.keys():
                kills.append(orig)
            elif dic2[fwd_tid] is None:
                kills.append(orig)
            else:
                rev_len,rev_tid = dic2[fwd_tid]
                if rev_tid != orig:
                    kills.append(orig)
    for k in kills:
        dic1.pop(k)

best_forward,best_reverse = load_best('B6_D2.coords')
require_reflection(best_forward,best_reverse)
require_reflection(best_reverse,best_forward)

def scrub(dic):
    new_dict = dict()
    for k in dic:
        (alen,tid) = dic[k]
        new_dict[k]=tid
    return new_dict

map_B_to_D = scrub(best_forward)
map_D_to_B = scrub(best_reverse)


def write_fasta(infilename,outfilename,allele):
    with open (infilename,'r') as fin:
        with open (outfilename,'w') as fout:
            defline_count = 0
            is_defline = False
            write_this = False
            for line in fin:
                line=line.strip()
                is_defline = line.startswith('>')
                if is_defline:
                    tid = line.split(' ')[0][1:]
                    # For both strains, write B tid plus D tid for tracking
                    write_this = False
                    if allele=='B6' and tid in map_B_to_D.keys():
                            match = map_B_to_D[tid]
                            defline = '>'+tid+' matches '+match+'\n'
                            write_this = True
                    if allele=='D2' and tid in map_D_to_B.keys():
                            match = map_D_to_B[tid]
                            defline = '>'+tid+' matches '+match+'\n'
                            write_this = True
                    if write_this:
                        if defline_count>0:
                            fout.write('\n') # end prev sequence
                        fout.write(defline) # start next sequence
                        defline_count += 1
                elif write_this:
                    fout.write(line)  # sequence continuation
            # end last sequence
            if write_this:
                fout.write('\n') # end the last sequence
    print('Wrote %d sequences' % defline_count)

if True:
    print('Writing B6...')
    B6_in = 'B6.cdna.all.fa'
    B6_out = 'B6.fasta'
    write_fasta(B6_in,B6_out,'B6')
    print('Writing D2...')
    D2_in = 'D2.cdna.all.fa'
    D2_out = 'D2.fasta'
    write_fasta(D2_in,D2_out,'D2')
    print('Wrote two fasta files.')
