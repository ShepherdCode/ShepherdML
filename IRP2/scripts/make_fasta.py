'''
Write two fasta files.
One for lyrata, one for halleri.
Both have the same number of sequences.
Each sequence in one has exactly one match in the other.
Both IDs are recorded in each defline.
For ease of processing, remove internal newlines from sequences.
'''

map_L_to_H = dict()
map_H_to_L = dict()
match_file = 'best_matches.lyrata_len_halleri.txt'

def load_matches(filename):
    with open (filename,'r') as fin:
        for line in fin:
            line=line.strip()
            fields=line.split('\t')
            L = fields[0]
            length = fields[1]
            H = fields[2]
            map_L_to_H[L]=H
            map_H_to_L[H]=L

load_matches(match_file)

#print('Just testing:')
#test_L = 'scaffold_81700001.1'
#print(test_L, 'maps to', map_L_to_H[test_L])
#test_H = 'g10924.t1'
#print(test_H, 'maps to', map_H_to_L[test_H])

def write_fasta(infilename,outfilename,allele):
    with open (infilename,'r') as fin:
        with open (outfilename,'w') as fout:
            first_line = True
            defline = None
            for line in fin:
                line=line.strip()
                if line.startswith('>'):
                    if defline is not None:
                        if first_line:
                            first_line=False
                        else:
                            fout.write('\n') # end prev sequence
                    defline = None
                    tid = line.split(' ')[0][1:]
                    if allele=='halleri':
                        if tid in map_H_to_L.keys():
                            # Write tid of best match in lyrata just for tracking
                            hid = tid
                            lid = map_H_to_L[tid]
                            defline = '>'+hid+' halleri matches lyrata:'+lid+'\n'
                    else:
                        if tid in map_L_to_H.keys():
                            # Use tid of best match in halleri to emphasize same gene
                            # Write original lyrata tid just for tacking
                            hid = map_L_to_H[tid]
                            defline = '>'+hid+' lyrata matches halleri:'+tid+'\n'
                    if defline is not None:
                        # We intentionally remove newlines from sequence
                        fout.write(defline) # start next sequence
                else:
                    if defline is not None:
                        fout.write(line)  # more sequence
            # end last sequence
            if defline is not None:
                fout.write('\n')
                defline = None

print('Writing lyrata...')
lyrata_in = 'Arabidopsis_lyrata.v.1.0.cdna.all.fa'
lyrata_out = 'lyrata.fasta'
write_fasta(lyrata_in,lyrata_out,'lyrata')
print('Writing halleri...')
halleri_in = 'Arabidopsis_halleri.Ahal2.2.cdna.all.fa'
halleri_out = 'halleri.fasta'
write_fasta(halleri_in,halleri_out,'halleri')
print('Wrote two fasta files.')
