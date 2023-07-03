'''
Read one FASTA file and write another.
Remove all but the ID from the defline.
Remove internal newlines from sequences.
'''

import os
import argparse
import gzip

def process_fasta(infilename,outfilename,debug=False):
    cnt = 0
    if infilename.endswith('.gz'):
        with gzip.open (infilename,'r') as fin:
            cnt = write_fasta(fin,outfilename,debug)
    else:
        with open (infilename,'r') as fin:
            cnt = write_fasta(fin,outfilename,debug)
    return cnt

def write_fasta(fin,outfilename,debug=False):
    defline_count = 0
    with open (outfilename,'w') as fout:
        for line in fin:
            # When input is gzipped, lines are class byte.
            # Class byte has a decode() method but str does not.
            try:
                line = line.decode()
            except:
                pass
            line=line.strip() # strip the newline
            if line.startswith('>'):
                if defline_count>0:
                    fout.write('\n') # end prev sequence
                # strip away all but the ID
                defline = line.split(' ')[0]
                fout.write(defline+'\n') # start next sequence
                defline_count += 1
            else:
                fout.write(line)  # sequence continuation without newline
        # end last sequence
        if defline_count>0:
            fout.write('\n') # end the last sequence
    return defline_count

def args_parse():
    global args
    parser = argparse.ArgumentParser(description='prety fasta file')
    parser.add_argument('input_fasta', help = 'filename, read-only', type = str)
    parser.add_argument('output_fasta', help = 'filename, overwrite', type = str)
    parser.add_argument('--debug', action = 'store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    try:
        args_parse()
        infile = args.input_fasta
        outfile = args.output_fasta
        seq_cnt = process_fasta(infile,outfile,args.debug)
        print('Wrote %d sequences.' % seq_cnt)
    except Exception as e:
        print('ERROR! {}'.format(e))
        if args.debug:
            print(traceback.format_exc())
        else:
            print('Run with --debug for traceback.')
