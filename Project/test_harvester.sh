#!/bin/sh

SRC="/Users/jasonmiller/Source/Python/ShepherdML/Strings"

date
echo "Get 6-mer counts nc"
time python3 ${SRC}/fasta_to_feature.py --kmax 6 --kmin 6 ncRNA.fasta test
date

echo "Get 4-mer counts nc"
time python3 ${SRC}/fasta_to_feature.py --kmax 4 --kmin 4 ncRNA.fasta test
date

echo "Get 2-mer counts nc"
time python3 ${SRC}/fasta_to_feature.py --kmax 2 --kmin 2 ncRNA.fasta test
date



