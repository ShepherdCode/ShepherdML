#!/bin/sh

SRC="/Users/jasonmiller/Source/Python/ShepherdML/Strings"

date
echo "Get 6-mer counts nc"
time python3 ${SRC}/fasta_to_feature.py --kmax 6 --kmin 2 ncRNA.fasta ncRNA
date

echo "Get 6-mer counts pc"
time python3 ${SRC}/fasta_to_feature.py --kmax 6 --kmin 2 --base 50000 pcRNA.fasta pcRNA
date

