#!/bin/sh

SRC="/Users/jasonmiller/Source/Python/ShepherdML/Strings"
echo "Input files (must be unzipped):"
ls -l gencode.*.fa

echo "Clean up lncRNA fasta"
python3 ${SRC}/gencode_preprocess.py gencode.v34.lncRNA_transcripts.fa ncRNA.fasta

echo "Clean up prot coding fasta"
python3 ${SRC}/gencode_preprocess.py gencode.v34.pc_transcripts.fa pcRNA.fasta
