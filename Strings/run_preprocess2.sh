#!/bin/sh
date
echo "oneline(non-coding)"
python3 fasta_to_oneline.py --delete_N --minlen 200 \
  gencode.v34.lncRNA_transcripts.fa ncRNA.gc34.oneline.fasta
echo "oneline(protein-coding)"
python3 fasta_to_oneline.py --delete_N --minlen 200 \
  gencode.v34.pc_transcripts.fa pcRNA.gc34.oneline.fasta

echo "test_set_aside(non-coding)"
python3 setaside.py ncRNA.gc34.oneline.fasta \
  test_genes.ncRNA.gc34.fasta ncRNA.gc34.unprocessed.fasta
echo "test_set_aside(protein-coding)"
python3 setaside.py pcRNA.gc34.oneline.fasta \
  test_genes.pcRNA.gc34.fasta pcRNA.gc34.unprocessed.fasta

echo "process(non-coding)"
python3 preprocess2.py ncRNA.gc34.unprocessed.fasta ncRNA.gc34.processed.fasta
echo "process(protein-coding)"
python3 preprocess2.py pcRNA.gc34.unprocessed.fasta pcRNA.gc34.processed.fasta
