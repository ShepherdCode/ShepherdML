#!/bin/sh

echo "extract lncRNA"
echo "transcript,gene,name,length" > gencode.v34.lncRNA_transcripts.csv
grep '^>' gencode.v34.lncRNA_transcripts.fa | cut -c 2- | tr '\|' ',' | cut -d ',' -f 1,2,5,7 >> gencode.v34.lncRNA_transcripts.csv
echo "sortest sequences..."
sort -t ',' -k4n gencode.v34.lncRNA_transcripts.csv | head
echo "longest sequences..."
sort -t ',' -k4nr gencode.v34.lncRNA_transcripts.csv | head

echo "extract pcRNA"
echo "transcript,gene,name,length" > gencode.v34.pc_transcripts.csv
grep '^>' gencode.v34.pc_transcripts.fa | cut -c 2- | tr '\|' ',' | cut -d ',' -f 1,2,5,7 >> gencode.v34.pc_transcripts.csv
echo "sortest sequences..."
sort -t ',' -k4n gencode.v34.pc_transcripts.csv | head
echo "longest sequences..."
sort -t ',' -k4nr gencode.v34.pc_transcripts.csv | head
