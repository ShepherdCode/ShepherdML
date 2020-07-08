#!/bin/sh

echo "Duplicates in ncRNA?"
grep -v '^>' ncRNA.fasta |\
    sort | uniq -c |\
    grep -v '^   1 ' | awk '{print substr($0,1,60);}'

echo "Duplicates in pcRNA?"
grep -v '^>' pcRNA.fasta |\
    sort | uniq -c |\
    grep -v '^   1 ' | awk '{print substr($0,1,60);}'
