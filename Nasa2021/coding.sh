#!/bin/bash
python3 GenFASTA.py 1000 True 1 coding_1000_1.fasta 
python3 GenFASTA.py 1000 True 2 coding_1000_2.fasta 
python3 GenFASTA.py 1000 True 3 coding_1000_3.fasta 
cat coding_1000_?.fasta > coding_3000.fasta 
rm coding_1000_1.fasta 
rm coding_1000_2.fasta 
rm coding_1000_3.fasta 
python3 GenFASTA.py 2000 True 1 coding_2000_1.fasta 
python3 GenFASTA.py 2000 True 2 coding_2000_2.fasta 
python3 GenFASTA.py 2000 True 3 coding_2000_3.fasta 
cat coding_2000_?.fasta > coding_6000.fasta 
rm coding_2000_1.fasta 
rm coding_2000_2.fasta 
rm coding_2000_3.fasta 
python3 GenFASTA.py 4000 True 1 coding_4000_1.fasta 
python3 GenFASTA.py 4000 True 2 coding_4000_2.fasta 
python3 GenFASTA.py 4000 True 3 coding_4000_3.fasta 
cat coding_4000_?.fasta > coding_12000.fasta 
rm coding_4000_1.fasta 
rm coding_4000_2.fasta 
rm coding_4000_3.fasta 
python3 GenFASTA.py 8000 True 1 coding_8000_1.fasta 
python3 GenFASTA.py 8000 True 2 coding_8000_2.fasta 
python3 GenFASTA.py 8000 True 3 coding_8000_3.fasta 
cat coding_8000_?.fasta > coding_24000.fasta 
rm coding_8000_1.fasta 
rm coding_8000_2.fasta 
rm coding_8000_3.fasta 
python3 GenFASTA.py 16000 True 1 coding_16000_1.fasta 
python3 GenFASTA.py 16000 True 2 coding_16000_2.fasta 
python3 GenFASTA.py 16000 True 3 coding_16000_3.fasta 
cat coding_16000_?.fasta > coding_48000.fasta 
rm coding_16000_1.fasta 
rm coding_16000_2.fasta 
rm coding_16000_3.fasta 
python3 GenFASTA.py 32000 True 1 coding_32000_1.fasta 
python3 GenFASTA.py 32000 True 2 coding_32000_2.fasta 
python3 GenFASTA.py 32000 True 3 coding_32000_3.fasta 
cat coding_32000_?.fasta > coding_96000.fasta 
rm coding_32000_1.fasta 
rm coding_32000_2.fasta 
rm coding_32000_3.fasta 
python3 GenFASTA.py 64000 True 1 coding_64000_1.fasta 
python3 GenFASTA.py 64000 True 2 coding_64000_2.fasta 
python3 GenFASTA.py 64000 True 3 coding_64000_3.fasta 
cat coding_64000_?.fasta > coding_192000.fasta 
rm coding_64000_1.fasta 
rm coding_64000_2.fasta 
rm coding_64000_3.fasta 
