python3 GenFASTA.py 16667 True 1 coding_16k_1.fasta
python3 GenFASTA.py 16667 True 2 coding_16k_2.fasta
python3 GenFASTA.py 16666 True 3 coding_16k_3.fasta
python3 GenFASTA.py 50000 False 3 noncod_validation.fasta
cat coding_16k_?.fasta > coding_validation.fasta
rm coding_16k_1.fasta
rm coding_16k_2.fasta
rm coding_16k_3.fasta
