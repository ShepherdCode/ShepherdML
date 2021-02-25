python GenFASTA.py 3000 True 1 coding_3000_1.fasta 
python GenFASTA.py 3000 True 2 coding_3000_2.fasta 
python GenFASTA.py 3000 True 3 coding_3000_3.fasta 
cat coding_3000_?.fasta > fake_noncod_9000.fasta 
rm coding_3000_1.fasta 
rm coding_3000_2.fasta 
rm coding_3000_3.fasta 
python GenFASTA.py 6000 True 1 coding_6000_1.fasta 
python GenFASTA.py 6000 True 2 coding_6000_2.fasta 
python GenFASTA.py 6000 True 3 coding_6000_3.fasta 
cat coding_6000_?.fasta > fake_noncod_18000.fasta 
rm coding_6000_1.fasta 
rm coding_6000_2.fasta 
rm coding_6000_3.fasta 
python GenFASTA.py 12000 True 1 coding_12000_1.fasta 
python GenFASTA.py 12000 True 2 coding_12000_2.fasta 
python GenFASTA.py 12000 True 3 coding_12000_3.fasta 
cat coding_12000_?.fasta > fake_noncod_36000.fasta 
rm coding_12000_1.fasta 
rm coding_12000_2.fasta 
rm coding_12000_3.fasta 
python GenFASTA.py 24000 True 1 coding_24000_1.fasta 
python GenFASTA.py 24000 True 2 coding_24000_2.fasta 
python GenFASTA.py 24000 True 3 coding_24000_3.fasta 
cat coding_24000_?.fasta > fake_noncod_72000.fasta 
rm coding_24000_1.fasta 
rm coding_24000_2.fasta 
rm coding_24000_3.fasta 
python GenFASTA.py 48000 True 1 coding_48000_1.fasta 
python GenFASTA.py 48000 True 2 coding_48000_2.fasta 
python GenFASTA.py 48000 True 3 coding_48000_3.fasta 
cat coding_48000_?.fasta > fake_noncod_144000.fasta 
rm coding_48000_1.fasta 
rm coding_48000_2.fasta 
rm coding_48000_3.fasta 
python GenFASTA.py 96000 True 1 coding_96000_1.fasta 
python GenFASTA.py 96000 True 2 coding_96000_2.fasta 
python GenFASTA.py 96000 True 3 coding_96000_3.fasta 
cat coding_96000_?.fasta > fake_noncod_288000.fasta 
rm coding_96000_1.fasta 
rm coding_96000_2.fasta 
rm coding_96000_3.fasta 
python GenFASTA.py 192000 True 1 coding_192000_1.fasta 
python GenFASTA.py 192000 True 2 coding_192000_2.fasta 
python GenFASTA.py 192000 True 3 coding_192000_3.fasta 
cat coding_192000_?.fasta > fake_noncod_576000.fasta 
rm coding_192000_1.fasta 
rm coding_192000_2.fasta 
rm coding_192000_3.fasta 
