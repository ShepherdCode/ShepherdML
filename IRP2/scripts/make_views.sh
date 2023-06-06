samtools view Sorted.bam | cut -f 2-9,12- > samview.tsv
python ../samview_to_csv.py < samview.tsv > alignment_summaries.txt
rm samview.tsv
