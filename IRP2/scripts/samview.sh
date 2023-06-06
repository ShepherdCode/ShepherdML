# Reduce filesize by cutting away each read info (ID, bases, base qualities).
# This takes about 15 minutes.
samtools view Sorted.bam | cut -f 2-9,12- > samview.tsv
