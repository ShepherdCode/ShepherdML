#!/bin/sh

echo Download cDNA transcriptomes

# wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_lyrata/cdna/Arabidopsis_lyrata.v.1.0.cdna.all.fa.gz
# wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_halleri/cdna/Arabidopsis_halleri.Ahal2.2.cdna.all.fa.gz
# gunzip *.fa.gz

echo Align transcripts

# nucmer Arabidopsis_halleri.Ahal2.2.cdna.all.fa Arabidopsis_lyrata.v.1.0.cdna.all.fa
# mv out.delta halleri_lyrata.delta

echo Select allele pairs

# Use show-coords from nucmer/mummer.
# Reduce the alignments to at most one per pair of alleles.
# Insist the alignment covers 75% of both transcripts.
# Sort the halleri length reverse, and keep the longest alignment per halleri id.
# Then do the same for lyrata.
# The result is 15093 pairs of IDs and each ID occurs once at most.
./show-coords -H -T -c halleri_lyrata.delta \
    | grep '\.t1' \
    | awk '{if ($7>=75 && $8>=75) print $0;}' \
    | awk '{printf "%10s\t%6d\t%s\n", $10, $5, $11;}' \
    | sort -r \
    | awk '{if ($1!=p) printf "%30s\t%6d\t%s\n", $3, $2, $1; p=$1;}' \
    | sort -r \
    | awk '{if ($1!=p) print $0; p=$1;}' \
    > best_matches.lyrata_len_halleri.txt

echo Create fasta files formatted for IRP scripts

python make_fasta.py

echo Create the synthetic diploid

cat lyrata.fasta | sed 's/>\(.*\) \(.*\)/>\1_lyrata \2/' > diploid.fasta
cat halleri.fasta | sed 's/>\(.*\) \(.*\)/>\1_halleri \2/' >> diploid.fasta

echo Compress files

gzip -v *.fa *.fasta

echo Done
