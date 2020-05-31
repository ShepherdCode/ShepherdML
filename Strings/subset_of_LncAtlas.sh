#/bin/sh

# Extract a subset from LncAtlas.
# + non-coding (nc)
# + one cell line (HeLa.S3)
# + one score type (cyto vs nuclear RCI)

cat LncAtlas.download_all_raw_data.csv | grep "nc,nc" | grep "HeLa.S3,CNRCI" | grep -v "CNRCI,NA," | tr ',' '\t' > non-coding.Hela.CNRCI.tsv
