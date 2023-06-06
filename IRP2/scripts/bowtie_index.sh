#! /bin/sh
# Use a module command to enable bowtie2.
# Run this once, then map reads from multiple libraries.
# Do not repeat this operation within every mapping job.
date
bowtie2-build consensus.fasta consensus
date
ls -l
echo 'Done'
