#!/bin/sh

# Process every bam file of non-hybrid reads mapped to simulated diploid transcriptome.
# Create text files of maps per read.
# The --debug flag generates stats on stderr.

echo "Submit MxM"
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/MxM_BR1
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/MxM_BR2
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/MxM_BR3
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/MxM_BR4

echo "Submit SxS"
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/SxS_BR1
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/SxS_BR2
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/SxS_BR3
sbatch --account=${ACCOUNT} run_diploid_diff.sh ../Diploid/SxS_BR4

echo "Done"
