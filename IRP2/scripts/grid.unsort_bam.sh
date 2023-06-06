#!/bin/sh

# In every subdirectory, recreate Primary.bam by unsorting Sorted.bam

echo "Submit MxM"
sbatch --account=${ACCOUNT} unsort_bam.sh MxM_BR1
sbatch --account=${ACCOUNT} unsort_bam.sh MxM_BR2
sbatch --account=${ACCOUNT} unsort_bam.sh MxM_BR3
sbatch --account=${ACCOUNT} unsort_bam.sh MxM_BR4

echo "Submit SxS"
sbatch --account=${ACCOUNT} unsort_bam.sh SxS_BR1
sbatch --account=${ACCOUNT} unsort_bam.sh SxS_BR2
sbatch --account=${ACCOUNT} unsort_bam.sh SxS_BR3
sbatch --account=${ACCOUNT} unsort_bam.sh SxS_BR4

echo "Done"
