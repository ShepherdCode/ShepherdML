#!/bin/sh

# We accidentally deleted Primary.bam
# Here, unsort Sorted.bam to recerate Primary.bam

#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=unsort
#SBATCH --time=04:00:00   
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads
set -o errexit # exit on errors

echo MODULES
module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module list

# export SRC=/cluster/home/jasonrm/Source/MOLBAR/src
export SRC='.'    # use local copy for now
echo 'PROCESS THIS DIRECTORY:'
echo $1
cd $1
pwd
date
echo 'SORTING...'
samtools sort -n -@ 4 -o Primary.bam Sorted.bam
echo "exit status" $?

# Tried collate. It did not work. I generatd empty bam files.
# collate is faster than sort - it groups by read in random order
# collate -f is faster still - cannot use it on a "best 2" bam file
#samtools collate -@ 4 -o Primary.bam Sorted.bam
echo "DONE"
date
