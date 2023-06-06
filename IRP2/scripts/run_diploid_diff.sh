#!/bin/sh

# Process bam files of non-hybrid reads mapped to simulated diploid transcriptome.
# Create text files of maps per read.
# The --debug flag generates stats on stderr.

#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=dip_dif
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
pwd
echo 'PROCESS THIS DIRECTORY:'
echo $1
date
echo differential_mapping.py
python3 ${SRC}/differential_mapping.py  ${1}/Primary.bam --debug 1> ${1}.tmp 2> ${1}.report

date
echo count_diffmap_per_gene.py
python3 ${SRC}/count_diffmap_per_gene.py  ${1}.tmp > ${1}.tsv

echo "You may zip or delete the tmp files."
echo "DONE"
date
