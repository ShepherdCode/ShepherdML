#!/bin/sh

#SBATCH --account=$ACCOUNT
#SBATCH --job-name=pilon
#SBATCH --time=8:00:00 # pilon takes about 5 hours
#SBATCH --mem-per-cpu=18G  # pilon can require 18 GB RAM
#SBATCH --cpus-per-task=1  # pilone multi-threading is inefficient
# Jobs submit with 'arrayrun 1-10' receive a $TASK_ID.
# Jobs submit with 'sbatch --array=1-10' recieve a $SLURM_ARRAY_TASK_ID.

## source /cluster/bin/jobsetup   # abel only
set -o errexit # exit on errors
# Pilon expected outputs
#savefile pilon.changes
#savefile pilon.fasta 

# For the very first pilon, if already have directories
# consensus_0 and map_to_consensus_0, then 
# the number of this round should be 1.

if [ $# -ne 2 ]; then
    echo "Please specify <MxM|SxS> and number of this round."
    exit 1
fi
GENOME=$1
echo GENOME $GENOME
ROUND=$2
echo ROUND $ROUND
#let PREV=${ROUND}-1  # works on some unix flavors but is not standard
PREV=$((${ROUND}-1))  # standard unix arithmetic
echo PREV $PREV
PREV_REF=consensus_${PREV}
PREV_MAP=map_to_consensus_${PREV}
echo PREV_REF $PREV_REF
echo PREV_MAP $PREV_MAP
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
echo

echo MODULES
module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load Pilon/1.23-Java-11
module list
# To execute Pilon run: java -Xmx8G -jar $EBROOTPILON/pilon.jar
# Pilon Requirement: fasta file, bam file, bai index
# Pilon Requirement: about 18 GB RAM
JARFILE=$EBROOTPILON/pilon.jar
echo PILON JAR; ls -l $JARFILE
echo HEADER; ls -l $HEADERFILE
echo "Java version is:"
JAVA=`which java`
${JAVA} -showversion |& head -n 4
PILON="${JAVA} -Xmx16G -jar ${JARFILE}"
echo "Pilon command is: ${PILON}"
echo

# Here are some interesting pilon options
OPTIONS="--diploid "  # assume diploid
OPTIONS="--vcf"       # output a VCF
OPTIONS="--variant"   # heuristic for variants not assembly
OPTIONS="--fix all"   # correct bases, indels, and local misassembly and write a FASTA file
OPTIONS="--changes"   # show changes to the FASTA
# Use this set of options to generate the homozygous genome informed by reads.
OPTIONS="--fix all --changes --output ${GENOME}"
echo "Pilon options set to: ${OPTIONS}"
echo

echo "SAMTOOLS MERGE"
echo ${GENOME}
samtools view -b -H ../${PREV_MAP}/${GENOME}_BR1/Sorted.bam > ${GENOME}.header.bam
samtools merge -h ${GENOME}.header.bam -@ 4 -O BAM ${GENOME}.bam \
	../${PREV_MAP}/${GENOME}_BR1/Sorted.bam \
	../${PREV_MAP}/${GENOME}_BR2/Sorted.bam \
	../${PREV_MAP}/${GENOME}_BR3/Sorted.bam \
	../${PREV_MAP}/${GENOME}_BR4/Sorted.bam 
echo

date
echo SAMTOOLS INDEX
samtools index -@ 4 ${GENOME}.bam
echo

echo "Link to old fasta"
ln -s ../${PREV_REF}/${GENOME}.fasta prev.${GENOME}.fasta
ls -l *.fasta
echo

date
echo START PILON
${PILON} --genome prev.${GENOME}.fasta --frags ${GENOME}.bam ${OPTIONS}
echo -n $?; echo " exit status"
echo DONE PILON
date
echo

echo START BOWTIE BUILD
bowtie2-build ${GENOME}.fasta ${GENOME}
date
echo DONE
