#!/bin/sh
#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=bowtie
#SBATCH --time=04:00:00   # Bowtie takes about 35 minutes
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads
## source /cluster/bin/jobsetup   ## abel only
set -o errexit # exit on errors
# Our python will generate smaller Aligned.bam which we keep.
#savefile *.bam
#savefile *.db
#savefile *.log
#savefile *.SN

echo MODULES
module --force purge
module load StdEnv 

module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
# Would conflict with Bowtie: Python/3.10.8-GCCcore-12.2.0 SAMtools/1.17-GCC-12.2.0
module list

echo
echo LD_LIBRARY_PATH $LD_LIBRARY_PATH
echo

which python3
python3 --version
#expect MOLBAR_HOME=/cluster/MOLBAR
echo MOLBAR_HOME ${MOLBAR_HOME}
ls ${MOLBAR_HOME}
echo

date
#echo COPY INPUTS TO GRID
#echo RUN THIS FROM THE DIRECTORY CONTAINING R1 AND R2.fastq.gz
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
#echo SCRATCH ${SCRATCH}
#cd ${SCRATCH}
#cp -vHR ${INITIALDIR}/*.fq.gz .
#cp -vHR ${INITIALDIR}/*.bt2 .
#cp -vHR ${INITIALDIR}/*.fasta .
#date
echo

THREADS=4
echo THREADS $THREADS

echo "Expect one fasta file"
FILENAME=*.fasta
BASENAME=` echo ${FILENAME} | sed 's/.fasta//' `
TARGET=${BASENAME}
echo "Found FASTA ${FILENAME}"
echo "Use TARGET ${TARGET}"

R1=*_R1_*.fq.gz
R2=*_R2_*.fq.gz
echo R1 ${R1}
echo R2 ${R2}
echo

ls -l
echo

date
echo run Bowtie
# assumes *.fasta and *.bt2 have same first name
# for homozygous, remove option: --heterozygous
python3 ${MOLBAR_HOME}/src/mapping.py ${TARGET} ${R1} ${R2} Aligned.out.sam --debug
echo -n $?
echo " exit status"
date

echo run python filter
echo "Filter for primary alignments only"
echo "Filter for unspliced alignments only"
python3 ${MOLBAR_HOME}/src/samfilter.py Aligned.out.sam Primary.bam
echo -n $?
echo " exit status"
date

echo sort bam file
# the -T option is critical: uses local directory rather than /tmp
samtools sort --threads $THREADS -T tmp --output-fmt BAM -o Sorted.bam Primary.bam

echo stats
samtools flagstat Sorted.bam > samtools.flagstat
samtools stats Sorted.bam | grep '^SN' | cut -f 2- > samtools.stats.SN

date
ls -l
echo DONE
