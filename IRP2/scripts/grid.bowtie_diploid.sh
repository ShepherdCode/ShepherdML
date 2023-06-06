#!/bin/sh

#Add path to directory with trimmed reads
TRIM_BASE='/cluster/projects/nn9525k/hybrids/jasonrm/Arenosa/Trimming/reads_combined'
#Add path where you want to run pilon
pwd

#Add lines for each subfolder based on names of reads
DATA1[0]="MxM_BR1"
DATA1[1]="MxM_BR2"
DATA1[2]="MxM_BR3"
DATA1[3]="MxM_BR4"
DATA1[4]="SxS_BR1"
DATA1[5]="SxS_BR2"
DATA1[6]="SxS_BR3"
DATA1[7]="SxS_BR4"

#Change 'seq x y' so x and y are first and last DATA indices.
for DD in `seq 0 7` ; do
    SAMPLE=${DATA1[${DD}]}
    DIR1="${SAMPLE}"
    rm -rf $DIR1
    mkdir $DIR1
    cd $DIR1
    # link to fasta, bowtie index bt2, traimmed reads fq.gz
    ln -s ../diploid.* .   
    ln -s ${TRIM_BASE}/${SAMPLE}*.fq.gz .
    cd ..
    #
done

echo
echo "GRID SUBMIT"

export ACCOUNT=nn9525k
export MOLBAR_HOME=/cluster/home/jasonrm/Source/MOLBAR

# Visit every subdirectory
for D in */;
do
    cd $D
    pwd
    sbatch --account=${ACCOUNT} ../bowtie_diploid.sh
    cd ..
done
echo DONE
cd ..
