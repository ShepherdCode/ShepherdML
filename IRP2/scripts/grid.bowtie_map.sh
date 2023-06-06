#!/bin/sh

# Command line parameter gives the number for this round.
# Give a zero for the initial round where we map to CDS
# Give a 1 to map to consensus_1, etc.

if [ $# -eq 0 ]; then
	echo "Please provide the number of this round."
	exit 1
fi
ROUND=$1
echo ROUND $ROUND

#Add path to directory with trimmed reads
TRIM_BASE='/cluster/projects/nn9525k/hybrids/jasonrm/Arenosa/Trimming/reads_combined'
#Add path where you want to run pilon
CONS_BASE='/cluster/work/users/jasonrm/Arenosa/Homozygous'
CONS_DIR="consensus_${ROUND}"
MAP_DIR="map_to_consensus_${ROUND}"
mkdir ${MAP_DIR}
cd ${MAP_DIR}
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
    if [ $DD -lt 4 ]; then
         ln -s ../../${CONS_DIR}/MxM.* .
    else
         ln -s ../../${CONS_DIR}/SxS.* .
    fi
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
    sbatch --account=${ACCOUNT} ../../bowtie_map.sh
    cd ..
done
echo DONE
cd ..
