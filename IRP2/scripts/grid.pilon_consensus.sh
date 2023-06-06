#!/bin/sh

# Launch pilon on the grid.

# Command line parameter 1 gives the number for this round.
# For example, on the very first pilon run,
# you should have consensus_0 and map_to_consensus_0 directories,
# and you should supply the parameter 1.
# Then, script would make directory consensus_1.

if [ $# -eq 0 ]; then
    echo "Please provide the number of this round."
    exit 1
fi
ROUND=$1
WORKDIR=consensus_${ROUND}
echo WORKDIR $WORKDIR

mkdir "$WORKDIR"
# As a precaution, refuse to overwrite an existing directory.
# The usual culprit is user provided the wrong round.
if [ $? -ne 0 ] ; then
    echo "FAIL: Directory already exists."
    exit 2
fi

export ACCOUNT="nn9525k"
# It is critical to launch the script in the subdirectory.
# The launch determines where the pilon.fasta result gets written.
# Launching two jobs in one directory leads to overrittedn files.

cd $WORKDIR
sbatch --account=${ACCOUNT} ../pilon_consensus.sh MxM $ROUND
sbatch --account=${ACCOUNT} ../pilon_consensus.sh SxS $ROUND
cd ..

