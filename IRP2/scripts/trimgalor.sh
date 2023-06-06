#!/bin/bash
#SBATCH --job-name=trim1
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=6
##SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
##SBATCH --partition=normal

module purge
#module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
module load Trim_Galore/0.6.7-GCCcore-10.3.0
trim_galore $1 $2 --paired -j 6
