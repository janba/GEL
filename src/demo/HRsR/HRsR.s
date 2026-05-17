#!/bin/bash

#SBATCH --job-name=HRsR
#SBATCH --output=HRsR-%J.out
#SBATCH --cpus-per-task=20
#SBATCH --time=36:00:00
#SBATCH --mem=256gb
#SBATCH --mail-user=ruicu@dtu.dk
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=titans
#SBATCH --export=ALL

## INFO
echo "Node: $(hostname)"
echo "Start: $(date +%F-%R:%S)"
echo -e "Working dir: $(pwd)\n"

SCRATCH=/scratch/$USER
if [[ ! -d $SCRATCH ]]; then
  mkdir $SCRATCH
fi

source ~/.bashrc
./build/HRsR

echo "Done: $(date +%F-%R:%S)"