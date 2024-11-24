#!/bin/bash

#SBATCH -J lindsay
#SBATCH -p general
#SBATCH -o lindsay_%j.txt
#SBATCH -e lindsay_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH -A r00213

# mkdir ./Generalist\ HP\ Mechanisms/$JB;
cd ./Specifically\ Evolved\ HP\ mechanisms/Every\ Circuit/$JB;
# echo "ND" $JB
time ../../../main $JB
cd ../../;

