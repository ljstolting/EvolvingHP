#!/bin/bash

#SBATCH -J lindsay
#SBATCH -p general
#SBATCH -o lindsay_%j.txt
#SBATCH -e lindsay_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=20:00:00
#SBATCH --mem=16G
#SBATCH -A r00213

# cd ./Specifically\ Evolved\ HP\ mechanisms/Every\ Circuit/$JB;
# cd ./No\ Timing\ Requirements/$JB;
# ../../mainHPparspace $JB;
# cd ../../;

mkdir ./Test3DHPonPyloricSolutions/predictedHPperformanceslice/$JB;
cd ./Test3DHPonPyloricSolutions/predictedHPperformanceslice/$JB
../../../mainavgpyloriccombo $JB;


