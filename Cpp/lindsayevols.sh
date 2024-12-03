#!/bin/bash

#SBATCH -J lindsay
#SBATCH -p general
#SBATCH -o lindsay_%j.txt
#SBATCH -e lindsay_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH -A r00213

# mkdir ./Generalist\ HP\ Mechanisms/$JB;
cd  ./Specifically\ Evolved\ HP\ mechanisms/Every\ Circuit/$JB; #./Generalist\ HP\ Mechanisms/$JB;
# echo "ND" $JB

# for ((n = 0; n < 5; n += 1));
# do
#     mkdir ./$n;
#     cd ./$n;
#     echo "ND" $JB;
#     time ../../../../main $JB;
#     cd ../;
# done

time ../../../main $JB

cd ../../../;

