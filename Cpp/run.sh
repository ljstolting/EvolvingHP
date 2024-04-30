#!/bin/bash

for ((i = 150; i < 182; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
