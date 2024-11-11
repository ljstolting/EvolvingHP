#!/bin/bash

for ((i = 194; i < 195; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
