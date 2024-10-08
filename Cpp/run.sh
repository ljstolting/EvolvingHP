#!/bin/bash

for ((i = 300; i < 310; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
