#!/bin/bash

for ((i = 25; i < 50; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
