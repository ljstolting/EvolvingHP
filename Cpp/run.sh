#!/bin/bash

for ((i = 185; i < 187; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
