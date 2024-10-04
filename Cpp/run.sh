#!/bin/bash

for ((i = 108; i < 109; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
