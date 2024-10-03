#!/bin/bash

for ((i = 6; i < 7; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
