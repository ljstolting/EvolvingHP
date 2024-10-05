#!/bin/bash

for ((i = 205; i < 206; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
