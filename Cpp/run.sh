#!/bin/bash

for ((i = 5; i < 50; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done