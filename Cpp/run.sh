#!/bin/bash

for ((i = 200; i < 210; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
