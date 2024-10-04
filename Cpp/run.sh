#!/bin/bash

for ((i = 100; i < 101; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 3 
done
