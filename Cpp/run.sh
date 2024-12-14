#!/bin/bash

for ((i = 10; i < 100; i += 1));
do
  sbatch --export=JB=$i lindsayevols.sh
  sleep 4 
done
