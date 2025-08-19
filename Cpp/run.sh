#!/bin/bash
#not run yet
for ((i = 20; i < 128; i += 1));
do
  sbatch --export=JB=$i lindsaygeneral.sh
  sleep 3
done
