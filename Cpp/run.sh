#!/bin/bash

for ((i = 10; i < 600; i += 1));
do
  sbatch --export=JB=$i lindsayrandom.sh
  sleep 3
done
