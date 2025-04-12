#!/bin/bash

for ((i = 55; i < 100; i += 1));
do
  sbatch --export=JB=$i lindsayrandom.sh
  sleep 3
done
