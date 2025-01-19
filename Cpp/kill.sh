#!/bin/bash

for ((i = 4268223; i <= 4268232; i += 1));
do
  scancel $i
  sleep 3 
done
