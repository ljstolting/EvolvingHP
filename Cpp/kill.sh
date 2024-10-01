#!/bin/bash

for ((i = 3759716; i <= 3759721; i += 1));
do
  scancel $i
  sleep 3 
done
