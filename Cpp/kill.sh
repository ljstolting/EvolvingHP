#!/bin/bash

for ((i = 2709597; i <= 2709646; i += 1));
do
  scancel $i
  sleep 3 
done
