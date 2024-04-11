#!/bin/bash

for ((i = 2702982; i <= 2702986; i += 1));
do
  scancel $i
  sleep 3 
done
