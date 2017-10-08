#!/bin/sh

for i in 1600 16000 160000 1600000 16000000;
do
  #for a in "NPO" "PRO" "PMJ_7";
  for a in "PMJ_6" "PRO";
  do
    echo "\n----------------- $i $a\n"
    numactl --localalloc ./src/mchashjoins -r $i -s 128000000 --basic-numa -n 1 -t 10 -a $a
    echo "\n----------------- $i $a\n"
  done;
done
