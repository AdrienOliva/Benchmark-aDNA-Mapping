#!/bin/bash
for i in *.sam #change if you work with bamfiles
do
samtools view -F 4 $i > $i.map
sort -u -t$'\t' -k1,1 $i.map > $i.map.unique
echo $i
done
