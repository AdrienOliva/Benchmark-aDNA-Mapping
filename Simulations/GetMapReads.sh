#!/bin/bash
for i in *.sam #was working with sam files at the time
do
samtools view -F 4 $i > $i.map
sort -u -t$'\t' -k1,1 $i.map > $i.map.unique
echo $i
done
