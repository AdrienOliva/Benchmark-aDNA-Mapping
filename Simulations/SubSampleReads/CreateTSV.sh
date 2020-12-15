#!/bin/bash
for i in *.sam.map.unique
do
cut -f1,10 $i > $i.tsv
echo $i
done
