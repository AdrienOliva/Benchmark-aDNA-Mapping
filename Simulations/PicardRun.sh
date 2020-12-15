#!/bin/bash
module load picard

IN=$1
OUT=$2
KEEP=$3

java -jar $EBROOTPICARD/picard.jar FilterSamReads I=$IN O=$OUT READ_LIST_FILE=$KEEP  FILTER=includeReadList
