#!/bin/bash

## make a dummy design file

READSET=ReadSet.txt

echo -e "sample\tContrast" > design.txt

tail -n +2 $READSET | cut -f1 |
while read line
do

	echo -e "$line\t1"

done>>design.txt

