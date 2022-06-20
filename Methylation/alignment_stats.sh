#!/bin/bash


files=$(ls alignments/*mapstats)
arr=($files)

join ${arr[@]:0:2} > tmp.tmp
for f in ${arr[@]:2}
do
	join tmp.tmp $f > tmpf
	mv tmpf tmp.tmp

done

echo "Cat" $files > align_tmp
head -n 25 tmp.tmp | grep -E -v "^pairs|^mapped|mate1|mate2" >> align_tmp

rm tmp.tmp
