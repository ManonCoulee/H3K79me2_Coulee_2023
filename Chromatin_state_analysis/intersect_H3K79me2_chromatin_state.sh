#! /bin/bash

#conda activate tools

dir=$1
H3K79me2=$2

files=$(ls $dir)

for state in $files
do
	echo $state
	IFS='_' read -r -a array <<< $state
	res=${array[0]}_${array[1]}_${array[2]}_H3K79me2_${array[3]}
	res_not=${array[0]}_${array[1]}_${array[2]}_not_H3K79me2_${array[3]}

## Bedtools intersect between chromatin state and H3K79me2
# Keep common region
	bedtools intersect -a "${dir}/${state}" -b ${H3K79me2} -wa | uniq > "${dir}/${res}"
# Keep region whitout H3K79me2 peaks
	bedtools intersect -a "${dir}/${state}" -b ${H3K79me2} -wa -v | uniq > "${dir}/${res_not}"
done
