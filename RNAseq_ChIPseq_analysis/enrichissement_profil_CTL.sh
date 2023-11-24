#/bin/bash

dir=$1

##compute matrix
computeMatrix scale-regions \
	-S ~/Documents/Mouse/ChIPseq/H3K79me2/${dir}/${dir}_1_H3K79me2.bw \
	~/Documents/Mouse/ChIPseq/H3K79me2/${dir}/${dir}_2_H3K79me2.bw \
	-R ${dir}_CTL_grp1.bed \
	${dir}_CTL_grp2.bed \
	${dir}_CTL_grp3.bed \
	${dir}_CTL_grp4.bed \
	${dir}_CTL_grp5.bed \
	${dir}_CTL_grp6.bed \
	${dir}_CTL_grp7.bed \
	${dir}_CTL_grp8.bed \
	${dir}_CTL_grp9.bed \
	${dir}_CTL_grp10.bed \
	-a 1500 -b 1500 -m 3000 \
	-p 14 \
	-out ${dir}_CTL_expression_matrix.gz

plotProfile -m ${dir}_CTL_expression_matrix.gz \
	-out ${dir}_CTL_expression_profil.png \
	-z grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8 grp9 grp10
