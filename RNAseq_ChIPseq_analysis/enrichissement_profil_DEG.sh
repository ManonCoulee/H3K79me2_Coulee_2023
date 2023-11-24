#/bin/bash

dir=$1

##compute matrix
computeMatrix scale-regions \
	-S ~/Documents/Mouse/ChIPseq/H3K79me2/${dir}/${dir}_1_H3K79me2.bw \
	~/Documents/Mouse/ChIPseq/H3K79me2/${dir}/${dir}_2_H3K79me2.bw \
	-R ${dir}_KO_DR.bed \
	${dir}_KO_UR.bed \
	-a 1500 -b 1500 -m 3000 \
	-p 14 \
	-out ${dir}_DEG_expression_matrix_H3K79me2_NR.gz

plotProfile -m ${dir}_DEG_expression_matrix_H3K79me2_NR.gz \
	-out ${dir}_DEG_expression_profil_H3K79me2_NR.png \
	-z DR UR \
	--averageType mean \
	--colors blue red grey

plotHeatmap -m ${dir}_DEG_expression_matrix_H3K79me2_NR.gz \
	-out ${dir}_DEG_expression_heatmap_H3K79me2_NR.png \
	--whatToShow 'heatmap only' \
	-z DR UR

