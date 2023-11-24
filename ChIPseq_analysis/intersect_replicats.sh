#!/bin/bash
# $1 : directory name with replicat
DIR=$1
# intersect replicat
bedtools intersect -a "$DIR/$DIR_1_H3K79me2_broad_peak.bed" \
  -b "$DIR/$DIR_2_H3K79me2_broad_peak.bed" \
  -wa > "$DIR/$DIR_H3K79me2_broad_peak_common.bed"
