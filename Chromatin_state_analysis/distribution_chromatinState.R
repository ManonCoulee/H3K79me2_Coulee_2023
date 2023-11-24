####################################################################################################
##                                          LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('devtools', 'ggplot2', 'GenomicRanges', 'parallel', 'devtools', 'stringr','ChromENVEE')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                        INITIALIZATION
####################################################################################################

## Path of data
wd = "~/Documents/Mouse/Comparison_H3K79me2_enhancers_state"
chromHMM = "~/Documents/Mouse/ChIPseq/ChromHMM"
path_H3K79me2 = "~/Documents/Mouse/ChIPseq/H3K79me2/"

data(colorTable)

list_dirs = c("Kit-","Kit_p","SC","RS")

####################################################################################################
##                                 Chromatin state distribution
####################################################################################################

## Load of all segments file
list_file_chromatin_state = paste0(chromHMM,"/model_18_makeSegmentation_AGSC_SC_RS/",list_dirs,
  "_18_segments.bed")
names(list_file_chromatin_state) = c("Kit_m","Kit_p","RS","SC")

## Separate file according to presence of absence of H3K79me2
list_file_H3K79me2 = list.files(paste0(chromHMM,"/",list_dirs),
  pattern = "18_H3K79me2_segments.bed", full.names = T)
list_file_notH3K79me2 = list.files(paste0(chromHMM,"/",list_dirs),
  pattern = "18_not_H3K79me2_segments.bed", full.names = T)
names(list_file_H3K79me2) = c("GS","RS","SC")
names(list_file_notH3K79me2) = c("GS","RS","SC")

## Transform the data frame to GRanges
createGRanges = function(name, files){
  file = files[[name]]
  df = read.table(file, sep = "\t", h = F)
  colnames(df) = c("chr","start","end","state")
  df$name = name
  df = makeGRangesFromDataFrame(df, keep.extra.columns = T)
  return(df)
}
list_data_chromatin_state = lapply(names(list_file_chromatin_state), createGRanges,
  files = list_file_chromatin_state)
list_data_H3K79me2 = lapply(names(list_file_H3K79me2), createGRanges,
  files = list_file_H3K79me2)
list_data_notH3K79me2 = lapply(names(list_file_notH3K79me2), createGRanges,
  files = list_file_notH3K79me2)

## Generate plot of chromatin state distribution
summary_chromatin_state = plotChromatinState(list_data_chromatin_state, merge = T, plot = T,
  colorTable = colorTable, filename = "./Kit_SC_RS_chromatin_state_distribution.png")
with = plotChromatinState(list_data_H3K79me2, merge = T, plot = T, colorTable = colorTable,
  filename = "./Kit_SC_RS_chromatin_state_distribution_with_H3K79me2.png")
write.table(with,"Kit_SC_RS_chromatin_state_distribution_with_H3K79me2.tsv", sep = "\t",
  col.names = T, row.names = F, quote = F)
without = plotChromatinState(list_data_notH3K79me2, merge = T, plot = T, colorTable = colorTable,
  filename = "./Kit_SC_RS_chromatin_state_distribution_without_H3K79me2.png")
write.table(without,"Kit_SC_RS_chromatin_state_distribution_without_H3K79me2.tsv", sep = "\t",
  col.names = T, row.names = F, quote = F)
