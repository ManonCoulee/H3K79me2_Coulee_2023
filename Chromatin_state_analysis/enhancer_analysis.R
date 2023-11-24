####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer','ChromENVEE', 'GenomicRanges',"stringr")
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                          INITIALIZATION
####################################################################################################

chromHMM = "~/Documents/Mouse/ChIPseq/ChromHMM"
path_RNAseq = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_SC_SCII_RS_spike/"

tt = theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=30, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 30),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      legend.position="none")

list_dirs = c("Kit-","SC","RS")

genomeFile = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed", sep = "\t",
  h = F)
colnames(genomeFile) = c("chr","start","end","strand","score","geneENS")
# Rada method
# genomeFile = read.table("~/Documents/Mouse/Comparison_DEG_mark/Kit_m/Kit/Kit-_DEG_H3K79me2_GB.tsv",
#   sep = "\t", h = T)[,c(2:4,1)]
# colnames(genomeFile) = c("chr","start","end","geneENS")
genomeFile = makeGRangesFromDataFrame(genomeFile, keep.extra.columns = T)

## Create color table
colorTable_H3 = colorTable
colorTable_H3$stateName = paste0(colorTable_H3$stateName,"_H3K79me2")
colorTable_nH3 = colorTable
colorTable_nH3$stateName = paste0(colorTable_nH3$stateName,"_not_H3K79me2")
color = rbind(colorTable_H3,colorTable_nH3)
col = getStateColor(colorTable = colorTable)

list_files_enhancer_all = list.files(paste0(chromHMM,"/",list_dirs), pattern = "Enh", full.names = T)
list_files_enhancer = list_files_enhancer_all[grepl("H3K79me2",list_files_enhancer_all)]
names_file_enhancer_all = list.files(paste0(chromHMM,"/",list_dirs), pattern = "Enh", full.names = F)
names_file_enhancer = names_file_enhancer_all[grepl("H3K79me2", names_file_enhancer_all)]
names_file_enhancer = unlist(strsplit(names_file_enhancer,"_segments.bed"))

## Order file in function their size
order_dsc_file_size = order(file.info(list_files_enhancer)$size, decreasing = F)
list_files_enhancer_order = list_files_enhancer[order_dsc_file_size]
names(list_files_enhancer_order) = names_file_enhancer[order_dsc_file_size]
# Rada method
# dir = "Kit-"
# list_files_enhancer_order = list_files_enhancer_order[grepl(dir,names(list_files_enhancer_order))]

####################################################################################################
##                                       Enhancer annotation
####################################################################################################

list_table_enhancer = lapply(list_files_enhancer_order, function(file) {
  table = read.table(file, sep = "\t", h = F)
  colnames(table) = c("chr","start","end","chromatinState")
  table$sample = names(list_files_enhancer_order[list_files_enhancer_order == file])
  table$sampleName = paste0(unlist(strsplit(names(list_files_enhancer_order[list_files_enhancer_order == file]),
    "_"))[-c(1:2)], collapse = "_")
  # Rada method
  # table$sampleName = paste0(unlist(strsplit(names(list_files_enhancer_order[list_files_enhancer_order == file]),
  #   "_"))[-c(1:3)], collapse = "_")
  table = makeGRangesFromDataFrame(table, keep.extra.columns = T)
  return(table)
})

list_table_enhancer_gene = lapply(list_table_enhancer, enhancerAnnotation, genome = genomeFile,
  interval = 500000, nCore = 14)

####################################################################################################
##                                       Enhancer expression
####################################################################################################

dir = "SC"
if(dir == "Kit-") path_RNAseq = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_Kit/"

gene_expression = read.table(paste0(path_RNAseq,dir,"_CTL_gene_expression.tsv"),sep = "\t",
  h = T)
colnames(gene_expression)[6:7] = c("geneExpression","geneENS")
# colnames(gene_expression)[4:5] = c("geneExpression","geneENS") ## GS
list_reduce = list_table_enhancer_gene[grepl(dir,names(list_table_enhancer_gene))]

## Associate gene expression to gene-enhancer association
list_table_enhancer_gene_expression = lapply(list_reduce, enhancerExpression,
  geneExpressionTable = gene_expression)

lapply(names(list_table_enhancer_gene_expression), function(name){
  write.table(list_table_enhancer_gene_expression[[name]],
  paste0(name,"_expression_distribution_100kb.tsv"), sep = "\t", col.names = T,
  row.names = F, quote = F)
})

plot_expression = plotEnhancerExpression(list_table_enhancer_gene_expression, scale = "log10",
  colorTable = color, distance = 100000)
plot_expression = plot_expression + theme(axis.text.x = element_text(size=15, angle = 0, hjust = 1))
ggsave(paste0(dir,"_enhancer_expression_distribution_100kb.png"),plot = plot_expression,
  device = "png",width = 8, height = 5, dpi = 600)

# Rada method
# lapply(names(list_table_enhancer_gene_expression), function(name){
#   write.table(list_table_enhancer_gene_expression[[name]],
#   paste0(name,"_expression_distribution_RADA.tsv"), sep = "\t", col.names = T,
#   row.names = F, quote = F)
# })
# plot_expression = plotEnhancerExpression(list_table_enhancer_gene_expression, scale = "log10",
#   colorTable = color)
# ggsave(paste0(dir,"_enhancer_expression_distribution_RADA.png"),plot = plot_expression,
#   device = "png", width = 8, height = 5, dpi = 600)

####################################################################################################
##                                  Distribution of DEG in enhancer
####################################################################################################

DEG = read.table(paste0(path_RNAseq,"DOT1L_KO_effect_within_",dir,"_cells_volcano_table_FDR.tsv"), h = T,
  sep = "\t")[,c(1,14)]
lim = 500000/2
limit = seq(0,lim,length.out = 6)

limit_label = unlist(lapply(1:length(limit), function(l) {
  lab = paste0(str_replace(as.character(as.integer(limit[l])),"000$","kb"),"-",
    str_replace(as.character(as.integer(limit[l+1])),"000$","kb"))
  if((l+1) > length(limit)) {
    lab = paste0(">",str_replace(as.character(limit[l]),"000$","kb"))
  }
  return(lab)
}))

table_DEG = do.call("rbind",lapply(names(list_table_enhancer_gene_expression),function(name) {
  x = list_table_enhancer_gene_expression[[name]]
  x = getInformation(x)
  x_DEG = merge(x,DEG, by.x = "gene_name", by.y = "gene_ENS")

  x_DEG$distance = as.numeric(x_DEG$distance)
  for(l in limit) {
    pos = which(x_DEG$distance > l)
    x_DEG[pos,"distance_red"] = limit_label[limit == l]
  }
  x_DEG$state = factor(x_DEG$chromatin_state, levels = colorTable$stateNumber,
    labels = colorTable$stateName)
  x_DEG$H3K79me2 = unlist(strsplit(unique(name),paste0(unique(x_DEG$state),"_")))[2]
  x_DEG$sample_name = name
  x_DEG = x_DEG[x_DEG$distance <= 100000,]
  return(x_DEG)
}))

write.table(table_DEG,paste0(dir,"_DEG_H3K79me2_ENH_100kb_FDR.tsv"), sep = "\t", col.names = T,
  row.names = F, quote = F)
