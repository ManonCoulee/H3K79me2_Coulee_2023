####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('UpSetR','ggplot2', 'gtable','stringr','RColorBrewer', 'edgeR', 'DESeq2',
  'VennDiagram','gridExtra','pheatmap','factoextra','ggalluvial')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                       INITIALIZATION
####################################################################################################

dir = c("../DE_analysis_3_factors_Kit","../DE_analysis_3_factors_SC_SCII_RS_spike")
rout = file.path(".")
samples = c("Kit-","Kit+","RS","SC","SCII")

####################################################################################################
##                                           DATA
####################################################################################################

list_files = list.files(c(dir),
  pattern = "cells_volcano_table_FDR.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})

chrom = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", sep = ";", h = T)[,1:2]
chrom$Chromosome = paste0("chr",chrom$Chromosome)
centro = read.table("~/Documents/Annotations/mouse_centromere_location.tsv",sep = "\t",
  h = F)[,c(2:4,8)]
colnames(centro) = c("chr","start","end","gl")


####################################################################################################
##                                DEG gene localisation
####################################################################################################

lapply(names(list_data), function(cell){
  data = list_data[[cell]]
  data_anno = data[data$gl != "Not regulated",c(2:4,14)]
  data_anno$gl = factor(data_anno$gl, labels = c("DR","UR"))
  data_anno = rbind(data_anno,centro)
  data_anno = merge(data_anno, chrom, by.x = "chr", by.y = "Chromosome", all.x = T)

  p = ggplot(data_anno) +
    geom_bar(aes(x = Total_length, y = factor(chr, levels = chrom$Chromosome)),
      stat = "identity", position = position_dodge(),
      fill = "white", color = "grey") +
    geom_point(aes(x = end, y = chr, color = gl), alpha = 0.5, size = 3) +
    scale_color_manual(values = c("blue","red","black","green")) +
    xlab("chromosome size (bp)") + ylab("") + labs(color = "") +
    theme_bw() + theme(
      text = element_text(size=35, angle = 0),
      legend.position = "bottom",
      axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank())

  ggsave(paste0(cell,'_DEG_gene_position_chr.png'),  plot = p, width = 10, height = 8,
    device = 'png', dpi = 600)
})

####################################################################################################
##                                      ALLUVIAL
####################################################################################################

## Keep only unique gene betwenn all cell type
data_gene = unique(rbind(list_data$SCI[,1:5],list_data$SCII[,1:5],list_data$RS[,1:5],
  list_data$'Kit-'[,1:5],list_data$'Kit+'[,1:5]))
rownames(data_gene) = data_gene$gene_ENS

## Create a function to process alluvial plot
aluvial_plot = function(matrix, fill_color = "") {
  plot = ggplot(matrix,
    aes(y = freq, axis1 = factor(`Kit-`), axis2 = factor(`Kit+`), axis3 = factor(SCI),
      axis4 = factor(SCII), axis5 = factor(RS), fill = factor(motif), color = "1")) +
    geom_alluvium() +
    geom_stratum(width = 1/12, fill = fill_color, color = fill_color) +
    scale_x_discrete(limits = c("Kit-","Kit+","SCI","SCII","RS"), expand = c(.1, .1)) +
    scale_fill_grey(start = 1, end = 0) +
    scale_color_manual(values = "black") +
    ylab("") +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=25, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none")
    return(plot)
}

## Create a comparison matrix with all conditions
comparison_matrix = expand.grid('Kit-' = c("Up-regulated","Down-regulated","Not regulated"),
  'Kit+' = c("Up-regulated","Down-regulated","Not regulated"),
  SCI = c("Up-regulated","Down-regulated","Not regulated"),
  SCII = c("Up-regulated","Down-regulated","Not regulated"),
  RS = c("Up-regulated","Down-regulated","Not regulated"))
comparison_matrix$motif = paste0(
  comparison_matrix$`Kit-`,
  comparison_matrix$`Kit+`,
  comparison_matrix$SCI,
  comparison_matrix$SCII,
  comparison_matrix$RS)
rownames(comparison_matrix) = comparison_matrix$motif
comparison_matrix$freq = 0

## Return deregulated gene level for each cell type
gene_table = do.call("cbind",lapply(names(list_data), function(x) {
  data = list_data[[x]]
  rownames(data) = data$gene_ENS
  data_gene$gl = "Not regulated"
  data_gene[data$gene_ENS,"gl"] = data$gl
  data_gene$name = x
  return(data_gene)
}))
gene_table = gene_table[,colnames(gene_table) == "gl"]
colnames(gene_table) = names(list_data)
gene_table = cbind(data_gene,gene_table)

## Generate motif corresponding to dynamic
gene_table$motif = paste0(
  gene_table$`Kit-`,
  gene_table$`Kit+`,
  gene_table$SCI,
  gene_table$SCII,
  gene_table$RS)

write.table(gene_table, "Kit_SC_SCII_RS_DEG_dynamic.tsv", col.names = T, row.names = F,
  sep  ="\t", quote = F)

## Count the frequence of each dynamic motif
frequence_motif = table(gene_table$motif)
comparison_matrix[names(frequence_motif),"freq"] = frequence_motif

## Remove empty motif
comparison_matrix = comparison_matrix[comparison_matrix$freq > 10 ,]
## Remove dynamic with only not deregulated genes
comparison_matrix = comparison_matrix[-nrow(comparison_matrix),]

## Alluvial plot
allu = aluvial_plot(comparison_matrix, fill_color = rep(c("grey","blue","red"),5))
ggsave(filename = "Kit_SC_SCII_RS_alluvial_plot_FDR.png", plot = allu, width = 15, height = 5,
  device = 'png', dpi = 450)
