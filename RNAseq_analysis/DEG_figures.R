####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('UpSetR','ggplot2', 'gtable','stringr','RColorBrewer', 'edgeR', 'DESeq2',
  'VennDiagram','gridExtra','pheatmap','factoextra','ggalluvial') #'clusterProfiler'
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
  pattern = "cells_volcano_table_FDR_PV.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})

gene_length = read.table("~/Documents/Annotations/mouse_gene_length.tsv",stringsAsFactors = F,
  sep = "\t", h = F)
annotation = read.table("~/Documents/Annotations/mouse_gene_annotation.tsv",stringsAsFactors = F,
  h = T, sep = '\t')
bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed",h=F, sep = '\t')
colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")


####################################################################################################
##                                  Distribution DEG barplot
####################################################################################################

gene = unique(c(list_data$SC$gene_ENS, list_data$RS$gene_ENS, list_data$SCII$gene_ENS,
  list_data$'Kit-'$gene_ENS, list_data$'Kit+'$gene_ENS))
# gene = unique(c(list_data$SC$gene_ENS, list_data$RS$gene_ENS, list_data$SCII$gene_ENS))

RS = list_data$RS[list_data$RS$PValue <= 0.05,]
SC = list_data$SC[list_data$SC$PValue <= 0.05,]
SCII = list_data$SCII[list_data$SCII$PValue <= 0.05,]
SSC_m = list_data$'Kit-'[list_data$'Kit-'$PValue <= 0.05,]
SSC_p = list_data$'Kit+'[list_data$'Kit+'$PValue <= 0.05,]

df = rbind(data.frame(table(SSC_m$gl)),data.frame(table(SSC_p$gl)),data.frame(table(SC$gl)),
  data.frame(table(SCII$gl)),data.frame(table(RS$gl)))
df$name = rep(samples[c(1:2,4,5,3)],each = 2)

## Distribution of DEG genes (barplot)
p = ggplot(df, aes(x = Var1, y = Freq, alpha = factor(name, levels = samples[c(1:2,4,5,3)]))) +
  geom_bar(stat = "identity", position = "dodge", fill = rep(c("blue","red"),5)) +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "bottom")
  # scale_fill_manual(values = c("blue","red"))
ggsave(filename = "DEG_distribution_barplot.png", plot = p,width = 10, height = 8,
  device = 'png', dpi = 150)

p2 = ggplot(df, aes(x = factor(name, levels = samples[c(1:2,4,5,3)]), y = Freq,fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "bottom") +
  scale_fill_manual(values = c("blue","red"))
ggsave(filename = "DEG_distribution_barplot_V2.png", plot = p2,width = 10, height = 8,
  device = 'png', dpi = 150)

####################################################################################################
##                                      Centromere distribution
####################################################################################################

RS = read.table("DOT1L_KO_effect_within_RS_cells_volcano_table_FDR.tsv",sep = "\t", h = T)
RS_anno = RS[RS$gl != "Not regulated",c(2:4,14)]
RS_anno$gl = factor(RS_anno$gl, labels = c("DR","UR"))
chrom = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", sep = ";", h = T)[,1:2]
chrom$Chromosome = paste0("chr",chrom$Chromosome)

centro = read.table("~/Documents/Annotations/mouse_centromere_location.tsv",sep = "\t", h = F)[,c(2:4,8)]
colnames(centro) = c("chr","start","end","gl")
RS_anno = rbind(RS_anno,centro)
RS_anno = merge(RS_anno, chrom, by.x = "chr", by.y = "Chromosome", all.x = T)

p = ggplot(RS_anno) +
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

ggsave('RS_DEG_gene_position_chr.png',  plot = p, width = 10, height = 8, device = 'png', dpi = 600)

