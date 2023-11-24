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

dir = "SC"
path_RNAseq = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_SC_SCII_RS_spike/"
path_ChIPseq = "~/Documents/Mouse/Comparison_DEG_mark/"
path_RE = "~/Documents/Mouse/Enhancer_Annotation/ChromHMM_Kit_SC_RS_DEG_SC_SCII_RS/"

RE = read.table(paste0(path_RE,dir,"_DEG_H3K79me2_ENH_Met1_100kb_FDR.tsv"), sep = "\t", h = T)
DEG = read.table(paste0(path_RNAseq,"DOT1L_KO_effect_within_",dir,"_cells_volcano_table_FDR.tsv"), h = T,
  sep = "\t")[,c(1:5,14)]
GB = read.table(paste0(path_ChIPseq,dir,"/SC_SCII_RS/",dir,"_gene_enrich_H3K79me2_SC_SCII_RS_FDR.tsv"),
  sep = "\t", h = T)
# GB = read.table(paste0(path_ChIPseq,"/Kit_m/Kit/",dir,"_DEG_H3K79me2_GB.tsv"),
#   sep = "\t", h = T) ## GSC

####################################################################################################
##                                          INITIALIZATION
####################################################################################################

GB$GB = 1
RE_GB = unique(merge(DEG,GB, by.x = "gene_ENS", by.y = "geneId", all = T)[,c(1:6,16)])
RE_GB[is.na(RE_GB$GB),"GB"] = 0

####################################################################################################
##                                          Methode 1
####################################################################################################

RE$RE = 0
RE[RE$H3K79me2 == "H3K79me2","RE"] = 1
RE_GB = unique(merge(RE_GB,RE, by.x = "gene_ENS", by.y = "gene_name", all = T)[,c(1:7,16)])
RE_GB[is.na(RE_GB$RE),"RE"] = 0

RE_GB$type = RE_GB$RE - RE_GB$GB
# type = 0 RE-GB | NO
# type = -1 GB
# type = 1 RE
pos = rownames(RE_GB[RE_GB$GB == 0 & RE_GB$RE == 0,])
nbline = nrow(RE_GB)
RE_GB[pos,"type"] = -2
RE_GB = RE_GB[1:nbline,]
RE_GB$type = factor(RE_GB$type, levels = c(-1, 1, 0, -2), labels = c("GB","RE","RE-GB", "NO"))

write.table(RE_GB,paste0(dir,"_distribution_RE_GB_DEG_FDR.tsv"), sep = "\t", col.names = T,
  row.names = F, quote = F)

p = ggplot(RE_GB, aes(x = gl.x, fill = type)) +
  geom_bar(stat = "count", position = "fill") +
  xlab("") + ylab("") + labs(fill = "") +
  scale_fill_manual(values = brewer.pal(10,"Set2")) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 35, angle = 0, hjust = 1),
    legend.position = "bottom")

ggsave(filename = paste0(dir,"_distribution_RE_GB_DEG_FDR.png"),plot = p, width = 12,
  height = 12,device = 'png',dpi = 300)
