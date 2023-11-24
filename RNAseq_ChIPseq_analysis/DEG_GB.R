####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

path_RNAseq = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_SC_SCII_RS_spike/"
path_H3K79me2 = "~/Documents/Mouse/ChIPseq/H3K79me2/"

dir = "SC"

table = read.csv(paste0(path_H3K79me2,dir,"/",dir,"_peak_annotate_reduce.tsv"), sep = "\t", h = T)
table_DEG = read.table(paste0(path_RNAseq,"DOT1L_KO_effect_within_",dir,"_cells_volcano_table.tsv"),
  sep= "\t", h = T)

####################################################################################################
##                                      ENRICHMENT MARK barplot
####################################################################################################

x = merge(table,table_DEG,by.x = "geneId",by.y = "gene_ENS")[,c(1:4,21,25:28,38)]
colnames(x) = c("geneId","chr_peak","start_peak","end_peak","gene_symbol","annotation","chr_gene",
	"start_gene","end_gene","gl")

write.table(x, "distribution_DEG_H3K79me2_GB_FDR.tsv", sep = "\t", col.names = T, row.names = F,
  quote = F)

x = unique(x[,c(1,7:10)])

## Plot
df1 = data.frame(Freq = c(table(x$gl),table(table_DEG$gl)-table(x$gl)),
	Var1 = rep(c("DR","NR","UR"),2),
  grp = rep(c("+H3K79me2","-H3K79me2"), each = 5))
write.table(df1, paste0(dir,"_DEG_H3K79me2_GB_FDR.tsv"), sep = "\t", col.names = T, row.names = F,
  quote = F)

p = ggplot(df1, aes(x = Var1, y = Freq, fill = grp)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("") + ylab("") + labs(fill = "") +
  scale_fill_manual(values = c("grey80","grey50")) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=45, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 45, angle = 0, hjust = 1),
    legend.position = "none")

ggsave(paste0(dir,"_DEG_H3K79me2_GB_PV.png"), plot = p, width = 5, height = 5,
  device = "png", dpi = 600)
