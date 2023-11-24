####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'UpSetR', 'ChIPseeker', 'clusterProfiler',
  'org.Mm.eg.db', 'GenomicFeatures','ggalluvial')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

load("~/Documents/Annotations/MouseTableBed")

####################################################################################################
##                                          INITIALIZATION
####################################################################################################

SamplePlan = read.table("SamplePlan.tsv", h = T, sep = "\t", stringsAsFactors = F)
SamplePlan$CellType = factor(SamplePlan$CellType, levels = SamplePlan$CellType)

list_peaks = lapply(SamplePlan$SamplePath[c(1:3)], function(name) {
  read.csv(name, sep = "\t", h = T)[,c(1:3,6,8)]
})

names(list_peaks) = SamplePlan$SampleID[c(1:3)]

## Create an empty table wit peaks information and sample
df = do.call("rbind",lapply(names(list_peaks), function(id) {
  print(id)
  peak = read.csv(SamplePlan[SamplePlan$SampleID == id,"SamplePath"], sep = "\t", h = T)[,c(1:3,6,8)]
  peak$grp = id
  return(peak)
}))
write.table(df, "../Comparison_DEG_mark/Switch_enrichment/KIT_SC_RS_peak.bed",sep = "\t", col.names = F,
  row.names = F, quote = F)

## Open file contains peaks which are merge together
open_region = read.table("Kit_SC_RS_merge_peaks_pos.bed", sep = "\t", h = F)

####################################################################################################
##                                   QUALITY CONTROL (supplemental)
####################################################################################################

## Size of merge peak
open_region$size = abs(open_region$V3-open_region$V2)
pct90common = quantile(open_region$size, probs = seq(0,1,0.1))[10]

## Size of peak
df$size = abs(df$end-df$start)
pct90 = quantile(df$size, probs = seq(0,1,0.1))[10]

## Plot merge peak with 9th percentile value
p = ggplot(open_region, aes(x = size)) +
  geom_histogram() +
  geom_vline(aes(xintercept = pct90, color = "9th quantile H3K79me2 size")) +
  geom_vline(aes(xintercept = pct90common, color = "9th quantile common H3K79me2 size")) +
  scale_fill_manual(values = "grey") +
  labs(color = "") + xlab("peak size (bp)") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=25, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right")

ggsave(filename = "Kit_SC_RS_common_regions_distribution_plot.png", plot = p, width = 8,
  height = 5, device = 'png', dpi = 450)

####################################################################################################
##                                      Dynamic analysis alluvial
####################################################################################################

## Return position of presence peak for each cell type
x = lapply(SamplePlan$SamplePath[c(1:3)], function(name) {
  peak = read.csv(name, sep = "\t", h = T)[,1:4]

  presence = unique(unlist(lapply(rownames(peak),function(x){
    chrP = peak[x,"seqnames"]
    SP = peak[x,"start"]
    EP = peak[x,"end"]

    ## Return only chr matching table
    bin_table = open_region[open_region$V1 == chrP,]

    sub_bin_table = bin_table[bin_table$V3 >= SP,]
    pos = sub_bin_table[sub_bin_table$V2 <= SP,]
    return(rownames(pos))

  })))
  return(presence)
})

## Associate each present position in the table according to cell type
names(x) = SamplePlan$SampleID[c(1:3)]
open_region[,SamplePlan$SampleID[c(1:3)]] = 0
open_region[as.numeric(x[[1]]),names(x)[1]] = 1
open_region[as.numeric(x[[2]]),names(x)[2]] = 1
open_region[as.numeric(x[[3]]),names(x)[3]] = 1

## Alluvial of present region
open_region_upset = open_region
open_region_upset[open_region_upset == 0] = "absent"
open_region_upset[open_region_upset == 1] = "present"

# Create a comparison matrix
comparison_matrix = expand.grid(GS = c("absent","present"),
  SC = c("absent","present"),
  RS = c("absent","present"))
comparison_matrix$motif = paste0(comparison_matrix$GS,
  comparison_matrix$SC,
  comparison_matrix$RS)
rownames(comparison_matrix) = comparison_matrix$motif
comparison_matrix$freq = 0
open_region_upset$motif = unlist(lapply(rownames(open_region_upset), function(line) {
  l = open_region_upset[line,4:6]
  l = paste0(l,collapse = "")
  return(l)
}))
write.table(open_region_upset, "Kit_SC_RS_open_region.tsv", col.names = T, row.names = F,
  sep  ="\t", quote = F)

## Count the frequence of motif
x = table(open_region_upset$motif)
comparison_matrix[names(x),"freq"] = x

## Create alluvial plot
allu = ggplot(comparison_matrix[-c(1,6),], aes(y = freq, axis1 = factor(GS),
  axis2 = factor(SC),
  axis3 = factor(RS),
  fill = factor(motif))) +
  geom_alluvium() +
  geom_stratum(width = 1/36, fill = rep(c("black",NA),3), color = rep(c("black",NA),3)) +
  scale_x_discrete(limits = c("GSC","SCI", "RS"), expand = c(.1, .1)) +
  scale_fill_manual(values = brewer.pal(9, 'Set1'),
    labels= c("RS spe","SCI spe","SCI-RS","GSC spe","GSC-SCI", "common")) +
  ylab("") + labs(fill = "") +
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
    legend.position = "bottom")
ggsave(filename = "GS_SC_RS_alluvial_plot_V3.png",plot = allu, width = 6, height = 3,
  device = 'png', dpi = 600)

####################################################################################################
##                                   Annotation of dynamic pattern
####################################################################################################

## Initialization
txmm = makeTxDbFromGFF(file = file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'),
  format = 'gtf')
colnames(open_region_upset) = c("chr","start","end","GS","SC","RS","motif")

## ChIPseeker annotation
tt_peak = makeGRangesFromDataFrame(open_region_upset, keep.extra.columns = TRUE)
tt_annotate = annotatePeak(tt_peak, tssRegion = c(-3000,3000), TxDb = txmm, annoDb = "org.Mm.eg.db")
tt_anno = tt_annotate@anno
tt_anno$gene_ENS = unlist(strsplit(tt_anno$geneId,".", fixed = T))[seq(1,length(tt_anno)*2,2)]
tt_anno$gene_name = lapply(tt_anno$gene_ENS, function(gene) {
  unique(MouseTableBed[MouseTableBed$ensembl_gene_id == gene,"mgi_symbol"])
})

## Reduce complexity of annotation
tt_reduce = as.data.frame(tt_anno)
tt_reduce$reduce = "Intragenic"
tt_reduce[grepl("Downstream",tt_reduce$annotation),"reduce"] = "Downstream"
tt_reduce[grepl("Promoter",tt_reduce$annotation),"reduce"] = "Promoter"
tt_reduce[grepl("Distal",tt_reduce$annotation),"reduce"] = "Distal Intergenic"
write.table(tt_reduce, "Kit_SC_RS_peaks_annotated.tsv",col.names = T, row.names = F, sep = "\t",
  quote = F)

## Barplot of annotation categories
ll = lapply(names(x)[-5], function(x) { ## Remove GSC-RS motif
  return(tt_reduce[tt_reduce$motif == x,])
})
df = do.call("rbind",ll)

p = ggplot(df, aes(x = factor(motif, levels = names(table(df$motif))[c(4,5,2,3,1,6)]), fill = reduce)) +
  geom_bar(stat = "count", position = position_fill()) +
  scale_fill_manual(values = c("#B2BABB","#A020F0","#E67E22","#2980B9")) +
  ylab("") + xlab("") +
  scale_x_discrete(labels= c("GSC spe","GSC-SCI","SCI spe","SCI-RS","RS spe","common")) +
  theme_minimal() +
      theme(text = element_text(size = 40),legend.position = "bottom",legend.title = element_blank()) +
      guides(fill = guide_legend(nrow=2, byrow = T))
ggsave(filename = "GS_SC_RS_annotation.png",plot = p, width = 15, height = 5, device = 'png',
  dpi = 600)

####################################################################################################
##                                  GO enrichment of dynamic pattern
####################################################################################################

## GO enrichment
list_table = lapply(names(x), function(x) {
  return(tt_reduce[tt_reduce$motif == x,"ENTREZID"])
})
names(list_table) = c("RS spe","SCI spe","SCI-RS","GSC spe","GSC-RS","GSC-SCI","common")
compareGO = compareCluster(geneCluster = list_table[c(4,6,2,3,1,7)], fun = "enrichGO",
  Org = "org.Mm.eg.db",ont = "BP", pAdjustMethod = "BH")

png("GS_SC_RS_GO_dotplot.png", height = 400, width = 800)
dotplot(compareGO, showCategory = 5)
dev.off()
