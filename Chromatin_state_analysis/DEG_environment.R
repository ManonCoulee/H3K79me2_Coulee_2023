####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('devtools', 'ggplot2', 'GenomicRanges', 'parallel', 'stringr','ChromENVEE')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                        INITIALIZATION
####################################################################################################

## Path of data
wd = "~/Documents/Mouse/Comparison_H3K79me2_enhancers_state"
chromHMM = "~/Documents/Mouse/ChIPseq/ChromHMM"
path_RNAseq = "~/Documents/Mouse/RNAseq/Spike_analysis_Kit_SC_SCII_RS/"
path_H3K79me2 = "~/Documents/Mouse/ChIPseq/H3K79me2/"

genomeFile = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed", sep = "\t", h = F)
colnames(genomeFile) = c("seqnames","start","end","strand","score","gene_ENS")
genomeFile$seqnames = as.character(genomeFile$seqnames)

col = getStateColor(colorTable = colorTable)

list_dirs = c("Kit-","Kit+","SC","RS")

## Theme for ggplot
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

####################################################################################################
##                               Chromatin state localisation
####################################################################################################

## Load of all segments file
list_file_chromatin_state = paste0(chromHMM,"/model_18_makeSegmentation_AGSC_SC_RS/",list_dirs,
  "_18_segments.bed")
names(list_file_chromatin_state) = list_dirs

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
names(list_data_chromatin_state) = list_dirs

####################################################################################################
##                            Chromatin state gene environment
####################################################################################################

dir = "RS"
stateOrderReduce = c("TSSA","TSSFlnk","TSSFlnk","Tx","Tx","EnhG","EnhG","EnhA","EnhWk","ZNF.Rpts",
  "Het","TssBiv","EnhBiv","ReprPC","ReprPC","Quies","Quies","Quies")

## Return chromatin state for cell type
table_chromHMM = list_data_chromatin_state[[dir]]
table_chromHMM$stateName = factor(table_chromHMM$state, levels = colorTable$stateNumber,
    labels = stateOrderReduce)

## Associate gene to gene localisation
table_DEG = read.table(paste0(path_RNAseq,dir,"_spike_volcano_table_FDR.tsv"),
  sep= "\t", h = T)[,-c(2:6)]
table_DEG_bed = merge(genomeFile, table_DEG, by = "gene_ENS")

## Define the coverage of chromatin state in gene environment
table_overlapping = geneEnvironment(table_DEG_bed, table_chromHMM, unique(stateOrderReduce))
rownames(table_overlapping) = table_overlapping$gene_ENS
table_overlapping$gl_reduce = factor(table_overlapping$gl, levels = unique(table_overlapping$gl),
  labels = c(1, 0, 2))
table_overlapping$gl_reduce = as.numeric(table_overlapping$gl_reduce)

## Define predominant chromatin state
result_umap = predominantState(table_overlapping, state = unique(stateOrderReduce),
  header = unique(stateOrderReduce) ,neighbors = 32, metric = "euclidean", dist = 0.5)
result_umap$gl = table_overlapping$gl
table = cbind(table_overlapping,result_umap[,c("state")])
write.table(table, paste0(dir,"_FDR_spike.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")

####################################################################################################
##                       Quality controle to define predominant chromatin state
####################################################################################################

colnames(table_overlapping)[32] = "state"
rownames(table_overlapping) = table_overlapping$gene_ENS
stateOrderReduce = c("TSSA","TSSFlnk","Tx","EnhG","EnhA","EnhWk","ZNF.Rpts","Het","TssBiv",
  "EnhBiv","ReprPC","Quies")


## Distribution of max coverage percentage
table_overlapping$max_pct = apply(table_overlapping[,19:30],1,max)

## Second max chromatin state
# table_overlapping$scd_pct = unlist(lapply(rownames(table_overlapping), function(gene) {
#   name = stateOrderReduce[order(table_overlapping[gene,stateOrderReduce], decreasing = T)[2]]
#   return(table_overlapping[gene,name])
# }))

## Representation
RSp1 = ggplot(table_overlapping, aes(y = sort(max_pct), x = 1:nrow(table_overlapping))) +
  geom_hline(aes(yintercept = median(max_pct)), color = "black", linetype = "dashed") +
  geom_text(aes(100,median(max_pct), label = round(median(max_pct),2), vjust = -1))+
  # geom_hline(aes(yintercept = quantile(max_pct)[2]), color = "black", linetype = "dashed") +
  # geom_text(aes(100,quantile(max_pct)[2], label = round(quantile(max_pct)[2],2), vjust = -1))+
  geom_hline(aes(yintercept = median(scd_pct)), color = "red", linetype = "dashed") +
  geom_text(aes(100,median(scd_pct), label = round(median(scd_pct),2), vjust = 2))+
  # geom_hline(aes(yintercept = quantile(scd_pct)[2]), color = "red", linetype = "dashed") +
  # geom_text(aes(100,quantile(scd_pct)[2], label = round(quantile(scd_pct)[2],2), vjust = -1))+
  geom_smooth(color = "black") +
  geom_smooth(aes(y = sort(scd_pct), x = 1:nrow(table_overlapping)), color = "red") +
  xlab("") + ylab("% predominant state") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=25, angle = 0),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 25, angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    legend.position = "bottom")
ggsave(filename = paste0(dir,"_max_coverage_distribution.png"), plot = RSp1, device = "png",
    width = 6, height = 6, dpi = 600)

####################################################################################################
##                            Distribution of predominant state
####################################################################################################

colnames(table_overlapping)[32] = "state"
rownames(table_overlapping) = table_overlapping$gene_ENS
stateOrderReduce = c("TSSA","TSSFlnk","Tx","EnhG","EnhA","EnhWk","ZNF.Rpts","Het","TssBiv",
  "EnhBiv","ReprPC","Quies")

## Barplot de distributon des états chromatinens
table(table_overlapping$gl)
x1 = data.frame(table(table_overlapping$state, table_overlapping$gl))
x1$percentage = (x1$Freq/(rep(table(table_overlapping$gl),
  each = length(unique(table_overlapping$state)))))*100
x2 = data.frame(table(table_overlapping$state))
x2$Var2 = "all genes"
x2$percentage = (x2$Freq/nrow(table_overlapping))*100
x = rbind(x1,x2)

write.table(x,paste0(dir,"_environment_percentage_spike.tsv"), quote = F, col.names = F,
  row.names = F, sep = "\t")

df = data.frame(table(table_overlapping$state,table_overlapping$gl))
df2 = data.frame(table(table_overlapping$state))
df2$Var2 = "ALL"
df_all = rbind(df,df2)

p = ggplot(x1, aes(x = Var2, y = Freq,
    fill = factor(Var1, level = unique(stateOrderReduce)))) +
  geom_bar(stat = "identity", position = position_fill(), size = 0.1) +
  scale_fill_manual(values = col$stateName) +
  scale_x_discrete(labels = c("Down","Not","Up")) +
  xlab("") + ylab("") + labs(fill = "", alpha = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 25, angle = 0, hjust = 1),
    axis.text.x = element_text(colour = c("blue","black","red")),
    legend.position = "none")
ggsave(filename = paste0(dir,"_environment_percentage_barplot_FDR.png"), plot = p, device = "png",
  width = 4, height = 4, dpi = 300)
