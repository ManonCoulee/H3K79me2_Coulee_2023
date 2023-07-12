####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'FactoMineR', 'reshape', 'TCseq', 'statmod',
  'GenomicFeatures','RUVseq')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                               INITIALIZATION
####################################################################################################

wd = '.'
sample_dir = "~/Documents/Mouse/RNAseq/Samples"
dirs = list.dirs(sample_dir, full.names = T, recursive = F)
files = list.files(dirs,"RPG")

rout = wd

LFC = log2(1.5)

CT = list.dirs(sample_dir, full.names = F, recursive = F)

## Load SamplePlan data descriptor
SamplePlan = read.table(paste0(wd,"/SamplePlan.tsv"),sep = "\t", h = T, row.names = 1)
SamplePlan$SampleType = factor(SamplePlan$SampleType)
SamplePlan$CellType = factor(SamplePlan$CellType,levels = CT)
SamplePlan$SamplePool = factor(SamplePlan$SamplePool)

####################################################################################################
##                                             GENE EXPRESSED
####################################################################################################

x = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)), function(path){
  read.table(paste0(path, '.RPG.tsv'), sep = "\t", header = F, skip = 4)[,4]
}, simplify = T),
  row.names = read.table(paste0(SamplePlan[1,"SamplePath"],"/",rownames(SamplePlan)[1],".RPG.tsv"),
  sep = "\t", header = F, skip = 4)[,1])

x2 = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)), function(path){
  read.table(paste0(path, '.STAR.idxstats'), sep = "\t", header = F, skip = 0)[,3]
}, simplify = T),
  row.names = read.table(paste0(SamplePlan[1,"SamplePath"],"/",rownames(SamplePlan)[1],".STAR.idxstats"),
  sep = "\t", header = F, skip = 0)[,1])
x2 = x2[grepl("ERCC", rownames(x2)),]

x3 = rbind(x,x2)

## Gene counts data normalization
min_sample = 2
min_expression = 1
xx = cpm(y = x3, normalized.lib.sizes = T, log = F)
## Filter-out low expressed genes (cpm > 1 in minimum 2 samples)
expressed_genes = rowSums(xx > min_expression) >= min_sample

genecounts = list(
  ## Raw read counts (as estimated by STAR)
  raw = x3[expressed_genes, ],
  ## For expression data visualization
  cpm = cpm(y = x3[expressed_genes, ], lib.size = colSums(x3), normalized.lib.sizes = T, log = F),
  rlog = assay(rlog(DESeqDataSetFromMatrix(countData = x3,
    colData = SamplePlan, design = ~1)))[expressed_genes, ]
)

####################################################################################################
##                                   DE ANALYSIS without spike-in scaling
####################################################################################################

cell = factor(SamplePlan$CellType)
phenotype = relevel(SamplePlan$SampleType, ref = "CTL")

## Create a model with all multifactor combinations
design = model.matrix(~phenotype+cell+phenotype:cell)
colnames(design) = c("all","KOvsCTL.Kitm","CTL.Kitp","CTL.RS","CTL.SC","CTL.SCII","KOvsCTL.Kitp",
  "KOvsCTL.RS","KOvsCTL.SC","KOvsCTL.SCII")

####################################################################################################
##                         DE ANALYSIS with spike-in normalization (RUV method)
####################################################################################################

genetable = genecounts$raw

DEG_spike = lapply(unique(SamplePlan$CellType), function(cell) {
  path = SamplePlan[SamplePlan$CellType == cell,]
  path$sample = paste0(path$SamplePath,"/",rownames(path))
  zfGenes = genetable[,path$sample]
  filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
  filtered <- zfGenes[filter,]
  genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
  spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
  x <- as.factor(path$SampleType)
  set <- newSeqExpressionSet(as.matrix(filtered),
    phenoData = data.frame(x, row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")
  set1 <- RUVg(set, spikes, k=1)
  design <- model.matrix(~x + W_1, data=pData(set1))
  y <- DGEList(counts=counts(set1), group=x)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)

  fit <- glmFit(y, design)
  qlf = glmTreat(glmfit = fit, coef = "xKO", lfc = log(1.5))
  # lrt <- glmLRT(fit)
  tt = with(topTags(object = qlf, n = NULL, sort.by = 'none'), table)

  tt$gl = 0
  tt$gl[tt$FDR <= 0.05  & tt$logFC < -log(1.5)] = 1
  tt$gl[tt$FDR <= 0.05  & tt$logFC > log(1.5)] = 2
  tt$gl[grepl("ERCC",rownames(tt))] = 3
  tt$gl = factor(tt$gl, levels = c(1, 0, 2, 3), labels = c('Down-regulated', 'Not regulated',
    'Up-regulated',"ERCC"))
  tt$cell = cell
  tt = tt[tt$gl != "ERCC",]
  write.table(tt,paste0("DOT1L_KO_effect_within_",cell,"_cells_volcano_table_FDR.tsv"),
    row.names = F, col.names = T, sep = "\t", quote = F)

  return(tt)
})

df = do.call("rbind",DEG_spike)

p = ggplot(df[df$gl != "Not regulated",],
  aes(x = factor(gl),
    alpha = factor(cell, levels = c("Kit_m_spike","Kit_p_spike","SC_spike","SCII_spike","RS_spike"),
    labels = c("KIT-","KIT+","SC","SCII","RS")),
    fill = gl)) +
  geom_bar(stat = "count", position = "dodge") +
  xlab("") + ylab("number of genes") + labs(fill = "", alpha = "") +
  scale_fill_manual(values = rep(c("blue","red"),5)) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "bottom") +
  guides(fill = "none")

ggsave('distribution_barplot_spike-in_RUV.png', plot = p, width = 10, height = 6, device = 'png', dpi = 300)
