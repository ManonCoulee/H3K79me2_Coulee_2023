####################################################################################################
##                                       LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'reshape', 'TCseq', 'statmod',
  'GenomicFeatures','RUVseq')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                    INITIALIZATION
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
##                                   Gene expression count
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

gene_annotation = data.frame(do.call(rbind, strsplit(x = gsub(
  pattern = "^>", replacement = '', x = unlist(system(paste('zgrep -P "^>"',
    file.path('~/Documents/Annotations/gencode.vM19.transcripts.fa.gz')), intern = T)),
  perl = T), split = '|', fixed = T)), stringsAsFactors = F)
names(gene_annotation) = c('tx_ENS', 'gene_ENS', 'gene_OTT', 'tx_OTT', 'tx_name', 'gene_name',
  'tx_length', 'gene_type')
gene_annotation[gene_annotation == '-'] = NA
gene_annotation$tx_length = as.numeric(gene_annotation$tx_length)

## Estimating gene_length by using GenomicFeatures
txdb = makeTxDbFromGFF(file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'),
  format = 'gtf')
exons.list.per.gene = exonsBy(x = txdb, by = 'gene')
## Same but much faster
exonic.gene.sizes = sum(width(reduce(exons.list.per.gene)))
genecounts$rpkm = rpkm(y = x3, gene.length = exonic.gene.sizes[match(names(exonic.gene.sizes),
  rownames(x3))], normalized.lib.sizes = T, log = F)[expressed_genes, ]

  m_export = unique(subset(gene_annotation[, -grep('tx_', x = names(gene_annotation))],
    select = -gene_type))
  rownames(m_export) = m_export$gene_ENS

####################################################################################################
##                         DE ANALYSIS with spike-in normalization (RUV method)
####################################################################################################

genetable = genecounts$raw

DEG_spike = lapply(unique(SamplePlan$CellType), function(cell) {
  ## Rescale data using RUVg method
  path = SamplePlan[SamplePlan$CellType == cell,]
  path$sample = paste0(path$SamplePath,"/",rownames(path))
  filtered = genetable[,path$sample]
  genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
  spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
  x <- as.factor(path$SampleType)
  set <- newSeqExpressionSet(as.matrix(filtered),
    phenoData = data.frame(x, row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")
  set1 <- RUVg(set, spikes, k=1)

  ## Make model design
  design <- model.matrix(~x + W_1, data=pData(set1))

  ## Process differential gene analysis
  y <- DGEList(counts=counts(set1), group=x)
  y <- calcNormFactors(y)
  # y <- estimateDisp(y = y, design = design)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)

  fit <- glmFit(y, design)
  qlf = glmTreat(glmfit = fit, coef = "xKO", lfc = log(1.5))
  tt = with(topTags(object = qlf, n = NULL, sort.by = 'none'), table)
  tt$gl = 0
  tt$gl[tt$FDR <= 0.05  & tt$logFC < -log(1.5)] = 1
  tt$gl[tt$FDR <= 0.05  & tt$logFC > log(1.5)] = 2
  tt$gl[grepl("ERCC",rownames(tt))] = 3
  tt$gl = factor(tt$gl, levels = c(1, 0, 2, 3), labels = c('Down-regulated in SCII', 'Not regulated',
    'Up-regulated in SCII',"ERCC"))
  tt$cell = cell
  tt = tt[tt$gl != "ERCC",]
  tt$gene_ENS = rownames(tt)

  ## Make a MD plot
  p = ggplot(data = tt) +
    geom_point(aes(x = logFC, y = -log10(PValue), color = gl), size = 5, alpha = .5) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), alpha = .5,
      color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), alpha = .5, color = 'grey50',
      linetype = 'dashed') +
    # scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Paired')[c(2,1)],
    #   brewer.pal(9,'Set1')[c(9)],brewer.pal(9,'Paired')[c(5,6)])) +
    scale_color_manual(name = 'Gene expression', values = c(brewer.pal(9,'Set1')[c(2,9,1)])) +
    labs(title = paste(paste0('FC>',1.5),paste0('& FDR<',100*0.05,'%'),
        'V plot'),
      x = 'Log Fold-change',y = '-Log10(P)') +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
  ggsave(file.path(rout, paste0(gsub(pattern = ' ',replacement = '_', x = cell),'_volcano_plot_FDR.png')),
    plot = p, width = 10, height = 8, device = 'png', dpi = 300)

  write.table(file.path(rout,paste0(gsub(pattern = ' ', replacement = '_', x = cell),
      '_volcano_table_FDR.tsv')),
      x = data.frame(m_export[rownames(tt),],tt, check.names = F), sep = "\t",
      row.names = F, quote = F)
  return(tt)
})

####################################################################################################
##                                DE ANALYSIS count barplot figure
####################################################################################################

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

ggsave('distribution_barplot_spike-in.png', plot = p, width = 10, height = 6, device = 'png', dpi = 300)
