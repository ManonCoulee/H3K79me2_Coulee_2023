####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'FactoMineR', 'reshape', 'TCseq', 'statmod',
  'GenomicFeatures', 'VennDiagram', 'pheatmap')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

wd = '.'
sample_dir = "~/Documents/Mouse/RNAseq/Samples"
dirs = list.dirs(sample_dir, full.names = T, recursive = F)
files = list.files(dirs,"RPG")

LFC = log2(1.5)

CT = list.dirs(sample_dir, full.names = F, recursive = F)[-c(3)]
nk = 8

## Load SamplePlan data descriptor
SamplePlan = read.table(paste0(wd,"/SamplePlan.tsv"),sep = "\t", h = T, row.names = 1)
SamplePlan$SampleType = factor(SamplePlan$SampleType)
SamplePlan$CellType = factor(SamplePlan$CellType,levels = CT)
SamplePlan$CellType = factor(SamplePlan$CellType)
SamplePlan$SamplePool = factor(SamplePlan$SamplePool)
SamplePlan$CellTypeRed = factor(SamplePlan$CellTypeRed)

####################################################################################################
##                                           ANNOTATION
####################################################################################################

## Gencode gene counts from STAR
x = as.data.frame(sapply(paste0(SamplePlan$SamplePath,"/",rownames(SamplePlan)), function(path){
  read.table(paste0(path, '.RPG.tsv'), sep = "\t", header = F, skip = 4)[,4]
}, simplify = T),
  row.names = read.table(paste0(SamplePlan[1,"SamplePath"],"/",rownames(SamplePlan)[1],".RPG.tsv"),
  sep = "\t", header = F, skip = 4)[,1])

write.table(x,paste0(rout,"/genecounts_raw.tsv"), quote  = F, sep = "\t", col.names = T, row.names = T)

## Gene counts data normalization
min_sample = 2
min_expression = 0
xx = cpm(y = x, normalized.lib.sizes = T, log = F)

write.table(xx,paste0(rout,"/genecounts_cpm.tsv"), quote  = F, sep = "\t", col.names = T, row.names = T)

## Filter-out low expressed genes (cpm > 1 in minimum 2 samples)
expressed_genes = rowSums(xx > min_expression) >= min_sample

genecounts = list(
  ## Raw read counts (as estimated by STAR)
  raw = x[expressed_genes, ],
  ## For expression data visualization
  cpm = cpm(y = x[expressed_genes, ], lib.size = colSums(x), normalized.lib.sizes = T, log = F),
  ## For performing PCA analysis
  rlog = assay(rlog(DESeqDataSetFromMatrix(countData = x,
    colData = SamplePlan, design = ~1)))[expressed_genes, ]
)

## Distribution read alignment in each cell type
read_distribution = genecounts$cpm
table_variance = data.frame()
distribution = lapply(CT, function(x) {
  print(x)
  read_distribution_sample = read_distribution[,grepl(x,colnames(read_distribution))]
  list_condition = lapply(c("KO","Ctl"), function(condition) {
    read_distribution_sample_condition = read_distribution_sample[,grepl(condition,
      colnames(read_distribution_sample), ignore.case = T)]
    read_distribution_all = unlist(as.list(read_distribution_sample_condition))
    read_distribution_variance = var(read_distribution_all)
    print(paste0(condition," : ",read_distribution_variance))

    return(read_distribution_all)
  })
  return(list_condition)
})

## Load gene annotations (from Genecode transcripts.fa)
gene_annotation = data.frame(do.call(rbind, strsplit(x = gsub(
  pattern = "^>", replacement = '', x = unlist(system(paste('zgrep -P "^>"',
    file.path('~/Documents/Annotations/gencode.vM19.transcripts.fa.gz')), intern = T)),
  perl = T), split = '|', fixed = T)), stringsAsFactors = F)
names(gene_annotation) = c('tx_ENS', 'gene_ENS', 'gene_OTT', 'tx_OTT', 'tx_name', 'gene_name',
  'tx_length', 'gene_type')
gene_annotation[gene_annotation == '-'] = NA
gene_annotation$tx_length = as.numeric(gene_annotation$tx_length)
bed_annotation = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed", sep = "\t",
  h = F)
colnames(bed_annotation) = c("chr","start","end","strand","score","gene_ENS")

gene_annotation = merge(bed_annotation,gene_annotation, by = "gene_ENS")

## Estimating gene_length by using GenomicFeatures
txdb = makeTxDbFromGFF(file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'),
  format = 'gtf')
exons.list.per.gene = exonsBy(x = txdb, by = 'gene')
## Same but much faster
exonic.gene.sizes = sum(width(reduce(exons.list.per.gene)))
genecounts$rpkm = rpkm(y = x, gene.length = exonic.gene.sizes[match(names(exonic.gene.sizes),
  rownames(x))], normalized.lib.sizes = T, log = F)[expressed_genes, ]

write.table(genecounts$rpkm,"genecounts_rpkm.tsv", quote  = F, sep = "\t", col.names = T,
  row.names = T)

####################################################################################################
##                                        FILES GENERATION
####################################################################################################

## Create file for next GSEA analysis
table = data.frame(genecounts$cpm)
colnames(table) = rownames(SamplePlan)
table$Name = rownames(genecounts$raw)
table$Description = NA

for(dir in unique(SamplePlan$CellTypeRed)) {
  name = rownames(SamplePlan[SamplePlan$CellTypeRed == dir,])
  print(name)

  x_CTL = table[,name[grepl("Ctl",name, ignore.case = T)]]
  x_CTL$gene_expression = rowSums(x_CTL)/length(x_CTL)
  x_CTL$gene_ENS = rownames(x_CTL)
  write.table(x_CTL,paste0(dir,"_CTL_gene_expression.tsv"),row.names = F, col.names = T, quote = F,
    sep = "\t")

  x_KO = table[,name[grepl("KO",name)]]
  x_KO$gene_expression = rowSums(x_KO)/length(x_KO)
  x_KO$gene_ENS = rownames(x_KO)
  write.table(x_KO,paste0(dir,"_KO_gene_expression.tsv"),row.names = F, col.names = T, quote = F,
    sep = "\t")

  xx = table[,c("Name","Description",name)]
  # xx = table[,c("Name","Description",name[c(2,1,3,4,6,5)])] # SPIN1

  write.table(xx,paste0(dir,"_genecounts.gct"),row.names = F, col.names = T, quote = F,
    sep = "\t")
}

## Create a file for TC analysis
colnames(genecounts$raw) = rownames(SamplePlan)
write.table(genecounts$raw,paste0(rout,"/genecounts_filtered_raw.tsv"),row.names = T, col.names = T, quote = F,
sep = "\t")

####################################################################################################
##                                      SAMPLE DISTRIBUTION
####################################################################################################

## Perform a PCA to get an overview of overall sample distribution
x = PCA(X = t(genecounts$rlog), scale.unit = TRUE, graph = FALSE, ncp = 3)
SamplePCA = list(
  map = x$ind$coord,
  var = as.numeric(format(x$eig[1:3, 2], digits = 2, nsmall = 2, trim = TRUE)),
  cor = x$var$cor
)
colnames(SamplePCA$map) = colnames(SamplePCA$cor) = names(SamplePCA$var) = paste('Comp.', 1:3)
write.table(file.path(rout,'PCA-map.tsv'),
  x = data.frame(SampleID = rownames(SamplePCA$map), SamplePCA$map, check.names = F),
  sep = "\t", row.names = F, quote = F)

m_export = unique(subset(gene_annotation[, grep('gene', x = names(gene_annotation))],
  select = -gene_type))
m_export = m_export[match(rownames(SamplePCA$cor), m_export$gene_ENS), ]
write.table(file.path(rout, 'PCA-cor.tsv'),x = data.frame(m_export, SamplePCA$cor, check.names = F),
  sep = "\t", row.names = F, quote = F)

## Plot Individuals Factor map (first three pairwize component comparison)
invisible(lapply(1:length((L = list(c(1, 2), c(2, 3), c(1, 3)))), function(i){
  m = as.data.frame(SamplePCA$map[, -setdiff(1:3,L[[i]])]); names(m)[1:2] = c('x', 'y')
  p = ggplot(data = m, aes(x = x, y = y, shape = SamplePlan$SampleType,color=SamplePlan$CellTypeRed)) +
    geom_point(size = 10, alpha = .5) +
    geom_vline(xintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
    geom_hline(yintercept = 0, alpha = .5, color = 'grey50', linetype = 'dashed') +
    geom_text(aes(x = x, y = y), label = SamplePlan$SamplePool, fontface = 'bold', size = 4,
      alpha = .5, hjust = .5, vjust = .5, show_guide = FALSE) +
    scale_color_manual(name = 'CellType', values = brewer.pal(9, 'Set1')) +
    scale_shape_discrete(name = 'SampleType') +
    labs(
      title = 'Individuals Factor Map (PCA)',
      subtitle = paste('Comp.', L[[i]][1], 'vs.', L[[i]][2],
        paste0('(rlog normalized read counts / cpm>',min_expression,' / n=',nrow(SamplePCA$cor),')')),
      x = paste('Comp.', L[[i]][1], paste0('(', SamplePCA$var[L[[i]][1]], '%)')),
      y = paste('Comp.', L[[i]][2], paste0('(', SamplePCA$var[L[[i]][2]], '%)'))
    ) +
    theme(
      plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
      plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
      legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5)),
      panel.border = element_rect(colour = 'grey50', fill = NA)
    )
  ggsave(file.path(rout, paste0('PCA', L[[i]][1], '-', L[[i]][2], '.png')), plot = p, width = 12,
    height = 8, device = 'png', dpi = 150)
}))


