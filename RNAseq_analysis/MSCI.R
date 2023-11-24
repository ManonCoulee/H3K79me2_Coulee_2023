####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2','stringr','RColorBrewer','gridExtra','pheatmap','factoextra','ggalluvial',
  'ggpubr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

dir = c("../DE_analysis_3_factors_SC_SCII_RS_spike","../DE_analysis_3_factors_Kit")
samples = c("KIT-","KIT+","RS","SC","SCII")
rout = file.path(".")

####################################################################################################
##                                                 DATA
####################################################################################################


list_files = list.files(c(dir),
  pattern = "cells_volcano_table.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})
names(list_data) = samples

gene_expression = read.table("Gan_2013_RPKM.csv",stringsAsFactors = F, sep = "\t", h = T)

annotation = read.table("~/Documents/Annotations/mouse_gene_annotation.tsv",stringsAsFactors = F,
  h = T, sep = '\t')
bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed",h=F, sep = '\t')
colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")

####################################################################################################
##                              DISTRIBUTION GENES DEREGULATED (heatmap)
####################################################################################################

create_distribution_genes_deregulate_plot = function(list_data,bed,type_deregulation,out, name) {
  table_count_ref = as.data.frame(table(bed$chr))
  colnames(table_count_ref) = c("chr","ref_count")
  table_count_ref = table_count_ref[table_count_ref$chr != "chrM",]

  chr_order = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
    "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
  cell_type_order = c("SC","SCII","RS")

  ## Merge bed file with each data file
  list_data_chr = lapply(names(list_data), function(x,list) {
    m = list[[x]][,c("gene_ENS","chr","gl","FDR")]
    m = m[m$FDR <= 0.05,]
    t = as.data.frame(table(m$chr,m$gl))
    t$Var1 = as.character(t$Var1)
    t = t[t$Var1 != "chrM",]
    t$cell_type = x
    t
  }, list = list_data)

  names(list_data_chr) = names(list_data)

  table_count_data_pct = lapply(list_data_chr, function(x, ref) {
    colnames(x) = c("chr","deregulate_level","count","cell_type")
    up = x[x$deregulate_level == type_deregulation,]
    up = merge(up,ref,by = "chr")
    up$freq = (as.numeric(up$count)*100)/up$ref_count
    return(up)
  }, ref = table_count_ref)

  table_pct = do.call("rbind",lapply(table_count_data_pct,"[",,))

  table_pct$cell_type = factor(table_pct$cell_type, levels = cell_type_order)
  table_pct$chr = factor(table_pct$chr, levels = chr_order)

  pal = ifelse(type_deregulation == "Up-regulated (FDR < 5%)","Reds","Blues")

  p = ggplot(table_pct, aes(x = cell_type, y = ordered(chr, levels = rev(chr_order)),
    fill = freq)) +
    geom_tile() +
    geom_text(aes(label = paste0(count,"/",ref_count))) +
    scale_fill_distiller(palette = pal, trans = "reverse") +
    ggtitle(paste0("Distribution ",type_deregulation," genes (p-value<0.05)")) +
    xlab("") + ylab("") +
    theme(
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)) +
    labs(fill = "% genes")
  ggsave(filename = paste0("heatmap_distribution_",name,"_genes.png"),
    plot = p,width = 10, height = 8, device = 'png', dpi = 150)
}

create_distribution_genes_deregulate_plot(list_data, bed, "Up-regulated (FDR < 5%)", rout,
  "UR_FDR")
create_distribution_genes_deregulate_plot(list_data, bed, "Down-regulated (FDR < 5%)", rout,
  "DR_FDR")


####################################################################################################
##                              Gene expression distribution
####################################################################################################


create_genes_deregulate_FC_plot = function(name,out, listData) {

  data = listData[[name]]
  data = data[data$PValue <= 0.05,c("chr","logFC","PValue","FDR","gl")]
  # data = data[data$gl != "Not regulated",c("chr","logFC","PValue","FDR","gl")]
  data$chromosome = "autosome"
  data[data$chr == c("chrX","chrY"),"chromosome"] = "XY"
  XY = data[data$chromosome == "XY",]
  autosome = data[data$chromosome == "autosome",]

  DR_AUT = autosome$logFC
  DR_XY = XY$logFC

  p = ggplot(data, aes(x = chromosome, y = logFC)) +
    geom_boxplot() +
    labs(color = "") + xlab("") + ylab("log2FC KOvsCTL") +
    geom_hline(yintercept = 1.5, linetype = "dashed") +
    geom_hline(yintercept = -1.5, linetype = "dashed") +
    theme_bw() + theme(
      text = element_text(size=30, angle = 0),
      legend.position = "right",
      axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank()) +
    stat_compare_means(comparisons = list(c("autosome","XY")), method = "wilcox.test")

  ggsave(filename = paste0(name,"_FC.png"), plot = p, width = 5, height = 5, device = 'png',
    dpi = 600)

  data$cellType = name
  return(data)
}

listFC = lapply(samples, create_genes_deregulate_FC_plot, listData = list_data)
tableFC = do.call('rbind',listFC)
tableFC$cellType = factor(tableFC$cellType, levels = c("KIT-","KIT+","SC","SCII","RS"))

p = ggplot(tableFC, aes(x = chromosome, y = logFC, fill = cellType)) +
  geom_boxplot() +
  # geom_boxplot(aes(facet.by = chromosome)) +
  labs(fill = "") + xlab("") + ylab("log2FC KOvsCTL") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = -1.5, linetype = "dashed") +
  scale_fill_manual(values = c("grey10","grey30","grey50","grey70","grey90")) +
  # scale_fill_manual(values = c("grey30","grey60","grey90")) +
  # scale_fill_manual(values = c("DodgerBlue3","orange")) +
  theme_bw() + theme(
    text = element_text(size=30, angle = 0),
    legend.position = "right",
    # legend.position = "bottom",
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank())
    #+ stat_compare_means(label = "p.format", method = "wilcox.test")

ggsave(filename = "KIT_SC_SCII_RS_FC_1.png", plot = p, width = 6, height = 5, device = 'png',
  dpi = 600)

list_data = list(read.table("./KO_effect_volcano_table_without_filtration.tsv",sep = "\t", h = T))
names(list_data) = "RS"
d = create_genes_deregulate_FC_plot("RS",listData = list_data)

p = ggplot(d, aes(x = chromosome, y = logFC, fill = chromosome)) +
  geom_boxplot() +
  # geom_boxplot(aes(facet.by = chromosome)) +
  labs(fill = "") + xlab("") + ylab("log2FC KOvsCTL") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = -1.5, linetype = "dashed") +
  # scale_fill_manual(values = c("grey10","grey30","grey50","grey70","grey90")) +
  # scale_fill_manual(values = c("grey30","grey60","grey90")) +
  scale_fill_manual(values = c("DodgerBlue3","orange")) +
  theme_bw() + theme(
    text = element_text(size=30, angle = 0),
    # legend.position = "right",
    legend.position = "bottom",
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank()) +
    stat_compare_means(label = "p.format", method = "wilcox.test")

ggsave(filename = "RS_FC.png", plot = p, width = 6, height = 5, device = 'png',
  dpi = 600)

####################################################################################################
##                              Ratio XY/autosome in histone mark
####################################################################################################

library(ggplot2)
library(ggpubr)
length_chromosome = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", h = T, sep = ';',
  stringsAsFactors = FALSE)
mark_list = c("H3K4me1","H3K4me3","H3K9me3","H3K27ac","H3K27me3","H3K79me2")
list_files = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
  full.names = T)[c(2:4,6:8,10,13,14,16:24)]
list_data = data.frame(ratio = unlist(lapply(list_files[16:18],function(x){
  print(x)
  df = read.csv(x,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
  df$seqnames = factor(df$seqnames,
    levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))

  distribution  = data.frame(table(df$seqnames))
  distribution$peak_length = unlist(lapply(unique(df$seqnames), function(chr) {
    count <- sum(df[df$seqnames == chr,"width"])
    return(count)
  }))
  distribution$chr_length = length_chromosome$Total_length
  distribution$ratio = distribution$Freq/distribution$chr_length * 100000
  # distribution$ratio = distribution$peak_length/distribution$chr_length * 100
  # print(distribution)
  # return((mean(distribution$ratio[20:21])/mean(distribution$ratio[1:19]))-1)

  p = ggplot(distribution, aes(x = Var1,y = ratio, fill = Var1)) +
   geom_bar(stat = "identity") +
   scale_x_discrete(labels = c(1:19,"X","Y")) +
   ylim(0,15) +
   scale_fill_manual(values = c(rep("#424949",19),rep("#e74c3c",2))) +
   ggtitle(paste0("Distribution of peak (n =",nrow(df),")")) +
   xlab("chromosome") + ylab("nb peaks/chr length") +
   theme_bw() + theme(strip.background  = element_blank(),
     text = element_text(size=30, angle = 0),
     panel.grid.major = element_line(colour = "grey80"),
     panel.border = element_blank(),
     axis.ticks = element_blank(),
     panel.grid.minor.x=element_blank(),
     panel.grid.major.x=element_blank(),
     legend.position = "none")
 ggsave(filename = paste0(x,"_distribution_peak_chromosome_barplot.png"),
   plot = p,width = 10, height = 8, device = 'png', dpi = 150)
})))

list_data$mark = rep(sort(mark_list), each = 3)
list_data$cell = c("RS","SC","GS","RS","SC","GS","GS","RS","SC","RS","SC","GS","GS","RS","SC","RS",
  "SC","GS")
p = ggplot(list_data[4:15,], aes(x = mark, y = ratio, shape = cell)) +
  geom_point(size = 5) +
  xlab("") + ylab("XY/autosome") + labs(shape = "") +
  scale_shape_manual(values = c(0,1,4)) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 20, angle = 90),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank())
ggsave(filename = "XY_ratio_length.png", plot = p,width = 6, height = 7, device = 'png',dpi = 600)

####################################################################################################
##                              Ratio XY in histone mark in celltype
####################################################################################################

library(ggplot2)
library(ggpubr)
chr_order = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
length_chromosome = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", h = T, sep = ';',
  stringsAsFactors = FALSE)
mark_list = c("H3K4me1","H3K4me3","H3K9me3","H3K27ac","H3K27me3","H3K79me2")
list_SCRS = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
  full.names = T)[c(2:3,6:7,13:14,16:17,20:23)] ## SC - RS
list_GSSC = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
  full.names = T)[c(3:4,7:9,14,17:18,19,21,23:24)] ## GS - SC
list_data = do.call("rbind",lapply(mark_list,function(mark){
  GSSC = list_GSSC[grepl(mark, list_GSSC)]
  if(grepl("Kit",GSSC, ignore.case = T)) {
    GSSC = GSSC[c(2,1)]
  }
  table_GSSC = do.call("cbind",lapply(GSSC, function(file) {
    df = read.csv(file,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
    distribution  = data.frame(table(df$seqnames))
    return(distribution)
  }))
  table_GSSC$ratio = table_GSSC[,2]/table_GSSC[,4]
  table_GSSC$mark = mark
  table_GSSC$grp = "SC/GS"
  SCRS = list_SCRS[grepl(mark, list_SCRS)]
  table_SCRS = do.call("cbind",lapply(SCRS, function(file) {
      df = read.csv(file,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
      distribution  = data.frame(table(df$seqnames))
      return(distribution)
  }))
  table_SCRS$ratio = table_SCRS[,2]/table_SCRS[,4]
  table_SCRS$mark = mark
  table_SCRS$grp = "RS/SC"
  tt = rbind(table_GSSC,table_SCRS)
  p = ggplot(tt[,c(1,5:7)], aes(x = factor(Var1, level = chr_order), y = ratio, shape = grp)) +
    geom_point(size = 5) +
    # facet_grid(grp~mark, scales = "free") +
    xlab("") + ylab("ratio") + labs(shape = "") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=20, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 10, angle = 0),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())
  ggsave(filename = paste0(mark,"_ratio.png"), plot = p,width = 10, height = 6, device = 'png',
    dpi = 600)
  return(tt[c(20:21,41:42),])
}))

p = ggplot(list_data[,c(1,5:7)], aes(x = factor(Var1, labels = c("X","Y")), y = ratio, shape = grp)) +
  geom_point(size = 5) +
  facet_grid(grp~mark, scales = "free") +
  xlab("") + ylab("ratio") + labs(shape = "") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=20, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 20, angle = 0),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank())
ggsave(filename = "XY_ratio.png", plot = p,width = 10, height = 6, device = 'png',dpi = 600)

####################################################################################################
##                                         MSCI ANALYSIS
####################################################################################################

RS_expression = list_data[["RS"]]
RS_expression = RS_expression[RS_expression$gl == "Down-regulated (FDR < 5%)",]
RS_expression = RS_expression[RS_expression$chr == "chrX",]
RS_expression$gene = unlist(strsplit(RS_expression$gene_ENS,".",
  fixed = T))[seq(1,nrow(RS_expression)*2, 2)]

gene_expression = gene_expression[gene_expression$Chromosome == "X",]
rownames(gene_expression) = gene_expression$`Gene.Id`
gene_expression = gene_expression[RS_expression$gene,]


fviz_nbclust(log(gene_expression[,9:15]+0.01), kmeans, method = "wss", k.max = 15,linecolor = "black")

x = kmeans(log(gene_expression[,9:15]+0.01), centers = 5)
cl = x$cluster

gene_expression[names(cl),"cluster"] = cl

write.table(gene_expression,"DR_gene_expression_cluster.tsv", sep = "\t",col.names = T,
  row.names = T, quote = F)

df = do.call("rbind",lapply(colnames(gene_expression[,9:15]), function(col) {
  tt = gene_expression[,c("Gene.Id",col)]
  colnames(tt) = c("gene","expression")
  tt[names(cl),"cl"] = cl
  tt$celltype = col
  return(tt)
}))

df2 = data.frame(expression = unlist(lapply(colnames(gene_expression[,9:15]), function(col) {
  tt = gene_expression[,c("Gene.Id",col)]
  colnames(tt) = c("gene","expression")
  tt[names(cl),"cl"] = cl
  df = unlist(lapply(1:5, function(x) {
    mean(tt[tt$cl == x,"expression"])
  }))
  return(df)
})))
df2$celltype = rep(colnames(gene_expression[,9:15]), each = 5)
df2$cl = rep(1:5, 7)


p = ggplot(df, aes(x = factor(celltype,levels = colnames(gene_expression[,9:15])),
    y = log(expression+0.01),
    group = gene, color = factor(cl))) +
  geom_point() +
  geom_line() +
  xlab("") + ylab("log(RPKM)") + labs(color = "cluster") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 20),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank())
ggsave(filename = "DR_gene_expression.png", plot = p,width = 10, height = 8, device = 'png',
  dpi = 150)

p = ggplot(df2, aes(x = factor(celltype,levels = colnames(gene_expression[,9:15])),
    y = log(expression+0.01),
    group = cl, color = factor(cl))) +
  geom_point() +
  geom_line() +
  xlab("") + ylab("log(RPKM)") + labs(color = "cluster") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 20),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank())
ggsave(filename = "DR_mean_gene_expression.png", plot = p,width = 10, height = 8, device = 'png',
  dpi = 150)



## Create a table which count gene expression level in each chromosome
count_gene_chr = function(cell,chr,initial_file_ctl) {
  tt = initial_file_ctl[,grepl(paste0(".",cell,"_spike."),colnames(initial_file_ctl),fixed = T)]
  tt2 = initial_file_ctl[,grepl(paste0(".",cell,"."),colnames(initial_file_ctl),fixed = T)]
  new_table = cbind(tt,tt2)

  new_table$gene_expression = rowSums(new_table > expression) >= sample
  new_table$mean_gene_expression = rowSums(new_table[,1:length(new_table)-1])/(length(new_table)-1)

  new_table$gene_ENS = rownames(new_table)

  new_table_merge = merge(new_table,bed, by.x = "gene_ENS", by.y = "gene_ENS")
  new_table_merge = new_table_merge[new_table_merge$gene_expression == T,]
  chr_X = new_table_merge[new_table_merge$chr == chr,]

  return(mean(chr_X$mean_gene_expression))
}

## Params use to filtration (same that DEG)
sample = 2
expression = 0

initial_file = read.table("genecounts_cpm.tsv", sep = "\t", h = T)
expressed_genes = rowSums(initial_file > expression) >= sample
filtered_file = initial_file[expressed_genes,]
filtered_file_ctl = filtered_file[,grepl("Ctl",colnames(filtered_file),ignore.case = T)]
filtered_file_ko = filtered_file[,grepl("KO",colnames(filtered_file))]

## Create a table gene expression in each chromosome
table_nb_gene_chrX = c()
for(chr in unique(bed$chr)[-22]) {
  print(chr)
  table_nb_gene_chrX = c(table_nb_gene_chrX,
    unlist(lapply(c("Kit_m","Kit_p","SC","SCII","RS"),
    count_gene_chr,chr,filtered_file_ctl)))
}

table_nb_gene_chrX = data.frame("V1" = table_nb_gene_chrX)
table_nb_gene_chrX$sample = c(rep(samples[c(1:2,4:5,3)],length(unique(bed$chr))-1))
table_nb_gene_chrX$sample_simplify = rep(1:5,length(unique(bed$chr))-1)
table_nb_gene_chrX$chr = rep(unique(bed$chr)[-22],each = 5)
table_nb_gene_chrX$grp = c(rep("autosome",95),rep(c("chrX","chrY"),each = 5))

## Mean expression gene in function chr
p = ggplot(table_nb_gene_chrX, aes(x = sample_simplify, y = V1, color = grp, group = chr)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_continuous(labels = unique(table_nb_gene_chrX$sample)) +
  scale_color_manual(values = c("grey","green","orange")) +
  xlab("") + ylab("mean(gene expression)") + labs(color = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    # panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x= element_text(size = 30, hjust = 1),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank())
ggsave(filename = "MSCI_gene_expression.png", plot = p,width = 10, height = 8,
  device = 'png', dpi = 150)

## Comparison gene expression Kitp and RS
## RS initialization
# RS = initial_file_ctl[,grepl("RS",colnames(initial_file_ctl))]
# colnames(RS) = c("R1","R2","R3")
# RS$gene_ENS = rownames(RS)
# RS$gene_expression = rowSums(RS[,1:3] > expression) >= sample
# RS$mean_gene_expression = rowSums(RS[,1:3])/3
# RS_merge = merge(RS,bed, by.x = "gene_ENS", by.y = "gene_ENS")
# RS_merge = RS_merge[RS_merge$gene_expression == T,]
# RS_chrX = RS_merge[RS_merge$chr == "chrX",]
#
# ## Kit+ initialization
# Kitp = initial_file_ctl[,grepl("Kit_p",colnames(initial_file_ctl))]
# colnames(Kitp) = c("R1","R2","R3")
# Kitp$gene_ENS = rownames(Kitp)
# Kitp$gene_expression = rowSums(Kitp[,1:3] > expression) >= sample
# Kitp$mean_gene_expression = rowSums(Kitp[,1:3])/3
# Kitp_merge = merge(Kitp,bed, by.x = "gene_ENS",by.y = "gene_ENS")
# Kitp_merge = Kitp_merge[Kitp_merge$gene_expression == T,]
# Kitp_chrX = Kitp_merge[Kitp_merge$chr == "chrX",]
#
# common_gene_X = merge(RS_chrX, Kitp_chrX, by = "gene_ENS")[,c(1,6,16)]
# colnames(common_gene_X) = c("gene_ENS","RS_expression","Kitp_expression")
# common_gene_X$ratio_Kit_RS = common_gene_X$Kitp_expression / common_gene_X$RS_expression
# common_gene_X$log_ratio = log2(common_gene_X$ratio_Kit_RS)
# common_gene_X = common_gene_X[order(common_gene_X$log_ratio),]
# common_gene_X$gene_nb = 1:nrow(common_gene_X)
# common_gene_X$gl = "Kit = RS"
# common_gene_X[which(common_gene_X$log_ratio > 0),"gl"] = "Kit > RS"
# common_gene_X[which(common_gene_X$log_ratio < 0),"gl"] = "Kit < RS"
#
# df_RS = data.frame(table(round(log(common_gene_X$RS_expression))),"sample" = "RS")
# df_Kitp = data.frame(table(round(log(common_gene_X$Kitp_expression))),"sample" = "Kitp")
# df = rbind(df_RS,df_Kitp)
# df$Var1 = as.numeric(unfactor(df$Var1))
#
# p = ggplot(df,aes(x = Var1, y = Freq, group = sample, color = sample)) +
#   geom_line() +
#   scale_color_manual(values = c("red","blue")) +
#   geom_vline(xintercept = log(22.5837959), color = "red", linetype = "dashed") + ## Kitp niveau d'expression moyen dans X
#   geom_vline(xintercept = log(6.6054732), color = "blue", linetype = "dashed") + ## RS niveau d'expression moyen dans X
#   geom_vline(xintercept = log(0.5), linetype = "dashed") + ## cpm = 0.5
#   xlab("log(mean(gene expression))") + ylab("number of gene") + ylim(0,max(df$Freq)) +
#   theme_bw() + theme(strip.background  = element_blank(),
#     text = element_text(size=35, angle = 0),
#     panel.grid.major = element_blank(),
#     axis.ticks = element_blank())
# ggsave(filename = "distribution_Kit-RS_gene_expression.png", plot = p,width = 10, height = 8,
#   device = 'png', dpi = 150)

## chrX DEG gene expression in each cell type
deg_type = c("Down-regulated","Up-regulated")
deg = "Up-regulated"

table_cluster = lapply(deg_type, function(deg) {
  RS_deg = list_data$RS[list_data$RS$gl == deg,]
  RS_deg_merge = merge(RS_deg,bed, by.x = "gene_ENS",by.y = "gene_ENS")
  RS_deg_chrX = RS_deg_merge[RS_deg_merge$chr == "chrX",]

  gene_deg_chrX = filtered_file_ctl[RS_deg_chrX$gene_ENS,]

  list_gene_deg_chrX = lapply(c("Kit_m","Kit_p","SC","SCII","RS"),function(cell, file) {
    table = file[,grepl(paste0(".",cell,"."),colnames(file), fixed = T)]
    table2 = file[,grepl(paste0(".",cell,"_spike."),colnames(file), fixed = T)]
    table = cbind(table,table2)
    table$gene_expression = rowSums(table[,1:length(table)])/(length(table))

    return(table$gene_expression)
  }, file = gene_deg_chrX)
  table_gene_deg_chrX_expression = data.frame(list_gene_deg_chrX)
  colnames(table_gene_deg_chrX_expression) = c("Kit_m","Kit_p","SC","SCII","RS")
  rownames(table_gene_deg_chrX_expression) = RS_deg_chrX$gene_name

  # table_gene_deg_chrX_expression$cluster = c(1,2,4,3,3,1,4,1,4,4,1,4,2,3,3,4,2,2,2,2,2,2,2,2,2,2)
  #
  # table = c()
  # for(cluster in 1:4) {
  #   tt = t(table_gene_deg_chrX_expression[table_gene_deg_chrX_expression$cluster == cluster,1:5])
  #   table = c(table,rowSums(log(tt + 0.01))/length(colnames(tt)))
  # }
  # table = data.frame(CPM = table)
  # table$cluster = rep(1:4,each = 5)
  # table$sample = rep(1:5,4)
  # table$condition = "KO"

  ## Clusterize gene in function profil
  fviz_nbclust(log(table_gene_deg_chrX_expression+0.01), kmeans, method = "wss", k.max = 10,
    linecolor = "black")

  k = ifelse(deg == "Up-regulated",4,4)

  cl = kmeans(log(table_gene_deg_chrX_expression+0.01),k)
  cl_df = data.frame(factor(cl$cluster))

  color = list(group = c("#388e3c","#a5d6a7","#7b1fa2", "#ce93d8"))
  names(color$group) = unique(factor(cl$cluster))

  # df = data.frame(group = table_gene_deg_chrX_expression$cluster)
  # rownames(df) = rownames(table_gene_deg_chrX_expression)
  #
  # pheatmap(log(table_gene_deg_chrX_expression[order(table_gene_deg_chrX_expression$cluster),1:5]+0.01), cluster_cols = F, cluster_rows = F, scale = "none",
  #   filename = paste0("gene_expression_chrX_RS_",deg,"_heatmap_in_KO.png"),
  #   annotation_row = df,
  #   annotation_colors = color)

  pheatmap(log(table_gene_deg_chrX_expression[order(cl_df),]+0.01), cluster_cols = F,
    cluster_rows = F, scale = "none",
    filename = paste0("gene_expression_chrX_RS_",deg,"_heatmap.png"),
    annotation_row = data.frame(group = factor(cl$cluster)),
    annotation_colors = color)

  table_gene_deg_chrX = data.frame("expression" = unlist(list_gene_deg_chrX))
  table_gene_deg_chrX$sample = rep(1:5,each = length(RS_deg_chrX$gene_name))
  table_gene_deg_chrX$gene = rep(RS_deg_chrX$gene_name,5)
  table_gene_deg_chrX$logCPM = log(table_gene_deg_chrX$expression+0.01)

  cl_df = cl$centers
  table_cl = data.frame("CPM" = unlist(lapply(colnames(cl_df), function(x){
    return(cl_df[,x])
  })))
  table_cl$cluster = rep(1:k, 5)
  table_cl$sample = rep(1:5, each = k)

  p = ggplot(table_cl, aes(x = sample, y = CPM, color = factor(cluster))) +
    geom_point(size = 2) +
    geom_line() +
    xlab("") + ylab("log(CPM)") + labs(color = "cluster") +
    scale_x_continuous(labels = unique(table_nb_gene_chrX$sample)) +
    scale_color_manual(values = c("#388e3c","#a5d6a7","#7b1fa2", "#ce93d8")) +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=35, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())

  ggsave(paste0("gene_expression_chrX_RS_",deg,"_profil.png"), plot = p, dpi = 300, height = 8,
    width = 8, device = "png")
  return(table_cl)
})

# up_table = up_table[[2]]
# up_table$condition = "CTL"
# table_cl = rbind(up_table, table)

names(table_cluster) = c("Down-regulated","Up-regulated")

CTL = filtered_file_ctl
CTL$gene_ENS = rownames(CTL)
CTL_merge = merge(bed,CTL, by.x = "gene_ENS",by.y = "gene_ENS")
CTL_chrX = CTL_merge[CTL_merge$chr == "chrX",]

list_CTL_chrX = lapply(c("Kit_m","Kit_p","SC","SCII","RS"),function(cell, file) {
  table = file[,grepl(paste0(".",cell,"."),colnames(file), fixed = T)]
  table2 = file[,grepl(paste0(".",cell,"_spike."),colnames(file), fixed = T)]
  table = cbind(table,table2)
  table$gene_expression = rowSums(table[,1:length(table)])/(length(table))
  return(table$gene_expression)
}, file = CTL_chrX)
table_CTL_chrX_expression = data.frame(list_CTL_chrX)
colnames(table_CTL_chrX_expression) = c("Kit_m","Kit_p","SC","SCII","RS")
rownames(table_CTL_chrX_expression) = CTL_chrX$gene_name

## Clusterize gene in function profil
fviz_nbclust(log(table_CTL_chrX_expression+0.01), kmeans, method = "wss", k.max = 10,
  linecolor = "black")

cl = kmeans(log(table_CTL_chrX_expression+0.01),4)
cl_df = data.frame(factor(cl$cluster))

color = list(group = c("#388e3c","#a5d6a7","#7b1fa2", "#ce93d8"))
names(color$group) = unique(factor(cl$cluster))

table_CTL_chrX = data.frame("expression" = unlist(list_CTL_chrX))
table_CTL_chrX$sample = rep(1:5,each = length(CTL_chrX$gene_ENS))
table_CTL_chrX$gene = rep(CTL_chrX$gene_ENS,5)
table_CTL_chrX$logCPM = log(table_CTL_chrX$expression+0.01)

cl_df = cl$centers
table_cl = data.frame("CPM" = unlist(lapply(colnames(cl_df), function(x){
  return(cl_df[,x])
})))
table_cl$cluster = rep(1:4, 5)
table_cl$sample = rep(1:5, each = 4)

p = ggplot(table_cl, aes(x = sample, y = CPM, color = factor(cluster))) +
    geom_point(size = 2) +
    geom_line() +
    xlab("") + ylab("log(CPM)") + labs(color = "cluster") +
    scale_x_continuous(labels = unique(table_nb_gene_chrX$sample)) +
    scale_color_manual(values = c("#388e3c","#a5d6a7","#7b1fa2", "#ce93d8")) +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=35, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())

ggsave("gene_expression_chrX_RS_profil.png", plot = p, dpi = 300, height = 8,width = 8,
  device = "png")

## Merge DEG cluster and CTL cluster
table_down = data.frame("sample" = 1:5)
for(i in 1:5){
  table_down[i,"CPM"] = mean(table_cluster[[1]][table_cluster[[1]]$sample == i,"CPM"])
}

table_up = data.frame("sample" = rep(1:5,2))
x = c()
for(i in 1:5) {
  x = c(x,mean(table_cluster[[2]][table_cluster[[2]]$sample == i,"CPM"][1:3]))
}
table_up$CPM = c(x,table_cluster[[2]][table_cluster[[2]]$cluster == 4,"CPM"])

table_ctl = data.frame("sample" = rep(1:5,3))
table_ctl$CPM = c(table_cl[table_cl$cluster == 2,"CPM"],
  table_cl[table_cl$cluster == 3,"CPM"],
  table_cl[table_cl$cluster == 4,"CPM"])

table_cluster = do.call(rbind,list(table_down, table_up, table_ctl))
table_cluster$level = c(rep("down-regulated",5),rep("up-regulated",10),rep("CTL",15))
table_cluster$grp = c(rep("re-expressed",5),
  rep(c("new expressed","re-expressed"),each =5),
  rep(c("re-expressed","new expressed","no expressed"),each =5))
table_cluster$grp_cluster = paste0(table_cluster$grp,table_cluster$cluster)

p = ggplot(table_cluster,
  aes(x = sample, y = CPM, color = level,linetype = grp)) +
  geom_point() +
  geom_line() +
  labs(color = "") + ylab("log(CPM)") + xlab("") +
  scale_color_manual(values = c("black","blue","red")) +
  scale_linetype_manual(values = c("solid","solid","solid")) +
  scale_x_continuous(labels = unique(table_nb_gene_chrX$sample)) +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=35, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank()) +
  guides(linetype = F)
ggsave("gene_expression_chrX_RS_profil.png", plot = p, dpi = 300, height = 8,width = 8,
  device = "png")
