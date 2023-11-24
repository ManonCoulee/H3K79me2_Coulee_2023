####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2','stringr','RColorBrewer','ggpubr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

dir = "~/Documents/Mouse/RNAseq/Spike_analysis_Kit_SC_SCII_RS/"

samples = c("Kit-","Kit+","RS","SCI","SCII")
rout = file.path(".")

## RNAseq value from Gan et al paper
gene_expression = read.table("Gan_2013_RPKM.csv",stringsAsFactors = F, sep = "\t", h = T)

chr_order = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
length_chromosome = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", h = T, sep = ';',
  stringsAsFactors = FALSE)
mark_list = c("H3K4me1","H3K4me3","H3K9me3","H3K27ac","H3K27me3","H3K79me2")

list_files = list.files(c(dir),
  pattern = "_volcano_table_FDR.tsv", full.names = T)
names(list_files) = samples
list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})
names(list_data) = samples

####################################################################################################
##                              Gene expression distribution
####################################################################################################

## Gene expression distribution according to FC value
create_genes_deregulate_FC_plot = function(name,out, listData) {
  data = listData[[name]]
  data = data[data$PValue <= 0.05,c("chr","logFC","PValue","FDR","gl")]
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
tableFC$cellType = factor(tableFC$cellType, levels = c("Kit-","Kit+","SCI","SCII","RS"))

## Boxplot with all cell type
p = ggplot(tableFC, aes(x = cellType, y = logFC, fill = chromosome)) +
  geom_boxplot(aes(facet.by = cellType)) +
  labs(fill = "") + xlab("") + ylab("log2FC KOvsCTL") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = -1.5, linetype = "dashed") +
  scale_fill_manual(values = c("DodgerBlue3","orange")) +
  theme_bw() + theme(
    text = element_text(size=30, angle = 0),
    legend.position = "bottom",
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank()) +
    stat_compare_means(label = "p.format", method = "wilcox.test")
ggsave(filename = "GSC_SC_SCII_RS_FC_spike.png", plot = p, width = 6, height = 5, device = 'png',
  dpi = 600)

####################################################################################################
##                              Ratio XY/autosome in histone mark
####################################################################################################

list_files = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
  full.names = T)[c(2:4,6:8,10,13,14,16:24)]
list_data = data.frame(ratio = unlist(lapply(list_files,function(x){
  print(x)
  df = read.csv(x,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
  df$seqnames = factor(df$seqnames,
    levels = chr_order)

  distribution  = data.frame(table(df$seqnames))
  distribution$peak_length = unlist(lapply(unique(df$seqnames), function(chr) {
    count <- sum(df[df$seqnames == chr,"width"])
    return(count)
  }))
  distribution$chr_length = length_chromosome$Total_length
  distribution$ratio = distribution$Freq/distribution$chr_length * 100000

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

## Plot contains all data
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

list_SCRS = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
  full.names = T)[c(2:3,6:7,13:14,16:17,20:23)] ## SC - RS
# list_GSSC = list.files(mark_list, pattern = "peak_annotate_reduce.tsv",
#   full.names = T)[c(3:4,7:9,14,17:18,19,21,23:24)] ## GS - SC
list_data = do.call("rbind",lapply(mark_list,function(mark){
  # GSSC = list_GSSC[grepl(mark, list_GSSC)]
  # if(grepl("Kit",GSSC, ignore.case = T)) {
  #   GSSC = GSSC[c(2,1)]
  # }
  # table_GSSC = do.call("cbind",lapply(GSSC, function(file) {
  #   df = read.csv(file,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
  #   distribution  = data.frame(table(df$seqnames))
  #   return(distribution)
  # }))
  # table_GSSC$ratio = table_GSSC[,2]/table_GSSC[,4]
  # table_GSSC$mark = mark
  # table_GSSC$grp = "SC/GS"
  SCRS = list_SCRS[grepl(mark, list_SCRS)]
  table_SCRS = do.call("cbind",lapply(SCRS, function(file) {
      df = read.csv(file,sep = "\t", h = T, stringsAsFactors = F)[,c(1,4)]
      distribution  = data.frame(table(df$seqnames))
      return(distribution)
  }))
  table_SCRS$ratio = table_SCRS[,2]/table_SCRS[,4]
  table_SCRS$mark = mark
  table_SCRS$grp = "RS/SC"
  tt = table_SCRS
  # tt = rbind(table_GSSC,table_SCRS)
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
  ggsave(filename = paste0(mark,"_ratio_SCRS.png"), plot = p,width = 10, height = 6, device = 'png',
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
