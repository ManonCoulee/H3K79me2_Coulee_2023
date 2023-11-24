####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ChIPseeker','stringr', 'Rgraphviz','ggplot2','topGO','GenomicFeatures', 'org.Mm.eg.db',
  'EnsDb.Mmusculus.v79', 'clusterProfiler')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                                  DATA
####################################################################################################

args = commandArgs(trailingOnly=TRUE)
wd = args

args = "ES"
wd = "."
txmm = makeTxDbFromGFF(file = file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'),
  format = 'gtf')
length_chromosome = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", h = T, sep = ';',
  stringsAsFactors = FALSE)

samples_dir = list.dirs(wd, full.names = F, recursive = FALSE)

# print(args)
if(args == "ES") {
  sample_files = list.files(samples_dir, pattern = "narrowPeak", full.names = T)
} else {
  sample_files = list.files(samples_dir, pattern = "broadPeak", full.names = T)
}

if(args == "Kit") {
  names(sample_files) = c("A_1ug_H3K79me2","A_2ug_H3K79me2","B_1ug_H3K79me2","B_2ug_H3K79me2")
} else {
  names(sample_files) = c(paste0(args,"_1_H3K79me2"),paste0(args,"_2_H3K79me2"))
}

names(sample_files) = str_replace(samples_dir,"_macs2","")
names(sample_files) = c("B","Undiff")


chr_order = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

load("~/Documents/Annotations/MouseTableBed")

####################################################################################################
##                                           INDIVIDUAL ANALYSIS
####################################################################################################

peak_anno_list = lapply(names(sample_files),function(name,list_file,txmm) {

  ## Initialization
  print(paste0("##########################  ",name))
  tt = read.table(list_file[[name]], h = F, sep = '\t', stringsAsFactors = F)
  if(args == "ES") {
    colnames(tt) =  c("chr","start","end","peak_name","score","strand","signalValue","p-value",
      "FDR","position")
  } else {
    colnames(tt) =  c("chr","start","end","peak_name","score","strand","signalValue","p-value",
      "FDR")
    # colnames(tt) =  c("chr","start","end","peak_name","score","strand","signalValue","p-value",
    #   "FDR","position")
  }

  ## Remove in "chr" column chromosome exogene
  tt = tt[grepl("chr",tt$chr),]
  tt = tt[!grepl("chrM",tt$chr),]
  tt$chr = factor(tt$chr,levels = chr_order)

  ## Remove lower signifiant FDR
  # tt = tt[tt$FDR < 5,]
  # subset(tt, 'FDR' < 3)
  # tt = subset(tt, 'FDR' > 0.05)

  ## File for merge analysis
  write.table(tt,paste0(wd,"/",name,"_peak.bed"),sep = '\t',quote = FALSE,
   col.names = TRUE, row.names = FALSE)

  ## Distribution of peak in chromosome
  print("Chromosome distribution")
  tt$chr = factor(tt$chr,levels = chr_order)
  distribution  = as.data.frame(table(tt$chr))
  colnames(distribution) = c("chr","length")
  distribution$chr_length = length_chromosome$Total_length
  distribution$ratio = distribution$length/distribution$chr_length*100000 ## 10-5
  p = ggplot(distribution, aes(x = chr,y = ratio, fill = chr)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(labels = c(1:19,"X","Y")) +
    ylim(0,5) +
    scale_fill_manual(values = c(rep("#424949",19),rep("#e74c3c",2))) +
    ggtitle(paste0("Distribution of peak (n =",nrow(tt),")")) +
    xlab("chromosome") + ylab("nb peaks/chr length") +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=30, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      legend.position = "none")
  ggsave(filename = paste0(wd,"/",name,"_distribution_peak_chromosome_barplot.png"),
    plot = p,width = 10, height = 8, device = 'png', dpi = 150)

  print("Annotation")
  tt_peak = makeGRangesFromDataFrame(tt, keep.extra.columns = TRUE)
  tt_annotate = annotatePeak(tt_peak, tssRegion = c(-3000,3000), TxDb = txmm, annoDb = "org.Mm.eg.db")
  tt_anno = tt_annotate@anno

  tt_anno$gene_ENS = unlist(strsplit(tt_anno$geneId,".", fixed = T))[seq(1,length(tt_anno)*2,2)]

  tt_anno$gene_name = lapply(tt_anno$gene_ENS, function(gene) {
    unique(MouseTableBed[MouseTableBed$ensembl_gene_id == gene,"mgi_symbol"])
  })

  write.table(as.data.frame(tt_anno),paste0(wd,"/",name,"_peak_annotate.tsv"),sep = '\t',quote = FALSE,
   col.names = TRUE, row.names = FALSE)

  print(paste0("Number of annotate peaks: ",tt_annotate@peakNum))
  png(paste0(wd,"/",name,"_annotation_distribution_piechart.png"),width = 600, height = 180)
  plotAnnoPie(tt_annotate, title = "Feature distribution")
  dev.off()

  ## Reduction of annotation name
  tt_reduce = as.data.frame(tt_anno)
  tt_reduce$reduce = "Intragenic"
  tt_reduce[grepl("Downstream",tt_reduce$annotation),"reduce"] = "Downstream"
  tt_reduce[grepl("Promoter",tt_reduce$annotation),"reduce"] = "Promoter"
  tt_reduce[grepl("Distal",tt_reduce$annotation),"reduce"] = "Distal Intergenic"
  print(paste0("Number of gene : ", length(unique(tt_reduce$gene_name))))
  df_reduce = data.frame(table(tt_reduce$reduce))
  print(df_reduce)
  print(round(df_reduce$Freq/tt_annotate@peakNum*100,2))

  write.table(tt_reduce,paste0(wd,"/",name,"_peak_annotate_reduce.tsv"), sep = "\t", quote = F,
    col.names = T, row.names = F)

  p = ggplot(df_reduce,aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", position = position_fill()) +
    coord_polar("y", start=0) + ylab("") + xlab("") +
    scale_fill_manual(values = c("#B2BABB","#BB8FCE","#E67E22","#2980B9")) +
    theme_minimal() +
    theme(text = element_text(size = 40),legend.position = "bottom",legend.title = element_blank(),
      axis.text.x = element_blank()) +
    guides(fill = guide_legend(nrow=2, byrow = T))
  ggsave(filename = paste0(wd,"/",name,"_annotation_distribution_piechart_reduce.png"),plot = p,
    width = 10, height = 5, device = 'png', dpi = 450)

  df_reduce$name = name
  return(df_reduce)

},list_file = sample_files, txmm = txmm)

####################################################################################################
##                                           ANNOTATION ANALYSIS
####################################################################################################

names(peak_anno_list) = names(sample_files)
p = ggplot(do.call("rbind",peak_anno_list), aes(y = name, x = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_fill()) +
  scale_fill_manual(values = c("#B2BABB","#E67E22","#2980B9","#BB8FCE")) +
  ylab("") + xlab("") +
  theme_minimal() +
      theme(text = element_text(size = 40),legend.position = "bottom",legend.title = element_blank(),
        axis.text.x = element_blank()) +
      guides(fill = guide_legend(nrow=2, byrow = T))
ggsave(filename = paste0(args,"/",args,"_annotation_distribution_barplot_reduce.png"),plot = p,
  width = 15, height = 5, device = 'png', dpi = 450)
