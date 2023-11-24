####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ChIPseeker','stringr', 'Rgraphviz','ggplot2','topGO','GenomicFeatures', 'org.Mm.eg.db',
  'clusterProfiler','VennDiagram')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

args = commandArgs(trailingOnly=TRUE)
wd = "."
current_dir = rev(unlist(strsplit(getwd(),"/")))[1]
txmm = makeTxDbFromGFF(file = file.path('~/Documents/Annotations/gencode.vM19.annotation.gtf'),
  format = 'gtf')
length_chromosome = read.table("~/Documents/Annotations/mouse_length_chromosome.tsv", h = T, sep = ';',
  stringsAsFactors = FALSE)

replicat = c("Kit","RS","SC")

common_files = list.files(replicat, pattern = "_peak_common.bed", full.names = T)
names(common_files) = replicat

names_common = names(common_files)
merge_files = list.files("ES", pattern = "_peak_merge.bed", full.names = T)
names(merge_files) = "ES"
names_merge = names(merge_files)

chr_order = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

load("~/Documents/Annotations/MouseTableBed")

####################################################################################################
##                                         INDIVIDUAL ANALYSIS
####################################################################################################

coverage = function(data) {
  data$len = data$end - data$start
  return(data)
}

chipseq_analysis = function(files) {
  lapply(names(files),function(name,list_file,txmm){

    ## Initialization
    print(paste0("##########################  ",name))
    tt = read.table(list_file[[name]], h = F,sep = "\t", stringsAsFactors = FALSE)
    if(name == "ES") {
      colnames(tt) =  c("chr","start","end","peak_name","score","strand","signalValue","p-value",
      "FDR","position")
    } else {
      colnames(tt) =  c("chr","start","end","peak_name","score","strand","signalValue","p-value",
        "FDR")
    }
    tt = unique(tt)

    ## Chromosome distribution
    print("Chromosome distribution")
    tt = tt[tt$chr != "chrMT",]
    tt$chr = factor(tt$chr,levels = chr_order)
    distribution  = as.data.frame(table(tt$chr))
    colnames(distribution) = c("chr","length")
    distribution$chr_length = length_chromosome$Total_length
    distribution$ratio = distribution$length/distribution$chr_length*100000 ## 10-5
    p = ggplot(distribution, aes(x = chr,y = ratio, fill = chr)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(labels = c(1:19,"X","Y")) +
      scale_fill_manual(values = c(rep("#424949",19),rep("#e74c3c",2))) +
      ggtitle(paste0("Distribution of peak (n =",nrow(tt),")")) +
      ylim(0,6) +
      xlab("chromosome") + ylab("nb peaks/chr length") +
      theme_bw() + theme(strip.background  = element_blank(),
        text = element_text(size=30, angle = 0),
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "none")
    ggsave(filename = paste0(wd,"/",name,"/",name,"_distribution_peak_chromosome_barplot.png"),
      plot = p, width = 10, height = 8, device = 'png', dpi = 150)

    ## Annotation feature
    print("Annotation")
    tt_peak = makeGRangesFromDataFrame(tt, keep.extra.columns = TRUE)
    tt_annotate = annotatePeak(tt_peak, tssRegion = c(-3000,3000), TxDb = txmm, annoDb = "org.Mm.eg.db")
    tt_anno = tt_annotate@anno

    tt_anno$gene_ENS = unlist(strsplit(tt_anno$geneId,".", fixed = T))[seq(1,length(tt_anno)*2,2)
    tt_anno$gene_name = lapply(tt_anno$gene_ENS, function(gene) {
      unique(MouseTableBed[MouseTableBed$ensembl_gene_id == gene,"mgi_symbol"])
    })

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

  },list_file = files,txmm = txmm)
}

common_files = chipseq_analysis(common_files)
merge_files = chipseq_analysis(merge_files)

common_merge_files = c(common_files, merge_files)
names(common_merge_files) = c(names_common,names_merge)
names(common_files) = names_common

####################################################################################################
##                                         COMMON ANALYSIS
####################################################################################################

p = ggplot(do.call("rbind",common_files), aes(x = Freq, fill = Var1,
  y = factor(name, levels = c("GS","SC","RS")))) +
  geom_bar(stat = "identity", position = position_fill()) +
  scale_fill_manual(values = c("#B2BABB","#BB8FCE","#E67E22","#2980B9")) +
  ylab("") + xlab("") +
  theme_minimal() +
      theme(text = element_text(size = 40),legend.position = "bottom",legend.title = element_blank(),
        axis.text.x = element_blank()) +
      guides(fill = guide_legend(nrow=2, byrow = T))
ggsave(filename = "GS_SC_RS_annotation_distribution_barplot_reduce.png",plot = p,
  width = 10, height = 5, device = 'png', dpi = 450)
