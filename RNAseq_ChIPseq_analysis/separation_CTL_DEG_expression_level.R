####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('gridExtra','ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'reshape',
  'TCseq', 'statmod', 'GenomicFeatures', 'VennDiagram', 'pheatmap')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

dir = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_SC_SCII_RS_spike"
dir = "~/Documents/Mouse/RNAseq/DE_analysis_3_factors_Kit"
dir = "."
CT = c('RS','SC')
CT = "Kit-"
bed = read.table("~/Documents/Annotations/gencode.vM19.annotation.bed",h=F, sep = '\t')
colnames(bed) = c("chr","start","end", "strand","score","gene_ENS")
chip = "~/Documents/Mouse/ChIPseq/H3K79me2/"

files_list = list.files(dir, pattern = "_CTL_gene_expression.tsv", full.name = T)[1:4]
names(files_list) = CT

####################################################################################################
##                                       CTL expression separation
####################################################################################################

## Separe bed in function CTL expression levels
separe_bed_genes_grp = function(name_file,bed,type) {
  file = read.table(name_file, sep = "\t", h = T)
  q = as.data.frame(quantile(file$gene_expression,probs = seq(0,1,0.1)))
  q$grp = seq(0:10)
  colnames(q) = c("quantile","grp")
  q[11,"quantile"] = q[11,"quantile"]+1

  file$grp = apply(file,1,function(x,q) {
    mean = x[["gene_expression"]]
    # print(mean)
    i = 1
    while(as.numeric(mean) >= q[i,"quantile"]) {
      res = paste0(q[i,"grp"],'-',q[i+1,"grp"])
      i = i+1
    }
    return(res)
  },q = q)

  bed_merge = merge(bed,file,by.x = "gene_ENS", by.y = "gene_ENS")
  write.table(bed_merge[,c(2:4,1,5:6)],paste0(type,"/",type,"_CTL.bed"),row.names = F, quote = F,
    col.names = F, sep = "\t")
  level = sort(unique(bed_merge$grp))
  for(i in level) {
    grp = unlist(strsplit(i,"-"))[[1]]
    write.table(bed_merge[bed_merge$grp == i,c(2:4,1,6,5)],paste0(type,"/",type,"_CTL_grp",grp,".bed"),
      row.names = F, quote = F, col.names = F, sep = "\t")
  }
  return(file)
}

SC_genecounts = separe_bed_genes_grp(files_list[["SC"]],bed,"SC")
RS_genecounts = separe_bed_genes_grp(files_list[["RS"]],bed,"RS")
Kit_m_genecounts = separe_bed_genes_grp(files_list[["Kit-"]],bed,"Kit-")
Kit_p_genecounts = separe_bed_genes_grp(files_list[["Kit_p"]],bed,"Kit_p")

####################################################################################################
##                                          DEG separation
####################################################################################################

## Separe bed in function DEG genes in RS and SC
separe_bed_gene_DEG = function(type,bed) {
  DE_table = read.table(paste0(dir,"/",type,"_spike_volcano_table_FDR.tsv"),
    h = T, sep = '\t', stringsAsFactors = F)
  # bed_merge = merge(bed,DE_table,by.x = "gene_ENS", by.y = "gene_ENS")[c(2:4,1,6,5,8,14)]
  level = unique(DE_table$gl)
  for(i in level) {
    tt = DE_table[DE_table$gl == i,]
    tt = tt[tt$strand == "+",]
    # write.table(tt[,1:6],paste0(type,"/",type,"_KO_",i,".bed"), sep = "\t", quote = F,
    write.table(tt[,c(2:6,1)],paste0(type,"_KO_",i,"_FDR.bed"), sep = "\t", quote = F,
      row.names = F, col.names = F)
  }
}

separe_bed_gene_DEG("RS",bed)
separe_bed_gene_DEG("SC",bed)
separe_bed_gene_DEG("Kit_m",bed)
separe_bed_gene_DEG("Kit_p",bed)
