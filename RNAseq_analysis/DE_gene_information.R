####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

dir = "~/Documents/Mouse/RNAseq/Spike_analysis_Kit_SC_SCII_RS"
files = list.files(dir,"volcano_table_FDR.tsv", full.names = T)

sample_dir = c("Kit-","Kit+","RS","SC","SCII")
names(files) = sample_dir

DEG_files = lapply(files, function(file) {
  read.table(file, sep = "\t", h = T)
})

gtf = read.table("~/Documents/Annotations/gencode.vM19.annotation.gtf", sep = "\t", h = F)
colnames(gtf) = c("chr","source","feature","start","end","score","strand","frame","attribute")

## Keep only gene
gtf = gtf[gtf$feature == "gene",]
gtf$size = abs(gtf$start - gtf$end)

## Return only information about gene
df_attribute = t(data.frame(lapply(strsplit(gtf$attribute,";"),function(gene){
  x = gene[1:3]
  list_info = unlist(strsplit(x," "))[c(2,5,8)]
  return(list_info)
})))
rownames(df_attribute) = 1:nrow(df_attribute)
colnames(df_attribute) = c("gene_ENS","gene_type","gene_name")
gtf = cbind(gtf,df_attribute)

## Rename GTF gene type
gtf$gene_type_red = gtf$gene_type
gtf[grepl("IG_",gtf$gene_type),"gene_type_red"] = "IG - TR"
gtf[grepl("TR_",gtf$gene_type),"gene_type_red"] = "IG - TR"
gtf[grepl("protein_coding",gtf$gene_type),"gene_type_red"] = "Protein coding"
gtf[grepl("pseudogene",gtf$gene_type),"gene_type_red"] = "Pseudogene"
gtf[grepl("lncRNA",gtf$gene_type),"gene_type_red"] = "lncRNA"
gtf[grepl("processed_transcript",gtf$gene_type),"gene_type_red"] = "lncRNA"
gtf[grepl("sense",gtf$gene_type),"gene_type_red"] = "lncRNA"
gtf[grepl("Mt",gtf$gene_type),"gene_type_red"] = "mtRNA"
gtf[grepl("scaRNA",gtf$gene_type),"gene_type_red"] = "snoRNA"

####################################################################################################
##                                      CROSS GTF / DEG
####################################################################################################

list_DEG_gtf = lapply(names(DEG_files), function(f) {
  file = DEG_files[[f]]
  DEG_gtf = merge(file,gtf, by = "gene_ENS")
  size_quantile = quantile(DEG_gtf$size,probs = seq(0, 1, 0.1))
  for(size in size_quantile[-11]){
    DEG_gtf[DEG_gtf$size >= size,"grp"] = size
  }
  write.table(DEG_gtf, paste0(f,"_gene_information.tsv"), quote = F, row.names = F, col.names = T,
    sep = "\t")
  return(DEG_gtf)
})
names(list_DEG_gtf) = names(DEG_files)

####################################################################################################
##                                      DEG in fct SIZE
####################################################################################################

lapply(names(list_DEG_gtf), function(name) {
  file = list_DEG_gtf[[name]]
  df = data.frame(table(file$grp,file$gl))
  gl_distribution = data.frame(table(file$gl))
  df = merge(df,gl_distribution, by.x = "Var2", by.y = "Var1")
  colnames(df) = c("gl","grp","count","gl_dist")
  df$pct = df$count / df$gl_dist
  df$p_value = unlist(lapply(1:nrow(df), function(x) {
    chisq.test(c(df[x,"count"], df[x,"gl_dist"] - df[x,"count"]), p = c(0.1, 0.9))$p.value
  }))
  p = ggplot(df[df$gl != "Not regulated",], aes(x = grp, y = pct, fill = gl, group = gl)) +
    # geom_point() +
    geom_bar(position = "dodge",stat = "identity") +
    facet_grid(.~factor(gl), switch = "x") +
    geom_line(y = 0.10) +
    xlab("gene size (bp)") + ylab("") + labs(fill = "") +
    scale_fill_manual(values = c("blue","red")) +
    theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=30, angle = 0, hjust = 1),
      axis.text.x = element_text(size=15, angle = 90, hjust = 1),
      legend.position = "none")
  ggsave(paste0(name,"_DEG_gene_size.png"), device = 'png', plot = p, width = 10, height = 8,
  dpi = 150)
  write.table(df, paste0(name,"_DEG_gene_size.tsv"), quote = F, row.names = F, col.names = T,
    sep = "\t")
})

####################################################################################################
##                                      DEG in fct GENE TYPE
####################################################################################################

lapply(names(list_DEG_gtf), function(name) {
  file = list_DEG_gtf[[name]]
  df = data.frame(table(file$gene_type_red,file$gl))
  gl_distribution = data.frame(table(file$gl))
  df = merge(df,gl_distribution, by.x = "Var2", by.y = "Var1")
  colnames(df) = c("gl","gene_type","count","gl_dist")
  df$pct = df$count / df$gl_dist
  df$p_value = unlist(lapply(1:nrow(df), function(x) {
    not = df[df$gl == "Not regulated",]
    type = df[x,"gene_type"]
    theorie = c(not[not$gene_type == type,"pct"], 1- not[not$gene_type == type,"pct"])
    obs = c(df[x,"count"], df[x,"gl_dist"] - df[x,"count"])
    chisq.test(obs, p = theorie)$p.value
  }))
  p = ggplot(df, aes(x = gl, y = pct, fill = gene_type)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab("percentage") + labs(color = "Gene type") +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20),
      panel.grid.major = element_line(colour = "grey90"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=20, angle = 90, hjust = 0.5),
      axis.text.x = element_text(size=20))

  ggsave(paste0(name,"_DEG_gene_type.png"), device = 'png', plot = p, width = 10, height = 8,
  dpi = 150)
  write.table(df, paste0(name,"_DEG_gene_type.tsv"), quote = F, row.names = F, col.names = T,
    sep = "\t")
})

####################################################################################################
##                                  GENE SIZE in fct GENE TYPE
####################################################################################################

lapply(names(list_DEG_gtf), function(name) {
  file = list_DEG_gtf[[name]]
  df = data.frame(table(file$gene_type_red,file$grp))
  grp_distribution = data.frame(table(file$grp))
  df = merge(df,grp_distribution, by.x = "Var2", by.y = "Var1")
  colnames(df) = c("gene_size","gene_type","count","grp_dist")
  df$pct = df$count / df$grp_dist
  p = ggplot(df, aes(x = gene_size, y = pct, color = gene_type, group = gene_type)) +
    geom_point() +
    geom_line() +
    xlab("gene size (bp)") + ylab("percentage") + labs(color = "Gene type") +
    scale_color_brewer(palette = "Paired") +
    theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20),
      panel.grid.major = element_line(colour = "grey90"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=30, angle = 90, hjust = 1),
      axis.text.x = element_text(size=20))

  ggsave(paste0(name,"_gene_size_type.png"), device = 'png', plot = p, width = 10, height = 8,
  dpi = 150)
  p2 = ggplot(df, aes(x = gene_size, y = pct, fill = gene_type)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab("percentage") + labs(color = "Gene type") +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20),
      panel.grid.major = element_line(colour = "grey90"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=20, angle = 90, hjust = 0.5),
      axis.text.x = element_text(size=20))

  ggsave(paste0(name,"_gene_size_type_barplot_spike.png"), device = 'png', plot = p2, width = 10,
  height = 8,dpi = 150)

})

####################################################################################################
##                               DEG GENE SIZE in fct GENE TYPE
####################################################################################################

lapply(names(list_DEG_gtf), function(name) {
  file = list_DEG_gtf[[name]]
  df = data.frame(table(file$gene_type_red,file$grp, file$gl))
  grp_distribution = data.frame(table(file$grp))
  gl_distribution = data.frame(table(file$gl))
  colnames(df) = c("gene_type","gene_size","gl","freq")
  df = merge(df,grp_distribution, by.x = "gene_size", by.y = "Var1")
  df = merge(df,gl_distribution, by.x = "gl", by.y = "Var1")
  p = ggplot(df, aes(x = freq, y = gene_size, fill = gene_type)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(.~factor(gl), switch = "x") +
    ylab("gene size (bp)") + xlab("") + labs(fill = "Gene type") +
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20),
      panel.grid.major = element_line(colour = "grey90"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(size=20, hjust = 1),
      axis.text.x = element_text(size=15))

  ggsave(paste0(name,"_DEG_gene_information_spike.png"), device = 'png', plot = p, width = 12, height = 8,
  dpi = 150)
})
