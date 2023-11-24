####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'edgeR', 'DESeq2', 'reshape', 'TCseq', 'statmod',
  'GenomicFeatures', 'VennDiagram', 'ggrepel')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

dir = "~/Documents/Mouse/RNAseq/Spike_analysis_Kit_SC_SCII_RS"

rout = file.path(".")
samples = c("Kit-","Kit+","RS","SC","SCII")

####################################################################################################
##                                                 DATA
####################################################################################################

list_files = list.files(dir,
  pattern = "spike_volcano_table_FDR.tsv", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})

repressor_file = c("negative_regulation","negative_transcription_of_regulation_RNA_polymerase_II",
  "repressor")
activator_file = c("positive_regulation_of_transcription_by_RNA_polII")
repressor_genes = data.frame(gene = unique(unlist(lapply(repressor_file, function(file) {
  read.table(file, h = F, sep = "\t")
}))))
activator_genes = data.frame(gene = unique(unlist(lapply(activator_file, function(file) {
  read.table(file, h = F, sep = "\t")
}))))

tt = theme(
    plot.title = element_text(size = rel(2), lineheight = 2, vjust = 1, face = 'bold'),
    plot.subtitle = element_text(size = rel(1.2), face = 'bold'),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5)),
    panel.border = element_rect(colour = 'grey50', fill = NA))


####################################################################################################
##                                          CORRELATION PLOT
####################################################################################################

lapply(names(list_data), function(name) {
  print(name)
  H3K79me2 = read.csv(paste0("~/Documents/Mouse/ChIPseq/H3K79me2/",name,"/",name,"_peak_annotate.tsv"),
    sep = "\t",h = T)
  data = list_data[[name]]
  rownames(data) = data$gene_ENS
  data_rep = merge(data,repressor_genes, by.x = "gene_name", by.y = "gene")
  rownames(data_rep) = data_rep$gene_ENS
  data_act = merge(data,activator_genes, by.x = "gene_name", by.y = "gene")
  rownames(data_act) = data_act$gene_ENS

  data$repressor = FALSE
  data[rownames(data_rep),"repressor"] = TRUE
  data[data$gl == "Not regulated", "repressor"] = FALSE

  data$activator = FALSE
  data[rownames(data_act),"activator"] = TRUE
  data[data$gl == "Not regulated", "activator"] = FALSE

  ## Presence of H3K79me2 (1) else 0
  data$H3K79me2 = 0
  data[unique(H3K79me2$geneId),"H3K79me2"] = 1
  data = data[!is.na(data$gl),]

  write.table(data, paste0(name,"_REP_ACT_genes_volcano_table_FDR.tsv"), quote = F, row.names = F,
    col.names = T, sep = "\t")
})
