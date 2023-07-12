####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2','RColorBrewer','stringr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                            INITIALIZATION
####################################################################################################

dir = "../GSEA_analysis_SC_SCII_RS_spike"
rout = file.path(".")
samples = c("CTL_GObp")

####################################################################################################
##                                                 DATA
####################################################################################################

list_files = list.files(dir,pattern = "interest", full.names = T)
names(list_files) = samples

list_data = lapply(list_files,function(x){
  read.table(x,sep = "\t", h = T, stringsAsFactors = F)
})

####################################################################################################
##                                                 PLOT
####################################################################################################

 color_table = c("black" = "grey","blue" = "#3498db","orange" = "#ffb74d","red" = "#e74c3c",
   "yellow" = "#58d68d", "purple" = "#9b59b6", "green" = "#e74c3c", "pink" =  "#f8bbd0",
   "grey" = "#000000")

lapply(samples,function(sample) {
  print(sample)
  data = list_data[[sample]]
  list_subdata = list(data[,c(1:4)],data[,c(5:8)],data[,c(9:12)],data[,c(13:16)], data[,c(17:20)])

  table_data = do.call(rbind,lapply(list_subdata,function(x) {
    name = toupper(unlist(strsplit(colnames(x)[1],"_GO"))[1])
    colnames(x) = c("pathway","ES","p-value", "grp")
    x = x[!is.na(x$ES),]
    x$pathway = str_replace_all(x$pathway,'_',' ')
    x$pathway = tolower(str_replace(x$pathway,'GO..',''))
    x$sample = name

    p = ggplot(x, aes(y = factor(pathway, levels = pathway[length(pathway):1]), x = abs(ES), fill = grp)) +
      geom_col(width = 0.5) +
      xlab("") + ylab("Enrichment Score (p.value < 0.05)") +
      scale_fill_manual(values = color_table) +
      theme_bw() + theme(
        axis.text.y = element_text(size = 40),
        text = element_text(size=50, angle = 0),
        panel.border = element_rect(colour = 'grey50', fill = NA),
        panel.grid.major.x = element_line(colour = "grey80"),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none")
    ggsave(plot = p, filename = paste0(sample,"_",name,"_summary.png"), device = "png", dpi = 300,
      height = 20,width = 20)

    return(x)
  }))

  ## Keep only interest pathway
  table_data = table_data[!table_data$grp == "black",]
  p = ggplot(table_data, aes(x = pathway, y = factor(sample, levels = unique(sample)),
    color = grp,size = abs(ES))) +
    geom_point() +
    coord_flip() +
    xlab("") + ylab("Enrichment Score (p.value < 0.05)") + labs(size = "ES") +
    scale_color_manual(values = color_table) +
    scale_size_continuous(range = c(5,10)) +
    theme_bw() + theme(
      axis.text.y = element_text(size = 40),
      text = element_text(size=40, angle = 0),
      panel.border = element_rect(colour = 'grey50', fill = NA),
      panel.grid.major.x = element_line(colour = "grey80"),
      axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom") +
    guides(color = "none")

  ggsave(plot = p, filename = paste0(sample,"_summary.png"), device = "png", dpi = 300,
      height = 25,width = 30)
})

code_name = c("pathways related to chromatin organization, chromosome",
  "pathways related to mitochondria activity",
  "pathways related to transcription and RNA/DNA regulation",
  "pathways related to spermatozoa formation",
  "pathways related to spermatogenesis",
  "pathways related to ribosome",
  "pathways related to apoptosis",
  "other")
color_code = c("#58d68d","#ffb74d","#9b59b6","#e74c3c","#3498db","#f8bbd0", "#000000","grey")

df = data.frame(color_code, code_name,expression = 1)

p = ggplot(df, aes(x = expression, y = code_name,
  fill = factor(code_name, levels = code_name),
  color = factor(code_name, levels = code_name))) +
  geom_col(width = 0.2,position = position_dodge(0.1)) +
  xlab("") + ylab("") +
  scale_fill_manual(values = color_code) +
  scale_color_manual(values = color_code) +
  theme_bw() + theme(
    text = element_text(size=60, angle = 0),
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank())

ggsave("legend.png", plot = p, device = "png",dpi = 300, height = 8, width = 20)
