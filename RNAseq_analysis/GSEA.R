####################################################################################################
##                                               LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('stringr', 'tidyverse')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                             INITIALIZATION
####################################################################################################

dirs = list.dirs(".",full.names = TRUE,recursive = FALSE)

wd = '.'

GObp = dirs[grepl(paste0("GObp.Gsea"),dirs)]

pattern_CTL = "gsea_report_for_CTL"

####################################################################################################
##                                          CREATE TABLE ANALYSE
####################################################################################################

## Function to merge in one table GSEA informations (name, ES score and p-value) for each samples
createTableGSEA = function(x, db){
  samplefiles = unlist(lapply(db,function(dir) {
    x1 = list.files(dir, pattern = "tsv", full.names = T)
    x2 = list.files(dir, pattern = x, full.names = T)
    x = intersect(x1,x2)
  }))#[c(2,3,1)]

  n = lapply(samplefiles,function(x) {unlist(strsplit(x,"/"))[length(unlist(strsplit(x,"/")))-1]})
  names(samplefiles) = lapply(n, function(x) {unlist(strsplit(x,".",fixed = T))[1]})
  # names(samplefiles) = unlist(lapply(names(samplefiles),function(x) {str_replace(x,"_GO","")}))

  list = lapply(samplefiles,read.table,h = T, sep = '\t', stringsAsFactors = F)
  list = lapply(list,function(x) {x[x$NOM.p.val < 0.05,]})
  list = lapply(list,function(x) {x[order(abs(x$ES), decreasing = T),]})
  list = lapply(list, function(x) {x[,c("NAME","ES","NOM.p.val")]})

  max = max(unlist(lapply(list, function(x) {nrow(x)})))
  addLineTable = function(x,max) {
    if (nrow(x) < max){
      x[(nrow(x)+1):max,] = NA
    }
    return(x)
  }

  list = lapply(list,addLineTable,max)
  table = bind_cols(list)

  colnames(table) = paste0(rep(names(samplefiles),each=3),rep(c("_NAME","_ES","_NOM.p.val"),
    length(names(samplefiles))),sep = "")
  return(table)
}

## Write in data frame
write.table(createTableGSEA(pattern_CTL, db = GObp),"summary_GSEA_CTL_GObp.tsv",sep = '\t',  col.names = T, row.names = F, quote = F)

