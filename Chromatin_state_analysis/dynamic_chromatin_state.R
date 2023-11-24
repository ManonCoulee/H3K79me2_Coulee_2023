####################################################################################################
##                                        LIBRARIES
####################################################################################################

rm(list = ls())
options(warn = -1, width = 150)

rlibs = c('ggplot2', 'RColorBrewer', 'ggalluvial', 'stringr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

####################################################################################################
##                                        INITIALIZATION
####################################################################################################

colorValue = c(rep(c("#B71C1C","#E65100","#E65100","#43A047","#1B5E20","#99FF66","#99FF66","#F5B041",
  "#FFEB3B","#48C9B0","#B39DDB","#880E4F","#666633","#424949","#424949","#D0D3D4","#D0D3D4",
  "#D0D3D4"), 2),"black","NA")
names(colorValue) = c(paste0("U",1:18),paste0("U",1:18,"-H3K79me2"),"+H3K79me2","-H3K79me2")

levels = c(paste0("U",1:18))
labels = c("active","weak","weak","active","weak","active","active","active","weak","het","het",
    "bivalent","bivalent","polycomb","polycomb","quiescent","quiescent","quiescent")
stateName = c("TSSA","TSSFlnk","TSSFlnk","Tx","TxWk","EnhG","EnhG","EnhA","EnhW","ZNF","Het",
    "TssBiv","EnhBiv","ReprPC","ReprPC","Quies","Quies","Quies")

####################################################################################################
##                                  CHARACTERIZATION DYNAMIC STATE
####################################################################################################

## Download files
region = read.table("Kit_SC_RS_open_region_GSC_SC_RS_V2.tsv", sep = "\t", h = F)
colnames(region) = c("chr","start","end","GSC_peak","SCI_peak","RS_peak","motif_peak","GSC_state",
  "GSC_coverage","SCI_state","SCI_coverage","RS_state","RS_coverage")

## Keep the predominant state in region (max coverage)
region$pos_peak = paste0(region[,1],region[,2],region[,3])
region_table = do.call("rbind",lapply(unique(region$pos_peak), function(pos) {
  df = region[region$pos_peak == pos,]
  df_region = df[1,1:7]
  GSC = aggregate(df$GSC_coverage, list(df$GSC_state), sum)
  df_region$GSC_state = GSC[GSC$x == max(GSC$x),"Group.1"][1]
  SC = aggregate(df$SCI_coverage, list(df$SCI_state), sum)
  df_region$SCI_state = SC[SC$x == max(SC$x),"Group.1"][1]
  RS = aggregate(df$RS_coverage, list(df$RS_state), sum)
  df_region$RS_state = RS[RS$x == max(RS$x),"Group.1"][1]
  return(df_region)
}))

## Create column to have state name
region_table$RS_stateName = factor(region_table$RS_state, levels = levels, labels = labels)
region_table$GSC_stateName = factor(region_table$GSC_state, levels = levels, labels = labels)
region_table$SCI_stateName = factor(region_table$SCI_state, levels = levels, labels = labels)

## Reduce complexity for similar chromatin state
region_table[region_table == "U3"] = "U2"
region_table[region_table == "U7"] = "U6"
region_table[region_table == "U15"] = "U14"
region_table[region_table == "U17"] = "U16"
region_table[region_table == "U18"] = "U16"
region_table[region_table$GSC_peak == "absent","GSC_state"] =
  paste0(region_table[region_table$GSC_peak == "absent","GSC_state"],"_NO")
region_table[region_table$SCI_peak == "absent","SCI_state"] =
  paste0(region_table[region_table$SCI_peak == "absent","SCI_state"],"_NO")
region_table[region_table$RS_peak == "absent","RS_state"] =
  paste0(region_table[region_table$RS_peak == "absent","RS_state"],"_NO")

region_table$motif_state = paste0(
  region_table$GSC_state,
  region_table$SCI_state,
  region_table$RS_state)

write.table(region_table,"GSC_SC_RS_H3K79me2_state_dynamic_table.tsv", sep = "\t", col.names = T,
  row.names = F, quote = F)

## Plotting data
tableRegion = do.call("rbind",lapply(unique(region_table$motif_state), function(motif){
  df = region_table[region_table$motif_state == motif,c(4:6,8:10,12)]
  nb = nrow(df)
  df = data.frame(motif = motif, freq = nb, grp = str_replace(colnames(df)[4:6],"_state",""),
    state = unlist(df[1,4:6]),H3K79me2 = unlist(df[1,1:3]))
  df$alluvium = "+H3K79me2"
  df[grepl("absent",df[,"H3K79me2"]),"alluvium"] = "-H3K79me2"
  return(df)
}))
tableRegion = tableRegion[tableRegion$freq > 100,]
tableRegion$stateName = factor(tableRegion$state, levels = levels, labels = stateName)

alluvial = ggplot(tableRegion, aes(x = factor(grp, levels = c("GSC","SCI","RS")), y = freq,
    stratum = factor(state, levels = names(colorValue)),
    alluvium = motif, label = as.numeric(freq))) +
  geom_flow(alpha = .2, stat = "alluvium", lode.guidance = "frontback",
    aes(color = state,fill = state)) +
  geom_stratum(alpha = .6, aes(fill = state, color = alluvium), size = .3) +
  scale_fill_manual(values = colorValue) +
  geom_text_repel(aes(label = ifelse(grp == "GSC" & H3K79me2 == "present", as.character(stateName), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5) +
  geom_text_repel(aes(label = ifelse(grp == "RS" & H3K79me2 == "present", as.character(stateName), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5) +
  # geom_text_repel(aes(label = ifelse(grp == "SCI" & H3K79me2 == "present", as.character(stateName), NA)),
  #   stat = "stratum", size = 4, direction = "y", nudge_x = .5) +
  scale_color_manual(values = colorValue) +
  xlab("") + ylab("") + labs(fill = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=15, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none") +
  guides(color = "none")
ggsave(filename = "GSC_SC_RS_H3K79me2_state_dynamic_V4.png",plot = alluvial, width = 6, height = 4,
  device = 'png',dpi = 600)

## Coutn the number of weak became active in SCI
table(region_table[region_table$GSC_peak == "present" &
  region_table$SCI_peak == "present" &
  region_table$GSC_stateName == "weak","SCI_stateName"])

####################################################################################################
##                      CHARACTERIZATION DYNAMIC STATE according to group
####################################################################################################

grp = "het"
region_bivalent = region_table[region_table$RS_stateName == grp,]

region_bivalent = region_table[region_table$GSC_stateName == "het" | region_table$RS_stateName == "het",]

tableRegionBivalent = do.call("rbind",lapply(unique(region_bivalent$motif_state), function(motif){
  df = region_bivalent[region_bivalent$motif_state == motif,c(4:6,8:11)]
  nb = nrow(df)
  df = data.frame(motif = motif, freq = nb, grp = str_replace(colnames(df)[4:6],"_state",""),
    state = unlist(df[1,4:6]),H3K79me2 = unlist(df[1,1:3]))
  df$alluvium = "+H3K79me2"
  df[grepl("absent",df[,"H3K79me2"]),"alluvium"] = "-H3K79me2"
  return(df)
}))
tableRegionBivalent = tableRegionBivalent[tableRegionBivalent$freq > 100,]
tableRegionBivalent$stateName = factor(tableRegionBivalent$state, levels = levels, labels = stateName)

allu = ggplot(tableRegionBivalent, aes(x = factor(grp, levels = c("GSC","SCI","RS")), y = freq,
    stratum = factor(state, levels = names(colorValue)),
    alluvium = motif)) +
  geom_flow(alpha = .3, stat = "alluvium", lode.guidance = "frontback",
    aes(color = state,fill = state)) +
  geom_stratum(alpha = .6, aes(fill = state, color = alluvium), size = .5) +
  scale_fill_manual(values = colorValue) +
  scale_color_manual(values = colorValue) +
  geom_text_repel(aes(label = ifelse(grp == "GSC" & H3K79me2 == "present", as.character(stateName), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5) +
  geom_text_repel(aes(label = ifelse(grp == "RS" & H3K79me2 == "present", as.character(stateName), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5) +
  xlab("") + ylab("") + labs(fill = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=15, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none") +
  guides(color = "none")
ggsave(filename = paste0("GSC_SC_RS_H3K79me2_",grp,"_dynamic.png"),plot = allu, width = 6, height = 4,
  device = 'png',dpi = 600)

####################################################################################################
##                        CHARACTERIZATION DYNAMIC STATE based on second state
####################################################################################################

region_bivalent$pos_peak = paste0(region_bivalent$chr,region_bivalent$start,region_bivalent$end)

region_test = merge(region_bivalent, region, by = "pos_peak")
region_table_test = do.call("rbind",lapply(unique(region_test$pos_peak), function(pos) {
  df = region_test[region_test$pos_peak == pos,]
  df_region = df[1,2:11]
  SC = aggregate(df$SCI_coverage, list(df$SCI_state.y), sum)
  order = SC$x[order(SC$x, decreasing = T)]
  if(length(order) > 1) {
    df_region$SC_state_new = SC[SC$x == order[2],"Group.1"][1]
  } else {df_region$SC_state_new = SC[SC$x == order,"Group.1"][1]}
  return(df_region)
}))

region_table_test[region_table_test == "U3"] = "U2"
region_table_test[region_table_test == "U7"] = "U6"
region_table_test[region_table_test == "U15"] = "U14"
region_table_test[region_table_test == "U17"] = "U16"
region_table_test[region_table_test == "U18"] = "U16"
region_table_test[region_table_test$SCI_peak.x == "absent","SC_state_new"] =
  paste0(region_table_test[region_table_test$SCI_peak.x == "absent","SC_state_new"],"_NO")
region_table_test$motif_state = paste0(
  region_table_test$GSC_state.x,
  region_table_test$SC_state_new,
  region_table_test$RS_state.x)

test_matrix = table(region_table_test$SCI_state.x, region_table_test$SC_state_new)
test_matrix$sum = 1
pheatmap(test_matrix, cluster_rows = F, cluster_cols = F, display_number = T,
  number_format = "%.0f")

tableRegionBivalent_test = do.call("rbind",lapply(unique(region_table_test$motif_state), function(motif){
  df = region_table_test[region_table_test$motif_state == motif,c(4:6,8:11)]
  nb = nrow(df)
  df = data.frame(motif = motif, freq = nb, grp = str_replace(colnames(df)[4:6],"_state",""),
    state = unlist(df[1,4:6]),H3K79me2 = unlist(df[1,1:3]))
  df$alluvium = "+H3K79me2"
  df[grepl("absent",df[,"H3K79me2"]),"alluvium"] = "-H3K79me2"
  return(df)
}))
tableRegionBivalent_test = tableRegionBivalent_test[tableRegionBivalent_test$freq > 100,]

t = ggplot(tableRegionBivalent_test, aes(x = factor(grp, levels = c("GSC.x","SCI.x","RS.x")), y = freq,
    stratum = factor(state, levels = names(colorValue)),
    alluvium = motif)) +
  geom_flow(alpha = .3, stat = "alluvium", lode.guidance = "frontback",
    aes(color = state,fill = state)) +
  geom_stratum(alpha = .6, aes(fill = state, color = alluvium), size = 1) +
  scale_fill_manual(values = colorValue) +
  scale_color_manual(values = colorValue) +
  xlab("") + ylab("") + labs(fill = "") +
  theme_bw() + theme(strip.background  = element_blank(),
    text = element_text(size=15, angle = 0),
    panel.grid.major = element_line(colour = "grey80"),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none") +
  guides(color = "none")
ggsave(filename = "GSC_SC_RS_H3K79me2_bivalent_dynamic_2eme.png",plot = t, width = 8, height = 5,
  device = 'png',dpi = 600)
