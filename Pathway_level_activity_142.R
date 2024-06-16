data("human_gem")
#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux_Young_Spleen)){
    activity_score[d]<-mean(abs(flux_Young_Spleen[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))


#heatmap 
mapal_Young_Spleen <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pheatmap::pheatmap(all_pathway_score,cluster_cols = F,color = rev(mapal_Young_Spleen),scale = "row")



