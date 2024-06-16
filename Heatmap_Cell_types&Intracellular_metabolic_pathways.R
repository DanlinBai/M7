# Heat map of cell types and intracellular metabolic pathways
#######################################################################################
### 1. Delete rows of data for extracellular exchange reactions in the flux results ###
#######################################################################################
# Use grep to find the line index that contains the line name "external_medium"
rows_to_remove <- grep("external_medium", rownames(flux_Young_Spleen))
rows_to_remove
# Remove these rows
flux_Young_Spleen_internal <- flux_Young_Spleen[-rows_to_remove, ]

rows_to_remove <- grep("external_medium", rownames(flux_Young_Aorta))
flux_Young_Aorta_internal <- flux_Young_Aorta[-rows_to_remove, ]

rows_to_remove <- grep("external_medium", rownames(flux_Aged_Spleen))
flux_Aged_Spleen_internal <- flux_Aged_Spleen[-rows_to_remove, ]

rows_to_remove <- grep("external_medium", rownames(flux_Aged_Aorta))
flux_Aged_Aorta_internal <- flux_Aged_Aorta[-rows_to_remove, ]

#############################################################################
### 2. calculate the mean of 50 possible solutions for each reaction flux ###
#############################################################################
# Calculate the mean of each row
YS_row_means <- apply(flux_Young_Spleen_internal, 1, mean)
YA_row_means <- apply(flux_Young_Aorta_internal, 1, mean)
AS_row_means <- apply(flux_Aged_Spleen_internal, 1, mean)
AA_row_means <- apply(flux_Aged_Aorta_internal, 1, mean)

# Create a new data box that contains only the mean of each row
YS_internal_flux_means <- data.frame(row_mean = YS_row_means)
YA_internal_flux_means <- data.frame(row_mean = YA_row_means)
AS_internal_flux_means <- data.frame(row_mean = AS_row_means)
AA_internal_flux_means <- data.frame(row_mean = AA_row_means)

# Create an empty data box to hold the new columns
internal_flux_means_11_celltypes <- data.frame(matrix(nrow = 13082, ncol = 11))

# Set a new column name
colnames(internal_flux_means_11_celltypes) <- paste0("celltype", 1:11)

# Fill data frame
for (i in 1:11) {
  celltype_rows <- grepl(paste0("celltype ", i," "), rownames(AA_internal_flux_means))
  internal_flux_means_11_celltypes[, i] <- AA_internal_flux_means$row_mean[celltype_rows]
}

# The filled row name are the intracellular reactions' ID
row.names(internal_flux_means_11_celltypes) <- human_gem$ID

# Use specific cell types as column names
new_colnames <- c("BCell", "CM", "CD3", "CD4", "CD8", "DNT", "γδ", "MΦ", "NKT", "SMC", "Unassigned")
colnames(internal_flux_means_11_celltypes) <- new_colnames

################################################################################################
### 3. create the heat map to find the relationship between cell types and metabolic pathway ###
################################################################################################
# Modify according to the script Pathway_level_activity_142.R
#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(internal_flux_means_11_celltypes)){
    activity_score[d]<-mean(abs(internal_flux_means_11_celltypes[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))
colnames(all_pathway_score) <- c("BCell", "CM", "CD3", "CD4", "CD8", "DNT", "γδ", "MΦ", "NKT", "SMC", "Unassigned")

#heatmap 
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pheatmap::pheatmap(all_pathway_score, cluster_cols = FALSE, color = rev(mapal), scale = "row",
                   cellwidth = 15, cellheight = 8)  












