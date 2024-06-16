# Installed METAFlux R package
devtools::install_github('KChen-lab/METAFlux')

# For single cell data analysis, we provide pipeline to work with Seurat V4
install.packages('Seurat')
packageVersion('Seurat')
# [1] ‘4.3.0.1’
remove.packages("Seurat")
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)

# Update the SeuratObject package to the desired version (>=5.0.2)
install.packages("SeuratObject")
packageVersion('SeuratObject')

install.packages("/rds/homes/d/dxb360/R/x86_64-pc-linux-gnu-library/4.3/Matrix-1.6-4.tar.gz", repos = NULL, type = "source")
packageVersion("Matrix")

######################
#For scRNA-seq sample#
######################

library(METAFlux)
#################
### one group ###
#################
#load the toy seurat example, this data only has one patient. For mutiple samples,see example code at the end of this section.
#data("sc_test_example")
obj.list[["Young_Spleen"]]

# medium file for human derived samples
data("human_blood")
# data("cell_medium")


dim(obj.list[["Young_Spleen"]]) # [1] 13544 18412

#check cell types
table(obj.list[["Young_Spleen"]]$Cell_type)
# BCell Cardiomyocyte           CD3           CD4           CD8           DNT    GammaDelta 
# 410             1           426          5375          7871          3175           364 
# Macrophage           NKT           SMC    Unassigned 
# 311           426             1            52 

#Calculate the mean expression for bootstrap samples from seurat object
#Using "Cell_type" in seurat as my grouping variable
#For the sake of demonstration, we only set the number of bootstraps to 3. 
#In real analysis, the number of bootstraps should be much higher(e.g. 100, 200....)
mean_exp=calculate_avg_exp(myseurat = obj.list[["Young_Spleen"]],myident = 'Cell_type',n_bootstrap=3,seed=1)

write.csv(mean_exp, "/rds/homes/d/dxb360/M7/data/METAFlux_data/one_group/mean_exp_Young_Spleen.csv", row.names = FALSE)

#calculate metabolic reaction scores
scores<-calculate_reaction_score(data=mean_exp)

write.csv(scores, "/rds/homes/d/dxb360/M7/data/METAFlux_data/one_group/scores_Young_Spleen.csv", row.names = FALSE)

#calculate the fractions of celltype/clusters
round(table(obj.list[["Young_Spleen"]]$Cell_type)/nrow(obj.list[["Young_Spleen"]]@meta.data),4)

# BCell Cardiomyocyte           CD3           CD4           CD8           DNT    GammaDelta 
# 0.0223        0.0001        0.0231        0.2919        0.4275        0.1724        0.0198 
# Macrophage           NKT           SMC    Unassigned 
# 0.0169        0.0231        0.0001        0.0028 


# calculate flux: here we used human blood as our medium, but please use cell line medium when samples are cell line.

#num_cell=number of cell types/clusters, here we have 11 cell types, thus input is 11. Make sure the sum of cell percentage is 1.The order of fraction should be the same as that of "Mean_exp" data
flux=compute_sc_flux(num_cell = 11,fraction =c(0.0223,0.0001,0.0231,0.2919,0.4275,0.1724, 0.0198,0.0169,0.0231,0.0001,0.0028),fluxscore=scores,medium = human_blood)

write.csv(flux, "/rds/homes/d/dxb360/M7/data/METAFlux_data/one_group/flux_Young_Spleen.csv", row.names = FALSE)

#optional: flux scores normalization
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux_nor = cbrt(flux)

write.csv(flux_nor, "/rds/homes/d/dxb360/M7/data/METAFlux_data/one_group/flux_nor_Young_Spleen.csv", row.names = FALSE)
######################
### multiple groups###
######################

#load your seurat object
seurat_obj_m2h_cga_nor <- readRDS("./METAFlux_data/seurat_object_m2h_cga_nor.rds")

seurat_obj_m2h_cga_nor
Layers(object = seurat_obj_m2h_cga_nor)


#sample code for 4 clusters in multiple patients
obj.list <- SplitObject(seurat_obj_m2h_cga_nor, split.by = "Group")

saveRDS(obj.list, file = "./METAFlux_data/obj_list.rds")


for (i in c(1:length(obj.list))){
  sc<-obj.list[[i]]
  mean_exp=calculate_avg_exp(myseurat = sc,myident = 'Cell_type',n_bootstrap=50,seed=1)
  scores<-calculate_reaction_score(data=mean_exp)
  g=round(table(sc$Cell_type)/nrow(sc@meta.data),4)
  # g=table(sc$Cell_type)/nrow(sc@meta.data)
  
  # Calculate the rounded sum
  sum_rounded <- sum(g)

  # Calculated adjustment
  difference <- 1 - sum_rounded
  
  # If the sum does not equal 1, adjust
  if (difference != 0) {
    # If the difference is greater than 0, increase the minimum value; Otherwise, reduce the maximum value
    if (difference > 0) {
      index <- which.min(g)
    } else {
      index <- which.max(g)
    }
    g[index] <- g[index] + difference
  }
  
  print(g)
  flux=compute_sc_flux(num_cell = 11,fraction =c(g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11]),fluxscore=scores,medium = human_blood)
  saveRDS(flux,paste0("/rds/homes/d/dxb360/M7/data/METAFlux_data/multiple_groups/object",i,".rds"))
}



