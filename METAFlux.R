# Installed METAFlux R package
devtools::install_github('KChen-lab/METAFlux')

# For single cell data analysis, we provide pipeline to work with Seurat V4
install.packages('Seurat')
packageVersion('Seurat')
# [1] ‘4.3.0.1’
remove.packages("Seurat")


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
##################
### one patient###
##################
#load the toy seurat example, this data only has one patient. For mutiple samples,see example code at the end of this section.
data("sc_test_example")
# medium file for human derived samples
data("human_blood")
# data("cell_medium")

dim(sc_test_example)

#check cell types
table(sc_test_example$Cell_type)
#B lymphocytes Epithelial cells    Myeloid cells    T lymphocytes 
#     37              120              120              120

#Calculate the mean expression for bootstrap samples from seurat object
#Using "Cell_type" in seurat as my grouping variable
#For the sake of demonstration, we only set the number of bootstraps to 3. 
#In real analysis, the number of bootstraps should be much higher(e.g. 100, 200....)
mean_exp=calculate_avg_exp(myseurat = sc_test_example,myident = 'Cell_type',n_bootstrap=3,seed=1)

write.csv(mean_exp, "/rds/homes/d/dxb360/M7/data/METAFlux_demo/mean_exp.csv", row.names = FALSE)

#calculate metabolic reaction scores
scores<-calculate_reaction_score(data=mean_exp)

write.csv(scores, "/rds/homes/d/dxb360/M7/data/METAFlux_demo/scores.csv", row.names = FALSE)

#calculate the fractions of celltype/clusters
round(table(sc_test_example$Cell_type)/nrow(sc_test_example@meta.data),1)
#   B lymphocytes Epithelial cells    Myeloid cells    T lymphocytes 
#       0.1              0.3              0.3              0.3 

# calculate flux: here we used human blood as our medium, but please use cell line medium when samples are cell line.

#num_cell=number of cell types/clusters, here we have 4 cell types, thus input is 4. Make sure the sum of cell percentage is 1.The order of fraction should be the same as that of "Mean_exp" data
flux=compute_sc_flux(num_cell = 4,fraction =c(0.1,0.3,0.3,0.3),fluxscore=scores,medium = human_blood)

write.csv(flux, "/rds/homes/d/dxb360/M7/data/METAFlux_demo/flux.csv", row.names = FALSE)

#optional: flux scores normalization
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux_nor = cbrt(flux)

write.csv(flux_nor, "/rds/homes/d/dxb360/M7/data/METAFlux_demo/flux_nor.csv", row.names = FALSE)

########################
### multiple patients###
########################

#load your seurat object
#sample code for 4 clusters in multiple patients
obj.list <- SplitObject(seurat, split.by = "patient_id")
for (i in c(1:length(obj.list))){
  sc<-obj.list[[i]]
  mean_exp=calculate_avg_exp(myseurat = sc_test_example,myident = 'Cell_type',n_bootstrap=50,seed=1)
  scores<-calculate_reaction_score(data=mean_exp)
  #g=round(table(sc$Cell_type)/nrow(sc@meta.data),3)
  g=table(sc$Cell_type)/nrow(sc@meta.data)
  print(g)
  flux=compute_sc_flux(num_cell = 4,fraction =c(g[1],g[2],g[3],g[4]),fluxscore=scores,medium = human_blood)
  saveRDS(flux,paste0("object",i,".rds"))
}

