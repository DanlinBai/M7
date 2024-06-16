# Load Seurat package
library(Seurat)

# View the current working directory
getwd()

# Set the working directory (to switch the working directory to the specified path)
setwd("/rds/homes/d/dxb360/M7/data/")



###########################################################
### 1. Merge cell types by cell name prefix in metadata ###
###########################################################

# read metadata.tsv.gz
metadata <- fread("./GSE210719/metadata.tsv.gz")

classify_cell_type <- function(cell_type) {
  if (grepl("BCell", cell_type)) {
    return("BCell")
  } else if (grepl("Cardiomyocyte", cell_type)) {
    return("Cardiomyocyte")
  } else if (grepl("CD3_", cell_type)) {
    return("CD3")
  } else if (grepl("CD4_", cell_type)) {
    return("CD4")
  } else if (grepl("CD8_", cell_type)) {
    return("CD8")
  } else if (grepl("DNT_", cell_type)) {
    return("DNT")
  } else if (grepl("GammaDelta_", cell_type)) {
    return("GammaDelta")
  } else if (grepl("Macrophage_", cell_type)) {
    return("Macrophage")
  } else if (cell_type == "NKT") {
    return("NKT")
  } else if (cell_type == "SMC") {
    return("SMC")
  } else {
    return("Unassigned")
  } 
}

metadata <- metadata %>%
  mutate(Cell_type = sapply(Graph_Based_Cluster_3, classify_cell_type))

write.table(metadata, file = "/rds/homes/d/dxb360/M7/data/Mouse_to_Human/metadata.tsv.gz", sep = "\t", quote = FALSE, row.names = TRUE)


#################################
### 2. Create a Seurat object ###
#################################

# create v3 assays
options(Seurat.object.assay.version = "v3")

# For reading 10x data, the data.dir parameter specifies the path to the file
seurat_data <- Read10X(data.dir = "./Mouse_to_Human")



seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                 project = "Mouse_to_Human",
                                 min.features = 200, # The minimum number of genes per cell, usually 200 to 500
                                 min.cells = 3) # The minimum number of cells per gene is usually 3 to 10

class(seurat_obj[["RNA"]])

# View basic information about Seurat objects
seurat_obj

########################################################
### 3. Add cell type, group, age columns to metadata ###
########################################################

# read metadata.tsv.gz
metadata <- fread("./Mouse_to_Human/metadata.tsv.gz")

# Suppose the first column of metadata is the cell name
cell_names <- metadata[[1]]
metadata <- metadata[, -1, with = FALSE]

# Make sure the cell name matches the cell name in the expression matrix
rownames(metadata) <- cell_names



# Add the Graph_Based_Cluster_3 information to the Cell_type in the metadata of the Seurat object
seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata[["Cell_type"]], col.name = "Cell_type")

# Add the group and ages in metadata of the Seurat object
seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata[["Group"]], col.name = "Group")
seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata[["Ages"]], col.name = "Ages")

# check metadata
head(seurat_obj@meta.data)

######################################################################
### 4. log normalization is performed on the representation matrix ###
######################################################################

seurat_obj_nor <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# write the human gene into the rownames of the matrix
rownames(seurat_obj_nor@assays[["RNA"]]@counts) <- H_gene_names
rownames(seurat_obj_nor@assays[["RNA"]]@data) <- H_gene_names

# Save the Seurat object
saveRDS(seurat_obj_nor, file = "./METAFlux_data/seurat_object_m2h_cga_nor.rds")


