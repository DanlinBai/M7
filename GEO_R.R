install.packages("Matrix")

library(Matrix) 

# read matrix.mtx.gz
expression_matrix <- readMM("/rds/homes/d/dxb360/M7/data/GSE210719/matrix.mtx.gz")

dim(expression_matrix)
# 18388 37260

# read metadata.tsv.gz
metadata <- fread("/rds/homes/d/dxb360/M7/data/GSE210719/metadata.tsv.gz")

# read features.tsv.gz
features <- fread("/rds/homes/d/dxb360/M7/data/GSE210719/features.tsv.gz", header = FALSE)

# read barcodes.tsv.gz
barcodes <- fread("/rds/homes/d/dxb360/M7/data/GSE210719/barcodes.tsv.gz", header = FALSE)

