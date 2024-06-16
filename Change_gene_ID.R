# The bioMart database can be freely connected and various genetic transformations can be performed  
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

# View the available BioMart database
listMarts()
#                biomart                version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 112
# 2   ENSEMBL_MART_MOUSE      Mouse strains 112
# 3     ENSEMBL_MART_SNP  Ensembl Variation 112
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 112


# Lists the data sets available in the selected BioMart
listDatasets(useMart("ENSEMBL_MART_ENSEMBL"))
#                            dataset
# 1     abrachyrhynchus_gene_ensembl
# 2         acalliptera_gene_ensembl
# 3       acarolinensis_gene_ensembl
# 4        acchrysaetos_gene_ensembl
#                             ......
# 79             hhucho_gene_ensembl
# 80           hsapiens_gene_ensembl
#                             ......
# 106          mmurinus_gene_ensembl
# 107         mmusculus_gene_ensembl
#                             ......

# Select the BioMart database you want to use and use the useMart function to create Mart object
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Use the lsitFilters() function to view the types available, select the type of annotation to retrieve, and the type of known annotation
listFilters(mart)
#                                      name
# 1                         chromosome_name
# 2                                   start
# 3                                     end
# 4                              band_start
# 5                                band_end
# 6                            marker_start

# The first step is to build the mart object
# Use the useMart function to connect to the specified BioMart database and data sets in the database
human <- useMart('ensembl', dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl', dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

listDatasets(human)
# 80           hsapiens_gene_ensembl                                     Human genes (GRCh38.p13)
listDatasets(mouse)
# 107         mmusculus_gene_ensembl                                         Mouse genes (GRCm39)

# Mapping mouse genes onto human genes
mouse.gene <- features[[2]]


m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                values = mouse.gene,mart = mouse,
                attributesL = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position"),
                martL = human,uniqueRows = T)
### attributes ### This represents the attribute parameters of the dataset we are trying to retrieve. 
# mgi_symbol represents the symbol name of the mouse gene. 
# Retrieve a list of possible attributes using the listAttributes function.
head(listAttributes(mouse))
#                            name                  description         page
# 1               ensembl_gene_id               Gene stable ID feature_page
# 2       ensembl_gene_id_version       Gene stable ID version feature_page
# 3         ensembl_transcript_id         Transcript stable ID feature_page
# 4 ensembl_transcript_id_version Transcript stable ID version feature_page
# 5            ensembl_peptide_id            Protein stable ID feature_page
# 6    ensembl_peptide_id_version    Protein stable ID version feature_page

### filter ### Parameter filter that should be used in the query. These filters are applied to the master data set.
### value ### represents the dataset we want to enter.
### mart ### refers to the mart object of input data, that is, the gene of the mouse
### attributesL ### represents another database that we need for homology transformation
### useMartL ### The parameter represents the Mart object we need to link to, which is human

head(m2h.g)
# MGI.symbol  Gene.stable.ID HGNC.symbol Chromosome.scaffold.name Gene.start..bp.
# 1    mt-Atp8 ENSG00000228253     MT-ATP8                       MT            8366
# 2    mt-Atp6 ENSG00000198899     MT-ATP6                       MT            8527
# 3     mt-Co2 ENSG00000198712      MT-CO2                       MT            7586
# 4     mt-Nd4 ENSG00000198886      MT-ND4                       MT           10760
# 5    mt-Cytb ENSG00000198727      MT-CYB                       MT           14747
# 6     mt-Nd5 ENSG00000198786      MT-ND5                       MT           12337

# read features.tsv
features_mapping <- fread("/rds/homes/d/dxb360/M7/data/GSE210719/features.tsv.gz", header = FALSE)

# Creates a new column and initializes it as NA
features_mapping$Human_GeneID <- NA
features_mapping$Human_GeneName <- NA

# Fill new column
for (i in 1:nrow(features_mapping)) {
  mouse_gene_name <- features_mapping$V2[i]
  # Find the corresponding human gene ID and gene name in the gene mapping file
  match_idx <- match(mouse_gene_name, m2h.g$MGI.symbol)
  if (!is.na(match_idx)) {
    features_mapping$Human_GeneID[i] <- m2h.g$Gene.stable.ID[match_idx]
    features_mapping$Human_GeneName[i] <- m2h.g$HGNC.symbol[match_idx]
  }
}

# Save as a new features mapping.tsv file
write.table(features_mapping, file = "/rds/homes/d/dxb360/M7/data/METAFlux_data/features_mapping.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# read expression matrix
expression_matrix_new <- readMM("/rds/homes/d/dxb360/M7/data/GSE210719/matrix.mtx.gz")
head(expression_matrix_new)

# Write the mouse gene name to the matrix name
gene_names<- features_mapping$V2
rownames(expression_matrix_new) <- gene_names
head(expression_matrix_new)

# Replace the content that does not contain letters or numbers with NA
features_mapping$Human_GeneName <- ifelse(grepl("[a-zA-Z0-9]", features_mapping$Human_GeneName), features_mapping$Human_GeneName, NA)

# Extrat the mouse gene names matching human gene as NA
# Extract the control information using a Boolean index
extracted_NA <- features_mapping$V2[is.na(features_mapping$Human_GeneName)]

# Removes the corresponding rows in the matrix and feature
mat_filtered <- expression_matrix_new[!(rownames(expression_matrix_new) %in% extracted_NA), ]
features_filtered <- features_mapping[!is.na(features_mapping$Human_GeneName), ]

# Write the human gene name to the matrix name
H_gene_names<- features_filtered$Human_GeneName
write.csv(H_gene_names, file = "/rds/homes/d/dxb360/M7/data/METAFlux_data/H_gene_names.csv", row.names = FALSE)

rownames(mat_filtered) <- H_gene_names
head(mat_filtered)


# Save as new files to cretate new Seurat Object
write.table(features_filtered, file = "/rds/homes/d/dxb360/M7/data/Mouse_to_Human/features.tsv.gz", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeMM(as(mat_filtered, "sparseMatrix"), "matrix.mtx",file = "/rds/homes/d/dxb360/M7/data/Mouse_to_Human/matrix.mtx.gz")
