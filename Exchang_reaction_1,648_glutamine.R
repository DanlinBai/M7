# Read rds files
flux_Young_Spleen <- readRDS("/rds/homes/d/dxb360/M7/data/METAFlux_data/multiple_groups/object1.rds")
flux_Aged_Spleen <- readRDS("/rds/homes/d/dxb360/M7/data/METAFlux_data/multiple_groups/object2.rds")
flux_Young_Aorta <- readRDS("/rds/homes/d/dxb360/M7/data/METAFlux_data/multiple_groups/object3.rds")
flux_Aged_Aorta <- readRDS("/rds/homes/d/dxb360/M7/data/METAFlux_data/multiple_groups/object4.rds")

# The total dimension of predicted flux data 
# (num_cell*13082+1648)*number_of_bootstrap
# Used Human-GEM[1] (consisting of 13082 metabolic reactions and 8378 metabolites) + 1648 exchange reactions
dim(flux_Young_Spleen)
dim(flux_Aged_Spleen)
dim(flux_Young_Aorta)
dim(flux_Aged_Aorta)
# [1] 145550     50
# (11*13082+1648)*50 = 145550*50 

# ‘Nutrient look up files’ contains 1648 exchange reactions that describe how cells uptake or release metabolites.
data("nutrient_lookup_files")

glutamine_Young_Spleen<-data.frame(glutamine=flux_Young_Spleen[grep("HMR_9063",rownames(flux_Young_Spleen)),])
glutamine_Young_Spleen
glutamine_Aged_Spleen<-data.frame(glutamine=flux_Aged_Spleen[grep("HMR_9063",rownames(flux_Aged_Spleen)),])
glutamine_Young_Aorta<-data.frame(glutamine=flux_Young_Aorta[grep("HMR_9063",rownames(flux_Young_Aorta)),])
glutamine_Aged_Aorta<-data.frame(glutamine=flux_Aged_Aorta[grep("HMR_9063",rownames(flux_Aged_Aorta)),])

#if needed we can apply cubic root normalization to normalize the scores
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

glutamine_Young_Spleen_nor = cbrt(glutamine_Young_Spleen)
glutamine_Aged_Spleen_nor = cbrt(glutamine_Aged_Spleen)
glutamine_Young_Aorta_nor = cbrt(glutamine_Young_Aorta)
glutamine_Aged_Aorta_nor = cbrt(glutamine_Aged_Aorta)

####################
### single image ###
####################

# One may look at the distribution of all bootstraps to different cell types nutrient uptake/release profile.
library(ggplot2)
glutamine_Aged_Aorta$celltype=rownames(glutamine_Aged_Aorta)
long_glutamine_Aged_Aorta=reshape2::melt(glutamine_Aged_Aorta,id.vars='celltype')
p <- ggplot(long_glutamine_Aged_Aorta,aes(y=value,x=celltype))+
  geom_boxplot()+
  ggtitle("glutamine uptake/release level for different cell types_Aged_Aorta")+
  xlab("")+ylab("glutamine uptake/release scores")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_AA_box.png", plot = p, width = 16.67, height = 12.5, dpi = 300)


### use the normalizational flux data
glutamine_Aged_Aorta_nor$celltype=rownames(glutamine_Aged_Aorta_nor)
long_glutamine_Aged_Aorta_nor=reshape2::melt(glutamine_Aged_Aorta_nor,id.vars='celltype')
p <- ggplot(long_glutamine_Aged_Aorta_nor,aes(y=value,x=celltype))+
  geom_boxplot()+
  ggtitle("glutamine uptake/release level for different cell types_Aged_Aorta_nor")+ 
  xlab("")+ylab("glutamine uptake/release scores")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_AA_box_nor.png", plot = p, width = 16.67, height = 12.5, dpi = 300)

####################
### Merge images ###
####################
library(reshape2)

### 1. flux data ###
glutamine_Young_Spleen$celltype=rownames(glutamine_Young_Spleen)
glutamine_Young_Aorta$celltype=rownames(glutamine_Young_Aorta)
glutamine_Aged_Spleen$celltype=rownames(glutamine_Aged_Spleen)
glutamine_Aged_Aorta$celltype=rownames(glutamine_Aged_Aorta)

long_glutamine_Young_Spleen=reshape2::melt(glutamine_Young_Spleen,id.vars='celltype')
long_glutamine_Young_Aorta=reshape2::melt(glutamine_Young_Aorta,id.vars='celltype')
long_glutamine_Aged_Spleen=reshape2::melt(glutamine_Aged_Spleen,id.vars='celltype')
long_glutamine_Aged_Aorta=reshape2::melt(glutamine_Aged_Aorta,id.vars='celltype')

# add the information of groups
long_glutamine_Young_Spleen$group <- "Young_Spleen"
long_glutamine_Aged_Spleen$group <- "Aged_Spleen"
long_glutamine_Young_Aorta$group <- "Young_Aorta"
long_glutamine_Aged_Aorta$group <- "Aged_Aorta"

# merge all the data
all_data <- rbind(long_glutamine_Young_Spleen, long_glutamine_Aged_Spleen, long_glutamine_Young_Aorta, long_glutamine_Aged_Aorta)

# add the names of cell type
all_data_celltype <- all_data

Boxplot_cell_type <- function(celltype) {
  if (grepl("celltype 1 ", celltype)) {
    return("BCell")
  } else if (grepl("celltype 2", celltype)) {
    return("CM") # Cardiomyocyte
  } else if (grepl("celltype 3", celltype)) {
    return("CD3")
  } else if (grepl("celltype 4", celltype)) {
    return("CD4")
  } else if (grepl("celltype 5", celltype)) {
    return("CD8")
  } else if (grepl("celltype 6", celltype)) {
    return("DNT")
  } else if (grepl("celltype 7", celltype)) {
    return("γδ") #GammaDelta
  } else if (grepl("celltype 8", celltype)) {
    return("MΦ") #Macrophage
  } else if (grepl("celltype 9", celltype)) {
    return("NKT")
  } else if (grepl("celltype 10", celltype)) {
    return("SMC")
  } else if (grepl("celltype 11", celltype)) {
    return("Unassigned")   
  } else {
    return("External")
  } 
}

all_data_celltype <- all_data_celltype %>%
  mutate(Cell_Type = sapply(celltype, Boxplot_cell_type))

# create boxplot
p <- ggplot(all_data_celltype, aes(x=celltype, y=value, fill=group)) +
  geom_boxplot() +
  facet_wrap(~Cell_Type, scales = "free_x", nrow = 1) +
  #facet_wrap(~group, scales = "free_x", nrow = 1) +
  ggtitle("glutamine uptake/release level for different cell types and groups") +
  xlab("Cell Type") + 
  ylab("glutamine uptake/release scores") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p)

#ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_group_box.png", plot = p, width = 16.67, height = 12.5, dpi = 300)
ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_celltype_box.png", plot = p, width = 16.67, height = 12.5, dpi = 300)



### 2. normaliztional flux data ###
glutamine_Young_Spleen_nor$celltype=rownames(glutamine_Young_Spleen_nor)
glutamine_Young_Aorta_nor$celltype=rownames(glutamine_Young_Aorta_nor)
glutamine_Aged_Spleen_nor$celltype=rownames(glutamine_Aged_Spleen_nor)
glutamine_Aged_Aorta_nor$celltype=rownames(glutamine_Aged_Aorta_nor)

long_glutamine_Young_Spleen_nor=reshape2::melt(glutamine_Young_Spleen_nor,id.vars='celltype')
long_glutamine_Young_Aorta_nor=reshape2::melt(glutamine_Young_Aorta_nor,id.vars='celltype')
long_glutamine_Aged_Spleen_nor=reshape2::melt(glutamine_Aged_Spleen_nor,id.vars='celltype')
long_glutamine_Aged_Aorta_nor=reshape2::melt(glutamine_Aged_Aorta_nor,id.vars='celltype')

# add the information of groups
long_glutamine_Young_Spleen_nor$group <- "Young_Spleen"
long_glutamine_Aged_Spleen_nor$group <- "Aged_Spleen"
long_glutamine_Young_Aorta_nor$group <- "Young_Aorta"
long_glutamine_Aged_Aorta_nor$group <- "Aged_Aorta"

# merge all the data
all_data_nor <- rbind(long_glutamine_Young_Spleen_nor, long_glutamine_Aged_Spleen_nor, long_glutamine_Young_Aorta_nor, long_glutamine_Aged_Aorta_nor)

# add the names of cell type
all_data_celltype_nor <- all_data_nor

all_data_celltype_nor <- all_data_celltype_nor %>%
  mutate(Cell_Type = sapply(celltype, Boxplot_cell_type))

# create boxplot
p <- ggplot(all_data_celltype_nor, aes(x=celltype, y=value, fill=group)) +
  geom_boxplot() +
  #facet_wrap(~Cell_Type, scales = "free_x", nrow = 1) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  ggtitle("glutamine uptake/release level for different cell types and groups_nor") +
  xlab("Cell Type") + 
  ylab("glutamine uptake/release scores") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p)

ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_group_box_nor.png", plot = p, width = 16.67, height = 12.5, dpi = 300)
#ggsave("/rds/homes/d/dxb360/M7/data/METAFlux_data/plots/glutamine_celltype_box_nor.png", plot = p, width = 16.67, height = 12.5, dpi = 300)


