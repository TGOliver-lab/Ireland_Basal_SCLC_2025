# WT Organoids versus RPM Basal Organoids and RPM Basal Allografts
# Ireland et al, 2025
# Related to *Old Fig. 2d-i and Extended Data Fig. 4a-f*
# Related to *Final Fig. 1i-n and Extended Data Fig. 4a-f*

# Load necessary packages
library(Seurat)
library(SeuratObject)
library(zellkonverter)
library(devtools)
library(SummarizedExperiment)

setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat/scanpy")
# If using Isilon share
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig2_ExtFig4_RPM_Organoids_Allo")

################################################################################################################################################################
### First, read in anndata of basal organoids (clustered WT and RPM transformed) for Extended Data Fig 4a-c
## Anndata object generated as in https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig2d_ExtFig4a-d_WT_RPM_Organoids_Allograft_Final_Clean.ipynb
################################################################################################################################################################

# Read in adata object as SCE 
adata<-readH5AD("RPM_organoids_only/050125_RPM_WT_CTpostCre_CTpreCre_TBOs_adata2.h5ad")

# Check object metadata
# "CT"= CellTag information

table(adata$Cre)
# CT_postCre  CT_preCre     No_Cre 
# 7196       1248       7393 

table(adata$leiden_scVI_1.2)
# 0    1    2    3    4    5    6    7 
# 3917 3420 2458 1819 1681 1270 1031  241 

table(adata$Batch)
# Org_CTpostCre  Org_CTpreCre    Org_No_Cre 
# 7196          1248          7393 


dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
length(rownames(adata))
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)

# Ensure sample barcodes are unique 
colnames(adata) <- make.unique(colnames(adata))

#Convert SCE to seurat
RPM_Orgs <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_Orgs
# An object of class Seurat 
# 55491 features across 15837 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

DefaultAssay(RPM_Orgs)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_Orgs[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_Orgs)<-'norm'

# Add embeddings of umap, and X_scVI_1.2 to assay norm
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_Orgs)
dim(test)
RPM_Orgs[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_Orgs[['umap']]@cell.embeddings

# Define color vector and use to generate UMAPs in Ext. Data Fig 4a

colors<-c('deepskyblue4', '#FC8D62','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#0aa6d8', 'lightblue','salmon1','lightgreen')

# Ext. Data Fig. 4a #

RPM_Orgs@meta.data$Sample<-factor(RPM_Orgs@meta.data$Batch,c("Org_No_Cre","Org_CTpreCre","Org_CTpostCre"))
table(RPM_Orgs@meta.data$Sample)
# Org_No_Cre  Org_CTpreCre Org_CTpostCre 
# 7393          1248          7196 

sample_cols<-c("orange","#D8BFD8","darkorchid4", # deep teal
               "turquoise", "#1B9E77" # dark orange
)

DimPlot(RPM_Orgs,group.by='Sample',cols=sample_cols, reduction='umap',label=FALSE,label.size=6, shuffle=TRUE)&NoAxes()
DimPlot(RPM_Orgs,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=TRUE,label.size=10)&NoAxes()


# Normalize and log counts for downstream signatures/visualization etc
DefaultAssay(RPM_Orgs)<-"RNA"
RPM_Orgs<-NormalizeData(RPM_Orgs)

## Perform Cell Cycle Analysis for Ext Data Fig 4c and visualize ##
# lowercase phase is from python, upper will be from seurat
library(devtools)
library(Seurat)

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

# Read in this table to convert human genes to mouse homologs
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")

s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse


DefaultAssay(RPM_Orgs)
RPM_Orgs <- CellCycleScoring(RPM_Orgs, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPM_Orgs@meta.data$Phase<-factor(RPM_Orgs@meta.data$Phase,levels=c("G1","S","G2M"))

# UMAP in Ext. Data. Fig. 4c
DimPlot(RPM_Orgs, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.05,cols=c("indianred3","green3","royalblue4"))+NoAxes()


# Group to plot bar graph in Ext. Data. Fig. 4c
RPM_Orgs@meta.data$Group<- ifelse(RPM_Orgs@meta.data$Sample %in% c("Org_No_Cre"), "WT_Organoid",
                                       ifelse(RPM_Orgs@meta.data$Sample %in% c("Org_CTpreCre","Org_CTpostCre"), "RPM_Organoid", "Other"))





RPM_Orgs@meta.data$Group<-factor(RPM_Orgs@meta.data$Group,c("WT_Organoid","RPM_Organoid"))
table(RPM_Orgs@meta.data$Group)
# WT_Organoid RPM_Organoid 
# 7393         8444 

# Plot bar graph in Ext. Data. Fig. 4c #
##### What % of cells per phase for WT vs RPM transformed organoids ####
RPM_Orgs@meta.data$Group
Idents(RPM_Orgs)<-'Group'
Idents(RPM_Orgs)

x<-table(Idents(RPM_Orgs),RPM_Orgs@meta.data$Phase)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")


ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20,angle=45, hjust = 1), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


saveRDS(RPM_Orgs,"05_2025_RPM_Orgs_Only_ExtFig4a-c.rds")




################################################################################################################################################################
## Next, cluster WT and transformed RPM basal organoids with resulting RPM allograft tumors (Allo1, Allo3) for Fig. 1i and Ext. Data. Fig. 4d ###
## Anndata object previously generated in Scanpy, https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig2d_ExtFig4a-d_WT_RPM_Organoids_Allograft_Final_Clean.ipynb
################################################################################################################################################################

# Read in adata object as SCE 
adata<-readH5AD("RPM_allo3_wt_transformed_TBOs/050325_RPM_Orgs_Allos_Updated_adata3.h5ad")

# Check anndata metadata to ensure correct dataset
table(adata$Cre)
# CT_postCre  CT_preCre     No_Cre 
# 14498       5108       6678 

table(adata$UnID)
# RPM_Allo_3_New RPM_Allo_Original       RPM_Org_Cre      WT_Org_NoCre 
# 7510              3890              8206              6678 

table(adata$Model)
# Allograft  Organoid 
# 11400     14884 

# Add count data as an assay 
dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
length(rownames(adata))
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)


# Ensure sample barcodes are unique 
colnames(adata) <- make.unique(colnames(adata))

#Convert SCE to seurat
RPM_Orgs_Allo <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_Orgs_Allo
# n object of class Seurat 
# 55491 features across 26284 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts
DefaultAssay(RPM_Orgs_Allo)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_Orgs_Allo[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_Orgs_Allo)<-'norm'

# Add embeddings of umap to assay norm
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_Orgs_Allo)
dim(test)
RPM_Orgs_Allo[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_Orgs_Allo[['umap']]@cell.embeddings

# Define color scheme and generate UMAP in Fig. 1i #
RPM_Orgs_Allo@meta.data$Sample<-factor(RPM_Orgs_Allo@meta.data$Batch,c("Org_No_Cre","Org_CTpreCre","Org_CTpostCre","RPM_Allo_Original","RPM_Allo_3_New"))
table(RPM_Orgs_Allo@meta.data$Sample)
# Org_No_Cre      Org_CTpreCre     Org_CTpostCre RPM_Allo_Original    RPM_Allo_3_New 
# 6678              1218              6988              3890              7510 

sample_cols<-c("orange","#D8BFD8","darkorchid4", # deep teal
               "turquoise", "#1B9E77" # dark orange
               )

# UMAP in Fig. 1i #
DimPlot(RPM_Orgs_Allo,group.by='Sample',cols=sample_cols, reduction='umap',label=FALSE,label.size=6, shuffle=TRUE, pt.size=0.01)&NoAxes()


# Normalize and log counts for downstream signatures/visualization etc
DefaultAssay(RPM_Orgs_Allo)<-"RNA"
RPM_Orgs_Allo<-NormalizeData(RPM_Orgs_Allo)

################## Define Groups by WT, RPM organoid, and RPM allograft to compare Cell Cycle ###############
RPM_Orgs_Allo@meta.data$Group<- ifelse(RPM_Orgs_Allo@meta.data$Sample %in% c("Org_No_Cre"), "WT_Organoid",
                                       ifelse(RPM_Orgs_Allo@meta.data$Sample %in% c("Org_CTpreCre","Org_CTpostCre"), "RPM_Organoid",
                                              ifelse(RPM_Orgs_Allo@meta.data$Sample %in% c("RPM_Allo_Original","RPM_Allo_3_New"), "RPM_Allograft", "Other")))





RPM_Orgs_Allo@meta.data$Group<-factor(RPM_Orgs_Allo@meta.data$Group,c("WT_Organoid","RPM_Organoid","RPM_Allograft"))
table(RPM_Orgs_Allo@meta.data$Group)
# WT_Organoid  RPM_Organoid RPM_Allograft 
# 6678          8206         11400 


## Perform Cell Cycle Analysis for Ext Data Fig 4d and visualize ##
########################################################################
# lowercase phase is from python, upper will be from seurat
library(devtools)
library(Seurat)

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

# Read in this table to convert human genes to mouse homologs
mouse94<-read.csv("Signatures/mouse94.csv")
s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse


DefaultAssay(RPM_Orgs_Allo)
RPM_Orgs_Allo <- CellCycleScoring(RPM_Orgs_Allo, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPM_Orgs_Allo@meta.data$Phase<-factor(RPM_Orgs_Allo@meta.data$Phase,levels=c("G1","S","G2M"))

# UMAP in Ext. Data Fig. 4d 
DimPlot(RPM_Orgs_Allo, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.05,cols=c("indianred3","green3","royalblue4"))+NoAxes()

# Plot barplot for Ext. Data Fig. 4d #
##### What % of cells per phase in each group? ####
Idents(RPM_Orgs_Allo)<-'Group'


x<-table(Idents(RPM_Orgs_Allo),RPM_Orgs_Allo@meta.data$Phase)
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20,angle=45, hjust = 1), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))

# Add basal and normal NE cell signature for Fig. 1i #

######## Cell Type signatures ###########
ct<-read.csv("Signatures/Montoro_Consensus.csv")
basal_noC<-ct$Basal_noC[1:41]
basal_noC
ne_noC<-ct$Neuroendocrine_noC[1:32]
ne_noC

RPM_Orgs_Allo<-AddModuleScore(
  object = RPM_Orgs_Allo,
  features = list(basal_noC),
  name = 'Basal_Consensus')

RPM_Orgs_Allo<-AddModuleScore(
  object = RPM_Orgs_Allo,
  features = list(ne_noC),
  name = 'NE_Consensus')


library(viridis)
library(ggplot2)

FeaturePlot(RPM_Orgs_Allo, features = c("NE_Consensus1"), pt.size=0.2, reduction='umap', order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(RPM_Orgs_Allo, features = c("Basal_Consensus1"), pt.size=0.2, reduction='umap',order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()


saveRDS(RPM_Orgs_Allo,"05_2025_RPM_Orgs_Allo_Fig2d.rds")


####################################################################################################
########### Last, analyze RPM allograft tumour cells only for Fig. 1j #######################
## FEWER CELLS BC DOWNSAMPLED TO KEEP HETEROGENEITY/representation of rare clusters like the P cluster #######
## Anndata object originally generated in scanpy/scVI here: https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig2e_ExtFig4e_RPM_Allos_Final_Clean.ipynb
####################################################################################################

# Read in adata object as SCE 
adata<-readH5AD("RPM_allograft_only/042725_RPM_TBOAllos_OriginalandAllo3_adata3_5kHVG_subsamplebycluster.h5ad")

# Check metadata to ensure correct dataset

table(adata$Cre)
# CT_post-Cre  CT_pre-Cre 
# 2836        1599 

table(adata$leiden_scVI_1.3)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 387 366 337 325 295 286 282 280 279 273 236 153 149 145 134 118 111 104  84  80  11 


dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)
adata

#Convert SCE to seurat
RPM_Allo <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_Allo
# An object of class Seurat 
# 55491 features across 4435 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

DefaultAssay(RPM_Allo)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_Allo
RPM_Allo[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_Allo)<-'norm'

# Add embeddings of umap, and X_scVI_1.3 to assay norm
adata

reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_Allo)
dim(test)
RPM_Allo[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_Allo[['umap']]@cell.embeddings


# Define color scheme and plot UMAP as in Fig. 1j #
my_colors <- c(
  "#E41A1C", # strong red
  "#377EB8", # medium blue
  "#4DAF4A", # green
  "#984EA3", # purple
  "#FF7F00", # orange
  "#FFFF33", # yellow
  "#A65628", # brown
  "#e7298a", # pink
  "#666666", # grey
  "lavender", # teal
  "#FC8D62", # salmon
  "#8DA0CB", # soft blue
  "#E78AC3", # soft pink (different from 8)
  "#A6D854", # light green (but yellowish tint, not green)
  "#FFD92F", # lemon yellow
  "#E5C494", # light brown
  "#B3B3B3", # light grey
  "#1B9E77", # deep teal
  "#D95F02", # dark orange
  "#7570B3", # strong purple
  "turquoise"  # olive green (NOT same green as before)
)

# Fig. 1j UMAP #
DimPlot(RPM_Allo,group.by='leiden_scVI_1.3',cols=my_colors, reduction='umap',label=TRUE,label.size=7)&NoAxes()

# Can visualize UMAP by individual sample
# Rename metadata with more experimental details
RPM_Allo@meta.data$Sample<- ifelse(RPM_Allo@meta.data$UnID %in% c("RPM_Allo3"), "New_RPM_CTpostCre",
                                       ifelse(RPM_Allo@meta.data$UnID %in% c("RPM_Allo_New"), "OG_RPM_CTpreCre", "Other"))


RPM_Allo@meta.data$Sample<-factor(RPM_Allo@meta.data$Sample, c("OG_RPM_CTpreCre","New_RPM_CTpostCre"))

# 'CT' = CellTag
table(RPM_Allo@meta.data$Sample)
# OG_RPM_CTpreCre New_RPM_CTpostCre 
# 1599              2836 

DimPlot(RPM_Allo,group.by='Sample',cols=c("lavender","darkorchid4"), reduction='umap',label=FALSE,label.size=7)&NoAxes()


######## Before applying gene signatures, normalize and log transform  count data ############
DefaultAssay(RPM_Allo)<-'RNA'
RPM_Allo<-NormalizeData(RPM_Allo)

####### Visualize A, N, P, Atoh1, Y, and P63 for Fig. 1j ##############

VlnPlot(
  RPM_Allo,
  features = c("Ascl1", "Neurod1", "Pou2f3", "Atoh1", "Yap1", "Trp63"),
  group.by = "leiden_scVI_1.3",
  cols = my_colors,
  alpha = 0.7,
  ncol = 6
) & 
  theme(
    axis.text.x = element_text(size = 8),                # smaller x-axis labels
    plot.title = element_text(face = "italic", size=30),          # italicize titles
    axis.title.x = element_text(size = 24),              # optional size tweaks
    axis.title.y = element_text(size = 24)
  ) &
  labs(
    y = "Expression",
    x = "Cluster"
  )

# Example code to perform Kruskal-Wallis tests on genes in violin plots in Fig. 2e #
library(dplyr)
library(rstatix)   # Optional, nice formatting of results


# Extract the data
expr_values <- FetchData(RPM_Allo, vars = c("Ascl1", "leiden_scVI_1.3"))

# Rename columns for convenience
colnames(expr_values) <- c("expression", "group")

# Run Kruskal-Wallis test
kruskal_result <- kruskal.test(expression ~ group, data = expr_values)
print(kruskal_result)

########################################################################
# Assign states for Fig. 1k
########################################################################
# Tuft- 20
# NE- 1, 2, 3, 5, 9, 11, 14, 16, 17, 18, 
# NE/Neuronal-0, 8, 15, 19
# Neuronal-10, 12, 13, 
# Atoh1-4, 7
# Basal-6

RPM_Allo@meta.data$Pheno<- ifelse(RPM_Allo@meta.data$leiden_scVI_1.3 %in% c("20"), "Tuft",
                                  ifelse(RPM_Allo@meta.data$leiden_scVI_1.3 %in% c("6"), "Basal",
                                         ifelse(RPM_Allo@meta.data$leiden_scVI_1.3 %in% c("0","8","15","19"), "NE_Neuronal",
                                                ifelse(RPM_Allo@meta.data$leiden_scVI_1.3 %in% c("1","2","3","5","9","11","14","16","17","18"), "NE",
                                                       ifelse(RPM_Allo@meta.data$leiden_scVI_1.3 %in% c("10","12","13"), "Neuronal", "Atoh1")))))




table(RPM_Allo@meta.data$Pheno)
# Atoh1       Basal          NE NE_Neuronal    Neuronal        Tuft 
# 655         282        2173         784         530          11 


RPM_Allo@meta.data$Pheno<-factor(RPM_Allo@meta.data$Pheno, levels=c("NE","NE_Neuronal","Neuronal","Atoh1","Tuft","Basal"))
pheno_col<-c("brown2","darkorchid4","dodgerblue","#66A61E","orange","turquoise4")

## UMAP by Fate in Fig. 1k #
DimPlot(RPM_Allo, group.by=c("Pheno"), cols=pheno_col, shuffle=TRUE, pt.size=0.6)+NoAxes()


######### Add Gene Signature Data ##################

## Look at ND1 vs ASCL1 vs POU2F3 vs ATOH1 ChIP targets as scores (Fig. 1n)

chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ATOH_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
atchip<-chip$ATOH1_Targets
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_mouse<-(pchip$gene_name)
pchip_mouse
atoh<-subset(mouse94, mouse94$human_homolog %in% atchip)
atoh_mouse<-(atoh$gene_name)
atoh_mouse


RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(achip),
  name = 'ASCL1_Targets')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(pchip_mouse),
  name = 'POU2F3_Targets')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(atoh_mouse),
  name = 'ATOH1_Targets')


## Assign YAP1 activity score in lieu of YAP1 target genes ## (Fig. 1n)

yap<-read.csv("Signatures/Yap-Taz_signatures.csv")
yap

y_activ<-yap$mYapTaz_Wang_genes
y_activ<-y_activ[1:22]


RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(y_activ),
  name = 'YAP1_Activity_Score')



######## Normal Cell Type signatures (Fig 1m) ###########
ct<-read.csv("Signatures/Montoro_Consensus.csv")
basal_noC<-ct$Basal_noC[1:41]
basal_noC
ne_noC<-ct$Neuroendocrine_noC[1:32]
ne_noC
tuft_noC<-ct$Tuft_noC[1:60]
tuft_noC

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(tuft_noC),
  name = 'Tuft_Consensus')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(basal_noC),
  name = 'Basal_Consensus')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(ne_noC),
  name = 'NE_Consensus')


## Look at SCLC archetype signatures by Leiden cluster (Ext Data Fig 4f) ##
arch<-read.csv("Signatures/Archetype_Sigs_Maddox.csv")
a<-arch$SCLC.A
a2<-arch$SCLC.A2
a2
n<-arch$SCLC.N
p<-arch$SCLC.P
y<-arch$SCLC.Y

a<-a[1:987]
a_sc<-subset(mouse94, mouse94$human_homolog %in% a)
a_sc_mouse<-(a_sc$gene_name)
a_sc_mouse

n_sc<-subset(mouse94, mouse94$human_homolog %in% n)
n_sc_mouse<-(n_sc$gene_name)
n_sc_mouse

a2
a2_sc<-subset(mouse94, mouse94$human_homolog %in% a2)
a2_sc_mouse<-(a2_sc$gene_name)
a2_sc_mouse

p
p_sc<-subset(mouse94, mouse94$human_homolog %in% p)
p_sc_mouse<-(p_sc$gene_name)
p_sc_mouse

y
y_sc<-subset(mouse94, mouse94$human_homolog %in% y)
y_sc_mouse<-(y_sc$gene_name)
y_sc_mouse


RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')


# Violin plots of Archetypes by leiden cluster for Ext. Data Fig. 4f #
VlnPlot(
  RPM_Allo,
  features = c("A_Archetype1","A2_Archetype1","N_Archetype1","P_Archetype1"),
  group.by = "leiden_scVI_1.3",
  cols = my_colors,
  alpha = 0.7,
  ncol = 2
) & 
  theme(
    axis.text.x = element_text(size = 8),                # smaller x-axis labels
    plot.title = element_text(face = "bold", size=0),          # italicize titles
    axis.title.x = element_text(size = 24),              # optional size tweaks
    axis.title.y = element_text(size = 24)
  ) &
  labs(
    y = "Expression",
    x = "Cluster"
  )


## Convert Seurat object to SCE to generate NE score for Fig. 1l #

# Convert seurat object to SCE
library(SingleCellExperiment)
library(SummarizedExperiment)
sce<-as.SingleCellExperiment(RPM_Allo)

# Import the NE/Non-NE data from the literature (Zhang et al, TLCR, 2018)
nedat_mouse<-read.csv("Signatures/nedat_mouse.csv")
nedat_mouse$gene_name
gene_idx <- rownames(sce) %in% nedat_mouse$gene_name
table(gene_idx)
X <- as.matrix(assay(sce, "logcounts")[gene_idx,])
X
# assure the gene order matches up
X <- X[nedat_mouse$gene_name,]
head(X)
dim(X)

# define the scoring function based on methods from Zhang et al TLCR 2018
ne_score <- function(cell, ne, notne, ...) { 
  # report missing value is no variation 
  if (sd(cell) == 0) { 
    NA
  } else { 
    (cor(x = cell, y = ne, ...) - cor(x = cell, y = notne, ...))/2
  }
}

sce$NE_spearman <- apply(X = X, 2, ne_score,  
                         ne = nedat_mouse$NE, 
                         notne = nedat_mouse$NonNE, 
                         method = "spearman")


# Add NE score to seurat object for Fig. 1l # 
ne_score<-sce$NE_spearman
RPM_Allo@meta.data$NE_spearman<-ne_score
RPM_Allo@meta.data$NE_spearman


# Fig 1l plots #
FeaturePlot(RPM_Allo, features = c("NE_spearman"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
VlnPlot(RPM_Allo,features = c("NE_spearman"), same.y.lims=FALSE,group.by=c("Pheno"),cols=pheno_col,pt.size=.01)


###################### Violin plotting signature scores using Scater #####################
###### For plots in Fig. 1l-n
# For example, NE score by phenotype #
library(scater)

# Plot ATOH1_Targets, for example
x<-scater::plotColData(sce, x = "Pheno", y = "ATOH1_Targets1", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)
x+ theme(axis.title.y = element_text(size = 24),axis.text.y = element_text(size = 16),
         axis.text.x = element_blank(),   # rotate x-axis labels
         plot.title = element_blank(),                      # remove plot title
         legend.position = "none"                           # remove legend
) +
  labs(
    x = "",             # custom x-axis title
    y = "Signature score"  # custom y-axis title
  )+  ylim(-.1,.3)+ geom_boxplot(fill=pheno_col, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)

### Example code to perform Kruskal-Wallis with post-hoc Dunn's on all signature scores, p-vals added manually in illustrator #
# Extract the data
expr_values <- FetchData(RPM_Allo, vars = c("NE_spearman", "Pheno"))

# Rename columns for convenience
colnames(expr_values) <- c("expression", "group")

# Run Kruskal-Wallis test
kruskal_result <- kruskal.test(expression ~ group, data = expr_values)
print(kruskal_result)

# If significant, do Dunn's test
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunnTest(expression ~ group, data = expr_values,method="bonferroni")
  print(dunn_result)
}

# Convert to data frame if needed
dunn_df <- dunn_result$res

# Function to assign p-value stars
get_p_stars <- function(p) {
  if (p < 0.0001) return("****")
  else if(p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# Apply function to create a new column
dunn_df$p_stars <- sapply(dunn_df$P.adj, get_p_stars)

# View result, Use p-val stars to annotate plots
print(dunn_df)

# Repeat above KW with post-Hoc Dunns for each signature set # 

# Save data #
saveRDS(RPM_Allo,"05_2025_RPM_AllograftOnly_Fig2e.rds")
# To read in later
# RPM_Allo<-readRDS("05_2025_RPM_AllograftOnly_Fig2e.rds")

## END of script. Any additional formatting occurred in Adobe Illustrator ##
