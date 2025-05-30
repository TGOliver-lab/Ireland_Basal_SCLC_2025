# RPM K5 Cre vs CGRP Cre Primary Tumor Analysis
# Ireland et al, 2025
# Related to Fig. 1 h-k and Extended Data Fig. 2e-i

setwd("/Users/abbieireland/Desktop/scRNAseq")
# If using Isilon share
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig1_RPMK5_CGRP_Primary")

# Load necessary packages
#load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(CellTagR)
  library(viridis)
  library(ggpubr)
  library(SeuratData)
  library(SeuratDisk)
  library(zellkonverter)
})


###############################################################################################################
############## Seurat conversion and analysis using leiden clustering from scanpy ##########
########################################################################################################

######### Now, bring into Seurat, from  anndata object processed in scanpy/scvi ##########################################################
## Analysis to generate anndata object is located at https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig1_ExtFig2_RPMK5vsCGRP_Final_Clean.ipynb

# Read in adata object as SCE
getwd()
adata<-readH5AD("021825_RPM_CGRPvK5_adata3.h5ad")
counts(adata)

# Check dataset specs
adata
table(adata$GenoCre)
# RPM_CGRP   RPM_K5 
# 7628    20046 
table(adata$leiden_scVI_1.3)
# 0    1    2    3    4    5    6    7    8    9   10 
# 7266 4993 3405 2299 2152 1839 1415 1366 1311 1116  512 

table(adata$UnID)
# RPM_CGRP1    RPM_CGRP2    RPM_CGRP3    RPM_CGRP4    RPM_CGRP5    RPM_CGRP6    RPM_CGRP7 
# 1271          189           63           96          184         2813         3012 
# RPM_K5_total      RPM_K5c      RPM_K5t 
# 7641             2049          10356 

#Convert SCE to seurat
RPM_K5vCGRP <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_K5vCGRP

# An object of class Seurat 
# 55491 features across 27674 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
RPM_K5vCGRP[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_K5vCGRP)<-'norm'
table(RPM_K5vCGRP@meta.data$UnID)

# Add embeddings of umap, and X_scVI_1.3 to assay norm
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_K5vCGRP)
dim(test)
RPM_K5vCGRP[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_K5vCGRP[['umap']]@cell.embeddings

# Plots to generate Fig 1h
unid_cols<-c("darkorchid4","orange")
DimPlot(RPM_K5vCGRP,split.by='Cre',group.by='Cre',cols=unid_cols,reduction='umap',shuffle=TRUE)+NoAxes()
DimPlot(RPM_K5vCGRP,group.by='Cre',cols=unid_cols,reduction='umap',shuffle=TRUE)+NoAxes()

# Plots to generate Fig 1i
# UMAP by Leiden cluster
colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#0aa6d8',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
DimPlot(RPM_K5vCGRP,group.by='leiden_scVI_1.3',cols=colors, reduction='umap',label=TRUE,label.size=6)&NoAxes()

# Barplot by Leiden cluster
### Look at distribution of leiden clusters in tumors by Cre (Fig. 1i)
library(ggplot2)
library(ggpubr)

##### What % of cells occupy what cluster for each sample ####
Idents(RPM_K5vCGRP)<-'Cre'
Idents(RPM_K5vCGRP)

x<-table(Idents(RPM_K5vCGRP),RPM_K5vCGRP@meta.data$leiden_scVI_1.3)
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =colors)+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=colors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


# For visualization of scores, log transform and normalize raw counts in Seurat
DefaultAssay(RPM_K5vCGRP)<-'RNA'
RPM_K5vCGRP<-NormalizeData(RPM_K5vCGRP)
RPM_K5vCGRP

# Visualize A, N, P expression by split violin plot for Fig. 1j
library(viridis)
VlnPlot(RPM_K5vCGRP,features = c("Ascl1","Neurod1","Pou2f3"), same.y.lims=FALSE,group.by=c("Genotype"),split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)

# Example code to generate wilcoxon rank sum stat information for Fig. 1j, adjust y.max higher until p value can be visualized
a<-VlnPlot(RPM_K5vCGRP,features = c("Ascl1"), group.by=c("Cre"),same.y.lims=FALSE,cols=c("darkorchid4","orange"),pt.size=.01,alpha=0.05, ncol=1, y.max=5)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9)) +stat_compare_means(method = "wilcox.test",label = "p.format",hide.ns = TRUE,comparisons = list(c("CGRP", "K5")))


################ Now, apply gene signatures ######################
# Read in this table to convert between mouse and human homologs
mouse94<-read.csv("Signatures/mouse94.csv")

## Archetype signature analysis (Fig. 1k and Extended Data Fig. 2g) ##

arch<-read.csv("Signatures/Archetype_Sigs_Maddox.csv")
a<-arch$SCLC.A
a2<-arch$SCLC.A2
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


RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')

#VlnPlot(RPM_K5vCGRP,features = c("A_Archetype1","A2_Archetype1","N_Archetype1","P_Archetype1"), group.by=c("Genotype"),same.y.lims=FALSE,split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=5)

# Violin plot of A2 archetype by leiden cluster for Fig. 1k
a<-VlnPlot(RPM_K5vCGRP,features = c("A2_Archetype1"), group.by=c("leiden_scVI_1.3"),same.y.lims=FALSE,cols=colors,pt.size=.01,alpha=0.05, ncol=1)
# Add mean point
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "red",
                 position = position_dodge(width = 0.9)) 


# Example code to perform KW test with post-Hoc dunn's on gene signatures by leiden cluster

######### KRUSKAL-WALLIS on VIOLIN PLOTS ############
### Perform kruskal-wallis test and post-hoc dunn's ##
library(dplyr)
library(FSA)       # For dunnTest
library(rstatix)   # Optional, nice formatting of results
library(patchwork)

# Extract the data
expr_values <- FetchData(RPM_K5vCGRP, vars = c("A2_Archetype1", "leiden_scVI_1.3"))

# Rename columns for convenience
colnames(expr_values) <- c("expression", "group")

# Run Kruskal-Wallis test
kruskal_result <- kruskal.test(expression ~ group, data = expr_values)
print(kruskal_result)

# If significant, do Dunn's test (if set method to "none", that would be uncorrected)
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunnTest(expression ~ group, data = expr_values, method = "none")
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


## Example code to generate violin plot by Cre for Ext. Data. Fig. 3g ##
# Adjust y.max per signature for best fit
a<-VlnPlot(RPM_K5vCGRP,features = c("A_Archetype1"), group.by=c("Cre"),same.y.lims=FALSE,cols=c("darkorchid4","orange"),pt.size=.01,alpha=0.05, ncol=1, y.max=.3)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9)) +stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("CGRP", "K5")))


#############################################################################
# Apply human A N P scRNA seq sigs from Chan et al (Extended Data Fig. 2h)
sc_sclc_sigs<-read.csv("Signatures/hSCLC_Chan_sigs.csv")
sc_sclc_sigs

sc_sclc_sigs$hSCLC_A
a_sc<-sc_sclc_sigs$hSCLC_A[1:67]
a_sc<-subset(mouse94, mouse94$human_homolog %in% a_sc)
a_sc_mouse<-(a_sc$gene_name)
a_sc_mouse

sc_sclc_sigs$hSCLC_N

n_sc<-sc_sclc_sigs$hSCLC_N[1:73]
n_sc<-subset(mouse94, mouse94$human_homolog %in% n_sc)
n_sc_mouse<-(n_sc$gene_name)
n_sc_mouse

sc_sclc_sigs$hSCLC_P

p_sc<-sc_sclc_sigs$hSCLC_P
p_sc<-subset(mouse94, mouse94$human_homolog %in% p_sc)
p_sc_mouse<-(p_sc$gene_name)
p_sc_mouse


RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a_sc_mouse),
  name = 'hSCLC_A')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(n_sc_mouse),
  name = 'hSCLC_N')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(p_sc_mouse),
  name = 'hSCLC_P')

library(viridis)
library(ggplot2)
library(gridExtra)

## Example code to generate violin plot by Cre for Ext. Data. Fig. 3g ##
# Adjust y.max per signature for best fit
a<-VlnPlot(RPM_K5vCGRP,features = c("hSCLC_A1"), group.by=c("Cre"),same.y.lims=FALSE,cols=c("darkorchid4","orange"),pt.size=.01,alpha=0.05, ncol=1, y.max=0.6)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9))+stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("CGRP", "K5")))

###############################################
## Apply A, N, P ChIP target genes as scores (Extended Data Fig. 1e)

chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_m<-(pchip$gene_name)
pchip_m

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(achip),
  name = 'ASCL1_Targets')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(pchip_m),
  name = 'POU2F3_Targets')

## Example code to generate violin plot by Cre for Ext. Data. Fig. 3e ##
# Adjust y.max per signature for best fit
a<-VlnPlot(RPM_K5vCGRP,features = c("ASCL1_Targets1"), group.by=c("Cre"),same.y.lims=FALSE,cols=c("darkorchid4","orange"),pt.size=.01,alpha=0.05, ncol=1, y.max=0.6)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9))+stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("CGRP", "K5")))

### Look at basal vs luminal hillock cell state sigs in CGRP vs K5 ####

hillock<-read.csv("Signatures/luminal_basal_hillock.csv")
basal<-hillock$Basal_hillock[1:170]
basal
lum<-hillock$Luminal_hillock
lum


RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(basal),
  name = 'Basal_hillock')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(lum),
  name = 'Luminal_hillock')

# Example code to generate Basal and Luminal Hillock signature violin plots for Fig. 1k with stats
a<-VlnPlot(RPM_K5vCGRP,features = c("Luminal_hillock1"), group.by=c("leiden_scVI_1.3"),same.y.lims=FALSE,cols=colors,pt.size=.05,alpha=0.05, ncol=1)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "red",
                 position = position_dodge(width = 0.9)) 

### Example code to perform kruskal-wallis test and post-hoc dunn's on violin plots in Fig. 1k##
# Extract the data
expr_values <- FetchData(RPM_K5vCGRP, vars = c("Luminal_hillock1", "leiden_scVI_1.3"))

# Rename columns for convenience
colnames(expr_values) <- c("expression", "group")

# Run Kruskal-Wallis test
kruskal_result <- kruskal.test(expression ~ group, data = expr_values)
print(kruskal_result)

# If significant, do Dunn's test (if set method to "none", that would be uncorrected)
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunnTest(expression ~ group, data = expr_values, method = "none")
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


######## Add NE Score for Ext. Data Fig. 2f ##########
### For NE score assignment, originally describe in Zhang, TLCR, 2017 ##

#### Converting to SCE to Add NE Score ####
library(SingleCellExperiment)
library(SummarizedExperiment)

sce<-as.SingleCellExperiment(RPM_K5vCGRP)

# import the NE/Non-NE data from the literature
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

# get the NE scores. note that 4 cells have missing values b/c 
# they had zero variation across the 50 genes. 

sce$NE_spearman <- apply(X = X, 2, ne_score,  
                         ne = nedat_mouse$NE, 
                         notne = nedat_mouse$NonNE, 
                         method = "spearman")


#Add NE score as metadata to Seurat object
ne_score<-sce$NE_spearman
RPM_K5vCGRP@meta.data$NE_spearman<-ne_score

# Generate violin plot with wilcoxon test for Ext. Data Fig. 2f
a<-VlnPlot(RPM_K5vCGRP,features = c("NE_spearman"), group.by=c("Cre"),same.y.lims=FALSE,cols=c("darkorchid4","orange"),pt.size=.05,alpha=0.3, ncol=1, y.max=1)
# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9)) +stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("CGRP", "K5")))



# Save resulting seurat object with signatures
saveRDS(RPM_K5vCGRP,"05_2025_RPMK5vCGRP_Norm_Allsigs_New.rds")

# To read in later
# RPM_K5vCGRP<-readRDS("05_2025_RPMK5vCGRP_Norm_Allsigs_New.rds")

## End of script, any additional formatting changes made in Adobe Illustrator ##
