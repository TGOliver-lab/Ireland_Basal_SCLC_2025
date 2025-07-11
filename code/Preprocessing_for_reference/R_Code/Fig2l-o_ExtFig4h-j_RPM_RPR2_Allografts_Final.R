# RPM vs RPR2 basal-organoid-derived allograft tumour analysis 
# Ireland et al, 2025
# Related to *Old Fig. 2l-o and Extended Data Fig. 4h-j*
# Related to *Final Ext. Data Fig. 4i-n*
# Load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(SeuratData)
  library(SeuratDisk)
  library(CellTagR)
  library(viridis)
  library(ggpubr)
  library(zellkonverter)
})

setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat/scanpy/RPM_RPR2")

# If using Isilon share
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig2_ExtFig4_RPM_RPR2_Allos")

#### Pull in existing anndata object from Scanpy into Seurat ####
## Anndata object previously generated as in: https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig2l-o_ExtFig4h-i_RPMvRPR2_Basal_Allografts_Final_Clean.ipynb

# Read in adata object as SCE
adata<-readH5AD("050925_RPM_TBOAllo_OriginalandAllo3_RPR2_adata2.h5ad")

# Check adata metadata to ensure correct dataset
table(adata$Genotype)
# RPM  RPR2 
# 10572  8792 
table(adata$leiden_scVI_1.2)
# 0    1    2    3    4    5    6    7    8    9 
# 3608 3049 2686 2558 2297 1858 1498 1406  287  117 


#Convert SCE to seurat
RPMvRPR2 <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPMvRPR2
# An object of class Seurat 
# 55491 features across 19364 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

DefaultAssay(RPMvRPR2)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))

RPMvRPR2[["norm"]] <- norm_assay


# Set assay norm as default assay for seurat object
DefaultAssay(RPMvRPR2)<-'norm'

# Add embeddings of umap, and X_scVI_1.2 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPMvRPR2)
dim(test)
RPMvRPR2[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPMvRPR2[['umap']]@cell.embeddings

# Define sample color vector and use to plot UMAP in Ext. Data Fig. 4i #
unid_cols<-c("turquoise3","maroon")
DimPlot(RPMvRPR2,group.by='Genotype',cols=unid_cols,reduction='umap',shuffle=TRUE,label=FALSE)+NoAxes()

# Define color vector and use to plot UMAP in Fig. 2m #
colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#4ba4d3',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
DimPlot(RPMvRPR2,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=TRUE,label.size=10)&NoAxes()&NoLegend()&ggtitle("Fig 2m clusters")



## Look at where RPM phenotypes are in this space for Ext. Data Fig. 4l ##
# If on Isilon share
#RPM_Allo<readRDS("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig2_ExtFig4_RPM_Organoids_Allo/05_2025_RPM_AllograftOnly_Fig2e.rds")

RPM_Allo<-readRDS("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat/scanpy/050625_RPM_Allo_AllSigs_Seurat.rds")

Idents(RPM_Allo)<-'Pheno'
table(Idents(RPM_Allo))
# NE NE_Neuronal    Neuronal       Atoh1        Tuft       Basal 
# 2173         864         530         575          11         282 

# Assign these as new metadata column 
RPMvRPR2@meta.data$RPM_Pheno<- ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "Basal", ]), "Basal",
                                      ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "NE", ]), "NE",
                                      ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "NE_Neuronal", ]), "NE_Neuronal",
                                      ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "Neuronal", ]), "Neuronal",
                                      ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "Atoh1", ]), "Atoh1",
                                      ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPM_Allo@meta.data[RPM_Allo@meta.data$Pheno == "Tuft", ]), "Tuft",
                                             ifelse(rownames(RPMvRPR2@meta.data) %in% rownames(RPMvRPR2@meta.data[RPMvRPR2@meta.data$Genotype == "RPR2", ]), "RPR2", "Undetermined")))))))
                                      
                                                 


table(RPMvRPR2@meta.data$RPM_Pheno)
# Atoh1        Basal           NE  NE_Neuronal     Neuronal         RPR2         Tuft 
# 575          281         2158          863          524         8792           11 
# Undetermined 
# 6160 

pheno_col<-c("brown2","darkorchid4","dodgerblue","#66A61E","orange","turquoise4","maroon","gray90")

RPMvRPR2@meta.data$RPM_Pheno<-factor(RPMvRPR2@meta.data$RPM_Pheno, c("NE","NE_Neuronal","Neuronal","Atoh1","Tuft","Basal","RPR2","Undetermined"))

# Ext. Data Fig. 4l UMAP #
DimPlot(RPMvRPR2,group.by='RPM_Pheno',cols=pheno_col,reduction='umap',shuffle=TRUE)+NoAxes()+ggtitle("RPM Cell State from Fig 2f")


######### For gene signature assignment and violin plots, apply log normalization to counts assay ###
DefaultAssay(RPMvRPR2)<-'RNA'
RPMvRPR2<-NormalizeData(RPMvRPR2)

# Ext. Data Fig. 4i #
## Example code to generate violin plots comparing RPM to RPR2 allos for Ext Data Fig. 4h + Wilcoxon rank-sum test ##
# Repeat for each gene of interest, adjusting y.max to visualize range and p-value per gene #
a<-VlnPlot(RPMvRPR2,features = c("Mycl"), group.by=c("Genotype"),same.y.lims=FALSE,cols=c("turquoise3","maroon"),
           pt.size=.01,alpha=0.5, ncol=1, y.max=4)

# Add mean point and do wilcoxon rank sum test 
a +stat_summary(fun = mean,
                geom = "point",
                shape = 18,
                size = 2,
                color = "yellow",
                position = position_dodge(width = 0.9)) +
  stat_compare_means(method = "wilcox.test", size=8,
                     label = "p.signif",
                     hide.ns = TRUE,
                     comparisons = list(c("RPM", "RPR2"))) +
  labs(title = NULL, x = NULL, y = "Expression") +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

# Apply various gene signatures #

# Read in this table to convert human gene signatures to mouse homologs #
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")

######################################################
## Look at SCLC-archetype signatures as scores (Ext Data Fig 4m)

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


RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')

# Extended Data Fig. 4i #
VlnPlot(RPMvRPR2, features = c("A_Archetype1","A2_Archetype1","N_Archetype1","P_Archetype1"),group.by=c("Genotype"), cols=c("turquoise3","maroon"),alpha=0.05,ncol=4)

#################################################################
## Look at A, N, P ChIP targets as scores (Ext Data Fig. 4n)
chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ATOH_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_m<-(pchip$gene_name)
pchip_m

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(achip),
  name = 'ASCL1_Targets')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(pchip_m),
  name = 'POU2F3_Targets')

# Also add MYC ChIP data
MYC_targets<-read.csv("Signatures/50_Conserved_MYC_Targets.csv")
MYC_targets<-MYC_targets$Gene
MYC_targets

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(MYC_targets),
  name = 'MYC_Targets')

## Ext Data Fig 4n ##
VlnPlot(RPMvRPR2, features = c("ASCL1_Targets1","NEUROD1_Targets1","POU2F3_Targets1","MYC_Targets1"),group.by=c("Genotype"), cols=c("turquoise3","maroon"),alpha=.06, ncol=4)


# Example code to visualize mean and also determine p-val by Wilcoxon rank sum for Ext. Fig 4m-n #
# Repeat for each signature, adjusting y.max accordingly

a<-VlnPlot(RPMvRPR2,features = c("A_Archetype1"), group.by=c("Genotype"),same.y.lims=FALSE,cols=c("turquoise3","maroon"),
           pt.size=.01,alpha=0.5, ncol=1, y.max=0.3)

# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "yellow",
                 position = position_dodge(width = 0.9)) +
  stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("RPM", "RPR2")))



######## Next, add NE Score for Ext. Data Fig. 4k ##########

#### Convert to SCE to add NE Score ###
# Convert seurat object to SCE
library(SingleCellExperiment)
library(SummarizedExperiment)
sce<-as.SingleCellExperiment(RPMvRPR2)


# Import the NE/Non-NE data from the literature
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


# Add NE score to seurat object
ne_score<-sce$NE_spearman
RPMvRPR2@meta.data$NE_spearman<-ne_score
RPMvRPR2@meta.data$NE_spearman

### UMAP for Ext. Data Fig. 4k ###
FeaturePlot(RPMvRPR2, features = c("NE_spearman"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

### Violin Plot for Ext. Data Fig. 4k ###
# Assign pheno/cell states 
RPMvRPR2@meta.data$Geno_C8<- ifelse(RPMvRPR2@meta.data$leiden_scVI_1.2 %in% "8", "Cluster_8",
                                    ifelse(RPMvRPR2@meta.data$Genotype %in% "RPM", "RPM",
                                    ifelse(RPMvRPR2@meta.data$Genotype %in% "RPR2", "RPR2","NA")))

RPMvRPR2@meta.data$Geno_C8<-factor(RPMvRPR2@meta.data$Geno_C8, c("RPM","RPR2","Cluster_8"))

table(RPMvRPR2@meta.data$Geno_C8)
# RPM      RPR2 Cluster_8 
# 10286      8791       287 

# Ext. Data Fig. 4k Violin plot final# 
x<-VlnPlot(RPMvRPR2,features = c("NE_spearman"), same.y.lims=FALSE,group.by=c("Geno_C8"),cols=c("turquoise3","maroon","#A0522D"),alpha=.2, ncol=1)
x
fortest<-x+ geom_boxplot(fill=c("turquoise3","maroon","#A0522D"),alpha=1, position = position_dodge(width = .2),size=0.5,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)+ggtitle("NE score")
fortest

### Perform kruskal-wallis test and post-hoc dunn's ##
library(dplyr)
library(FSA)
library(rstatix)   # Optional, nice formatting of results

# Extract the data
expr_values <- FetchData(RPMvRPR2, vars = c("NE_spearman", "Geno_C8"))
head(expr_values)
# Rename columns for convenience
colnames(expr_values) <- c("expression", "group")

# Run Kruskal-Wallis test
kruskal_result <- kruskal.test(expression ~ group, data = expr_values)
print(kruskal_result)

# If significant, do Dunn's test (if set method to "none", that would be uncorrected)
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunnTest(expression ~ group, data = expr_values)
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

# View result
print(dunn_df)


saveRDS(RPMvRPR2,"05_2025_RPMvRPR2_Allografts_Final.rds")

# To read in later
#RPMvRPR2<-readRDS("05_2025_RPMvRPR2_Allografts_Final.rds")




