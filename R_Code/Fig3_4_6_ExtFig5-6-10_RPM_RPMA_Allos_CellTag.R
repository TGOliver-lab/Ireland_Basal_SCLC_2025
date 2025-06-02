# RPM (CellTagged pre-Cre and CellTagged post-Cre) and RPM basal-organoid-derived allograft tumour analysis & CellTag analyses
# Ireland et al, 2025
# Related to Fig. 3e-l, Fig. 4c-f,h, Fig. 6d,e,g,h
# Also related to Extended Data Fig. 5e,f, Extended Data Fig. 6a-e, and Extended Data Fig. 10e


############################## CELLTAG EXPERIMENT KEY ##########################
# Matches metadata on GEO deposit GSE279200

# CellTag Experiment "A"= RPM matching organoids and resulting allograft, CellTagged pre-Cre

# CellTag Experiment "B"= RPM matching organoids and resulting allografts (n=2), CellTagged post-Cre; 

# CellTag Experiment "C"= RPMA matching organoids and resulting allograft, CellTagged pre-Cre.
############################################################################################################## First, analyze CellTag data from RPM organoids and allografts CellTaged post-Cre (Experiment B)  "#############################################


######################################################################################################################################################################################################
########## Moving into Seurat with four distinct CellTagged allograft tumour samples (RPM CellTagged pre-Cre Allo, RPM CellTagged post-Cre (Allo3, Allo4), and RPMA Allo (CellTagged pre-Cre)) ###############
######################## Bringing anndata object in from Scanpy/scVI analysis ######################
# Anndata object previously generated as in: https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig3-4_ExtFig5-6_RPM_RPMA_Allografts_CellTag_CellRank_Final_Clean.ipynb

# Load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(viridis)
  library(ggpubr)
  library(Seurat)
  library(SeuratData)
  library(SeuratDisk)
  library(zellkonverter)
})


setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat")
## If using Isilon share...
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig3-4-6_ExtFig5-6-10_RPM_RPMA_Allos_UMAP_FA_CellTag/RPM_RPMA_Seurat_CellTag")

# Read in adata object as SCE
adata<-readH5AD("042725_RPM_RPMA_TBOAllo_CellTagAnalysis_New_1.2.h5ad")

# Check object metadata to ensure correct dataset
table(adata$Cre)
# CT_post-Cre  CT_pre-Cre 
# 13578       13040 

table(adata$Genotype)
# RPM  RPMA 
# 16458 10160 

table(adata$GenoCT)
# RPMA_CTpreCre RPM_CTpostCre  RPM_CTpreCre 
# 10160         13578          2880 

table(adata$leiden_scVI_1.2)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 3020 2872 2726 2468 2284 1783 1343 1325 1312 1050  990  961  939  906  849  672  422  277  193  171   55 

table(adata$UnID)
# "RPM_Allo_New" in this metadata is the "CellTag pre-Cre" RPM allograft

# RPMA_Allo    RPM_Allo3    RPM_Allo4  RPM_Allo_New 
# 10160         6718         6860         2880 

dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)

#Convert SCE to seurat
TBO_seurat <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
TBO_seurat
# An object of class Seurat 
# 55491 features across 26618 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

TBO_seurat
TBO_seurat[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(TBO_seurat)<-'norm'

# Add embeddings of umap, and X_scVI_1.2 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(TBO_seurat)
dim(test)
TBO_seurat[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
TBO_seurat[['umap']]@cell.embeddings


# Define color vector to plot UMAP by leiden cluster for Fig. 3f #
# Fig. 3f #

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
  "#66C2A5", # teal
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
  "#66A61E"  # olive green (NOT same green as before)
)

# Fig. 3e UMAP #
DimPlot(TBO_seurat,group.by='Genotype',cols=c("darkorchid4","orange"), reduction='umap',label=FALSE,label.size=7)&NoAxes()

# Plot by individual samples if desired
DimPlot(TBO_seurat,group.by='UnID',cols=my_colors, reduction='umap',shuffle=TRUE,label=FALSE,label.size=6) & NoAxes()

# Fig. 3f UMAP by Leiden #
DimPlot(TBO_seurat,group.by='leiden_scVI_1.2',cols=my_colors, reduction='umap',label=TRUE,label.size=7)&NoAxes()


### Look at distribution leiden clusters per sample for Fig. 3f barplot ##
Idents(TBO_seurat)<-'leiden_scVI_1.2'

x<-table(TBO_seurat@meta.data$Genotype,Idents(TBO_seurat))
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=my_colors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


# Assign phenotype for all tumors cells for Fig. 3h (based on gene expression and signature expression patterns generated below) #
#Basal- 17
#NE/A- 1, 7, 11, 14, 18, 19, 20
#N/NE- 0, 8
#N- 4, 9, 15
#Atoh-6
#P - 16
#Triple-Neg/Mixed/Subtype-Low (SL) - 13, 2, 3, 5, 10, 12

TBO_seurat@meta.data$Pheno<- ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("16"), "Tuft",
                                    ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("1","7","11","14","18","19","20"), "NE",
                                           ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("4","9","15"), "Neuronal",
                                                  ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("0","8"), "NE/Neuronal",
                                                         ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("17"), "Basal",
                                                                ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("6"), "ATOH1",
                                                                       ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("13","2","3","5","10","12"), "Triple-Neg","NA")))))))




table(TBO_seurat@meta.data$Pheno)
TBO_seurat@meta.data$Pheno<-factor(TBO_seurat@meta.data$Pheno, levels=c("NE","NE/Neuronal","Neuronal","ATOH1","Tuft","Triple-Neg","Basal"))

# Define fate color vector and use to plot by cell state for Fig. 3h #
pheno_col<-c("brown2","darkorchid4","dodgerblue","#66A61E","orange","turquoise4","turquoise")
DimPlot(TBO_seurat,group.by='Pheno',cols=pheno_col, reduction='umap',label=FALSE,label.size=6) & NoAxes()

##### What % of cells occupy what cell fate per genotype? For Fig. 3h barplot. ###
Idents(TBO_seurat)<-'Pheno'

x<-table(TBO_seurat@meta.data$Genotype,Idents(TBO_seurat))
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")


ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =pheno_col)+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=pheno_col)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


########################################################################
# Before signature assessment, normalize and log count data 
########################################################################
DefaultAssay(TBO_seurat)<-'RNA'
TBO_seurat<-NormalizeData(TBO_seurat)

## Perform Cell Cycle Analysis #
########################################################################
# Read in mouse94 to convert human genes to mouse homologs #
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

TBO_seurat <- CellCycleScoring(TBO_seurat, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
TBO_seurat@meta.data$Phase<-factor(TBO_seurat@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(TBO_seurat, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()

## Visualize individual genes in violin plots for Extended Data Fig. 5 ##

#### Calculate fractions of + vs - for Ext Data Fig. 5c ###
## Example code for fraction + for Neurod1 ##
## Replace Neurod1 with any gene of interest and repeat ##

# Extract expression and metadata
gene <- "Neurod1"
data <- FetchData(TBO_seurat, vars = c(gene, "Genotype"))
head(data)

# Create a binary column: 1 for Positive (Expression > 0.01), 0 for None (Expression < 0.01)
data <- data %>%
  mutate(Expression_Status = ifelse(data$Neurod1 > 0.01, "Positive", "Negative"))


## Then, plot

x<-table(data$Genotype,data$Expression_Status)
x
proportions <- as.data.frame(prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Fraction")


ggbarplot(proportions, x="Sample", y="Fraction", fill = "Sample", group = "Sample", ylab = "Fraction", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)
proportions
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Fraction, x=Cluster)) +
  geom_bar(position="stack", stat="identity")
p <- ggplot(proportions,
            aes(fill = Sample, y = Fraction, x = Cluster)) +
  geom_bar(position = "stack", stat = "identity")
p + scale_fill_manual(values=c("blue4","red3"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=0), 
                                                                           plot.title = element_text(size = 18, face = "plain", hjust = 0.5))+labs(y = "Fraction of Neurod1-hi cells (>0.01)")



# Extended Data Fig. 5e #
a<-VlnPlot(TBO_seurat, features = c("Ascl1"),group.by="leiden_scVI_1.2", cols=my_colors,pt.size=.01,alpha=0.3, ncol=1)+ theme(axis.text.x = element_text(size = 8), legend.position = "none", axis.title.x = element_text(size = 20),
                                                                                                                              axis.title.y = element_text(size =20))+ggtitle("")+labs(x = "Leiden Cluster", y = "Expression")
b<-VlnPlot(TBO_seurat, features = c("Neurod1"),group.by="leiden_scVI_1.2", cols=my_colors,pt.size=.01,alpha=0.3, ncol=1)+ theme(axis.text.x = element_text(size = 8), legend.position = "none", axis.title.x = element_text(size = 20),
                                                                                                                                axis.title.y = element_text(size =20))+ggtitle("")+labs(x = "Leiden Cluster", y = "Expression")
c<-VlnPlot(TBO_seurat, features = c("Pou2f3"),group.by="leiden_scVI_1.2", cols=my_colors,pt.size=.01,alpha=0.3, ncol=1)+ theme(axis.text.x = element_text(size = 8), legend.position = "none", axis.title.x = element_text(size = 20),
                                                                                                                               axis.title.y = element_text(size =20))+ggtitle("")+labs(x = "Leiden Cluster", y = "Expression")
library(cowplot)
plot_grid(a,b,c, ncol = 3)


####### Kruskal Wallis and post-hoc Dunns on gene expression Ext. Data Fig. 5e ######
library(ggpubr)

# Define the genes
genes <- c("Ascl1","Neurod1","Pou2f3")

# Create violin plots with Kruskal-Wallis test
plots <- lapply(genes, function(gene) {
  p <- VlnPlot(TBO_seurat,
               features = gene,
               group.by = "leiden_scVI_1.2",
               cols = my_colors,
               alpha = 0.5) +  # smaller dots
    stat_compare_means(method = "kruskal.test", 
                       label = "p.format", 
                       label.x = 2,size = 5) +
    ggtitle("") +
    theme(plot.title = element_text(face = "italic"), legend.position = "none",  # Remove legend
          axis.title.x = element_blank(),axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 14)) +  # Remove x-axis label
    labs(y = "Expression")  # Change y-axis label to "Expression"
  return(p)
})

# Arrange plots
patchwork::wrap_plots(plots, ncol = 3)


########################################################################
######### Add Signature Data for Fig. 3,6 and Ext. Data Fig. 10 ##################
########################################################################
library(viridis)
library(ggplot2)
library(ggpubr)

## Look at A, N, P, ATOH1 ChIP targets and YAP activity as scores (Fig. 3j) ##
chip<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/ASCL1_NEUROD1_POU2F3_ATOH_ChIP_Targets.csv")
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


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(achip),
  name = 'ASCL1_Targets')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nchip),
  name = 'NEUROD1_Targets')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(pchip_mouse),
  name = 'POU2F3_Targets')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(atoh_mouse),
  name = 'ATOH1_Targets')

# Also add YAP1/TAZ activity score #
yap<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Yap-Taz_signatures.csv")

# This is from Wang 
yap<-yap$mYapTaz_Wang_genes[1:22]
yap

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(yap),
  name = 'YAP_Activity')


# UMAPs in Fig. 3j #
FeaturePlot(TBO_seurat, features = c("ASCL1_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("NEUROD1_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("ATOH1_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("POU2F3_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("YAP_Activity1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()


##### Now add VlnPlots by genotype with wilcoxon rank sum tests for Fig. 3j ###
## Example code for ATOH1 Targets, replace with target gene signature of interest and repeat ##
a <- VlnPlot(
  TBO_seurat,
  features = c("ATOH1_Targets1"),
  group.by = "Genotype",
  cols = c("darkorchid4", "orange"),
  pt.size = 0.01,
  alpha = 0.05,
  ncol = 1) + theme(axis.title.y = element_text(size = 24),
    axis.text.x = element_text(angle = 0, hjust = 0.5,size=24),  # rotate x-axis labels
    plot.title = element_blank(),                      # remove plot title
    legend.position = "none"                           # remove legend
  ) +
  labs(
    x = "",             # custom x-axis title
    y = "Expression"  # custom y-axis title
  )+  ylim(-.1, .32)

# Add mean point and do wilcoxon rank sum test 
a + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "red",
                 position = position_dodge(width = 0.9)) +stat_compare_means(method = "wilcox.test",label = "p.signif",size=8,hide.ns = TRUE,comparisons = list(c("RPM", "RPMA")))



######## Normal cell Type signatures for Fig. 3k and Ext. Data. Fig. 10e ###########
ct<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Montoro_Consensus.csv")
basal_noC<-ct$Basal_noC[1:41]
basal_noC
ne_noC<-ct$Neuroendocrine_noC[1:32]
ne_noC
tuft_noC<-ct$Tuft_noC[1:60]
tuft_noC

ct<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Montoro_Iono_Human.csv")
ct$Montoro_Iono_Human
ionoh<-subset(mouse94, mouse94$human_homolog %in% ct$Montoro_Iono_Human)
ionoh_hum<-(ionoh$gene_name)

ct<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Montoro_Iono_Extended.csv")
iono_mouse<-ct$Montoro_iono_Ext
iono_mouse


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(tuft_noC),
  name = 'Tuft_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(basal_noC),
  name = 'Basal_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ne_noC),
  name = 'NE_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ionoh_hum),
  name = 'Iono_Human')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(iono_mouse),
  name = 'Iono_Mouse_Ext')


# UMAPs for Fig. 3k #
FeaturePlot(TBO_seurat, features = c("NE_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )
FeaturePlot(TBO_seurat, features = c("Basal_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )
FeaturePlot(TBO_seurat, features = c("Tuft_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )


# UMAPs for Extended Data Fig. 10e #
FeaturePlot(TBO_seurat, features = c("Iono_Mouse_Ext1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )

FeaturePlot(TBO_seurat, features = c("Iono_Human1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )

#########################################################
## Look at SCLC Archetype signatures for Fig. 3l ##
arch<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Archetype_Sigs_Maddox.csv")
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


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')


######## Add NE Score (Fig. 3i) ##########

#### Converting to SCE to Add NE Score and create Violin plots ###

library(SingleCellExperiment)
library(SummarizedExperiment)
sce<-as.SingleCellExperiment(TBO_seurat)


# import the NE/Non-NE data from the literature
nedat_mouse<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/nedat_mouse.csv")
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


###################### Violin plotting signatures of interest #####################
library(scater)

# Plot NE spearman, for example, by cell fate (as in Fig. 3i) #
# Adjust ylim per signature to capture range # 

x<-scater::plotColData(sce, x = "Pheno", y = "NE_spearman", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)
x+ theme(axis.title.y = element_text(size = 24),axis.text.y = element_text(size = 16),
         axis.text.x = element_blank(),   # rotate x-axis labels
         plot.title = element_blank(),                      # remove plot title
         legend.position = "none"                           # remove legend
) +
  labs(
    x = "",             # custom x-axis title
    y = "Signature score"  # custom y-axis title
  )+  ylim(-1,1)+ geom_boxplot(fill=pheno_col, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)



# Perform KRUSKAL-WALLIS and post-Hoc Dunns on same data # 
df <- sce@colData[, c("NE_spearman", "Pheno")]
colnames(df) <- c("score", "group")
df <- na.omit(df)
df

kruskal.test(score ~ group, data = df)

library(FSA)
dunn_results <- dunnTest(score ~ group, data = df, method = "bonferroni")
# Extract results
dunn_df <- dunn_results$res
dunn_df
# Create stat.test data frame for ggpubr::stat_pvalue_manual
# Use separate() to split the Comparison column
dunn_df <- dunn_df %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-")

# Create stat.test data frame for ggpubr::stat_pvalue_manual
stat.test <- data.frame(
  group1 = dunn_df$group1,
  group2 = dunn_df$group2,
  p.adj = dunn_df$P.adj,
  y.position = seq(0.6, 1, length.out = nrow(dunn_df)),  # adjust if needed
  label = ifelse(dunn_df$P.adj < 0.0001, "****",
                 ifelse(dunn_df$P.adj < 0.001, "***",
                        ifelse(dunn_df$P.adj < 0.01, "**",
                               ifelse(dunn_df$P.adj < 0.05, "*", "ns"))))
)

print(stat.test)

# Repeat the above for any signatures of interest including for normal cell type signatures (Fig. 3k, Ext. Data. Fig. 10e) and SCLC-archetype signatures (Fig. 3l)
# Manually add p-values to plot using Illustrator



# Plot NE score by genotype for Fig. 3i
x<-scater::plotColData(sce, x = "Genotype", y = "NE_spearman", colour_by = "Genotype")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=c("darkorchid4","orange"))
x<-x+ theme(axis.title.y = element_text(size = 24),axis.text.y = element_text(size = 16), axis.text.x = element_blank(), plot.title = element_blank(),
            legend.position = "none") + labs(x = "",  y = "Signature score")+  ylim(-1,1)+ geom_boxplot(fill=c("darkorchid4","orange"), alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)

x

## Perform wilcoxon rank-sum test for RPM vs RPMA on the NE score data (Fig. 3i) ##
x + stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2,
                 color = "red",
                 position = position_dodge(width = 0.9)) +stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = TRUE,comparisons = list(c("RPM", "RPMA")))

# Manually add p-value using Illustrator #


## Move back to Seurat object to add NE_spearman score (Fig. 3i) ##
TBO_seurat@meta.data$NE_spearman<-sce$NE_spearman

## Fig. 3i UMAP ##
FeaturePlot(TBO_seurat, features = c("NE_spearman"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()+
  theme(
    legend.position = "none",    # remove legend
    plot.title = element_blank()  # remove title
  )


# Save Seurat object as is before moving into adding human signatures for Fig. 6 #
saveRDS(TBO_seurat,"05_2025_RPM_RPMA_Allos_Seurat.rds")

# To read in later #
setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat")
TBO_seurat<-readRDS("05_2025_RPM_RPMA_Allos_Seurat.rds")


####################################################################################
# Add human signatures for Fig. 6 #
####################################################################################

########## Adding  CARIS signatures (for Fig. 6d) ############

caris<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Caris_Top100_HumanSCLC.csv")
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")
a<-caris$Caris_SCLC.A[1:100]
a
n<-caris$Caris_SCLC.N[1:100]
n
p<-caris$Caris_SCLC.P[1:100]
p
y<-caris$Caris_SCLC.Y[1:100]
y
mixed<-caris$Caris_SCLC.Mixed[1:100]
tn<-caris$Caris_SCLC.TN[1:100]

# Convert to mouse
a<-subset(mouse94, mouse94$human_homolog %in% a)
a<-a$gene_name
n<-subset(mouse94, mouse94$human_homolog %in% n)
n<-n$gene_name
p<-subset(mouse94, mouse94$human_homolog %in% p)
p<-p$gene_name
y<-subset(mouse94, mouse94$human_homolog %in% y)
y<-y$gene_name
mixed<-subset(mouse94, mouse94$human_homolog %in% mixed)
mixed<-mixed$gene_name
tn<-subset(mouse94, mouse94$human_homolog %in% tn)
tn<-tn$gene_name


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a),
  name = 'Caris_A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n),
  name = 'Caris_N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p),
  name = 'Caris_P')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y),
  name = 'Caris_Y')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mixed),
  name = 'Caris_Mixed')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(tn),
  name = 'Caris_TN-')

# Use FeaturePlot to plot UMAPs of each human Caris signature (for Fig. 6d)
# TN- is Lin-
FeaturePlot(TBO_seurat, features = c("Caris_TN-1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()


# Apply NMF-3 signature (SCLC-inflammatory) from Liu et al to our data (Fig 6g) #
liu<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Liu_NMF_Sigs.csv")
nmf3<-liu$NMF3

nmf3<-subset(mouse94, mouse94$human_homolog %in% nmf3)
nmf3_m<-(nmf3$gene_name)
nmf3_m

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nmf3_m),
  name = 'NMF3-I')

# Use FeaturePlot to plot UMAP of Liu NMF3 signature (for Fig. 6g)
FeaturePlot(TBO_seurat, features = c("NMF3-I1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

####### Add Antigen presentation signature (used in Gay et al) (Fig 6g) #######
inf<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Human_inflamed_Gay.csv")
mhc<-inf$MHC

mhc<-subset(mouse94, mouse94$human_homolog %in% mhc)
mhc_m<-mhc$gene_name
mhc_m


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mhc_m),
  name = 'MHC_Sig_Gay')


# Fig 6g feature plot (MHC_Sig_Gay1="Antigen presentation") #
FeaturePlot(TBO_seurat, features = c("MHC_Sig_Gay1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

### Gay et al SCLC-I signature for Fig. 6g ###
inf<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Gay_SCLC-I.csv")
t<-inf$SCLC_I

t<-subset(mouse94, mouse94$human_homolog %in% t)
t_lm<-t$gene_name
t_lm

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(t_lm),
  name = 'Gay_SCLC-I')

FeaturePlot(TBO_seurat, features = c("Gay_SCLC-I1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

#### Converting to SCE to create Violin plots using scater package (for Fig. 6g)
# Convert seurat object to SCE
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)
sce<-as.SingleCellExperiment(TBO_seurat)
TBO_seurat$`Gay_SCLC-I1`
###################### Violin plotting human signatures of interest (for Fig. 3g) #####################
# Plot Gay_SCLC-I for example, replace with any signature of interest, adjust ylim to fit range of each signature 

pheno_col<-c("brown2","darkorchid4","dodgerblue","#66A61E","orange","turquoise4","turquoise")

x<-scater::plotColData(sce, x = "Pheno", y = "Gay_SCLC.I1", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)
x+ theme(axis.title.y = element_text(size = 24),axis.text.y = element_text(size = 16),
         axis.text.x = element_blank(),   # rotate x-axis labels
         plot.title = element_blank(),                      # remove plot title
         legend.position = "none"                           # remove legend
) +
  labs(
    x = "",             # custom x-axis title
    y = "Signature score"  # custom y-axis title
  )+  ylim(-.05,.3)+ geom_boxplot(fill=pheno_col, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)



# Perform Kruskal-Wallis (KW) and post-Hoc Dunn's # 
df <- sce@colData[, c("Gay_SCLC.I1", "Pheno")]
colnames(df) <- c("score", "group")
df <- na.omit(df)
df

kruskal.test(score ~ group, data = df)

library(FSA)
dunn_results <- dunnTest(score ~ group, data = df, method = "bonferroni")
# Extract results
dunn_df <- dunn_results$res
dunn_df
# Create stat.test data frame for ggpubr::stat_pvalue_manual
# Use separate() to split the Comparison column
dunn_df <- dunn_df %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-")

# Create stat.test data frame for ggpubr::stat_pvalue_manual
stat.test <- data.frame(
  group1 = dunn_df$group1,
  group2 = dunn_df$group2,
  p.adj = dunn_df$P.adj,
  y.position = seq(0.6, 1, length.out = nrow(dunn_df)),  # adjust if needed
  label = ifelse(dunn_df$P.adj < 0.0001, "****",
                 ifelse(dunn_df$P.adj < 0.001, "***",
                        ifelse(dunn_df$P.adj < 0.01, "**",
                               ifelse(dunn_df$P.adj < 0.05, "*", "ns"))))
)

print(stat.test)

# Add p-vals of interest for post-Hoc Dunns using Illustrator


########## Add SCLC subtype signatures from 3 independent datasets to compare to Caris (George, Liu, Lissa) (Fig. 6e) ###########

# First, signatures from George et al. 

george<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/George_ANPYMixed_Sigs.csv")
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")
a<-george$L2FC_AvsNYP[1:147]
n<-george$L2FC_NvsAYP
p<-george$L2FC_PvsAYN[1:79]
y<-george$L2FC_YvsANP[1:108]
mixed<-george$L2FC_MixedvsPY[1:114]


# Convert to mouse
a<-subset(mouse94, mouse94$human_homolog %in% a)
a<-a$gene_name
n<-subset(mouse94, mouse94$human_homolog %in% n)
n<-n$gene_name
p<-subset(mouse94, mouse94$human_homolog %in% p)
p<-p$gene_name
y<-subset(mouse94, mouse94$human_homolog %in% y)
y<-y$gene_name
mixed<-subset(mouse94, mouse94$human_homolog %in% mixed)
mixed<-mixed$gene_name


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a),
  name = 'George_A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n),
  name = 'George_N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p),
  name = 'George_P')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y),
  name = 'George_Y')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mixed),
  name = 'George_Mixed')


# Second, signatures from Liu et al. 
liu<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Liu_ANP_ClassSigs.csv")
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")
a<-liu$log2fcAvsNP[1:70]
n<-liu$log2fcNvsAP
p<-liu$log2fcPvsAN[1:27]
mixed<-liu$log2fcMixedvP[1:114]

# Convert to mouse
a<-subset(mouse94, mouse94$human_homolog %in% a)
a<-a$gene_name
n<-subset(mouse94, mouse94$human_homolog %in% n)
n<-n$gene_name
p<-subset(mouse94, mouse94$human_homolog %in% p)
p<-p$gene_name
mixed<-subset(mouse94, mouse94$human_homolog %in% mixed)
mixed<-mixed$gene_name


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a),
  name = 'Liu_A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n),
  name = 'Liu_N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p),
  name = 'Liu_P')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mixed),
  name = 'Liu_Mixed')

# Third, signatures from Lissa et al
lissa<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/Lissa_SubtypeManual_Sigs.csv")
mouse94<-read.csv("/Users/abbieireland/Desktop/scRNAseq/Signatures/mouse94.csv")
a<-lissa$log2fc_AvNYLN[1:101]
n<-lissa$log2fc_NvAYLN[1:52]
y<-lissa$log2fc_YvANLN
ln<-lissa$log2fc_LNvANLNY[1:72]
mixed<-lissa$log2fc_MixedvYLN[1:32]

# Convert to mouse
a<-subset(mouse94, mouse94$human_homolog %in% a)
a<-a$gene_name
n<-subset(mouse94, mouse94$human_homolog %in% n)
n<-n$gene_name
ln<-subset(mouse94, mouse94$human_homolog %in% ln)
ln<-ln$gene_name
y<-subset(mouse94, mouse94$human_homolog %in% y)
y<-y$gene_name
mixed<-subset(mouse94, mouse94$human_homolog %in% mixed)
mixed<-mixed$gene_name


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a),
  name = 'Lissa_A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n),
  name = 'Lissa_N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ln),
  name = 'Lissa_LN')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y),
  name = 'Lissa_Y')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mixed),
  name = 'Lissa_Mixed')

###### Now lets see how signatures from Caris, George, Liu, and Lissa correlate with each other, and normal cell type signatures/ChIP targets in our mouse data ###
## Correlation plot for Fig. 6e ##
library(ggplot2)
library(pheatmap)

# Extract signature scores from the Seurat object
signature_matrix <- FetchData(TBO_seurat, vars = c("Caris_A1", "George_A1", "Liu_A1","Lissa_A1","ASCL1_Targets1",
                                                   "Caris_N1","George_N1", "Liu_N1","Lissa_N1", "NEUROD1_Targets1",
                                                   "Caris_P1","George_P1", "Liu_P1","POU2F3_Targets1",
                                                   "Caris_Y1","George_Y1", "Lissa_Y1", "YAP_Activity1",
                                                   "T_Cell_Inflamed_Gay1","MHC_Sig_Gay1",
                                                   "NE_Consensus1","Basal_Consensus1","Tuft_Consensus1"))


signature_matrix

# Compute Pearson correlation matrix
correlation_matrix <- cor(signature_matrix, method = "pearson")


# Define color scale from -1 to 1
breaks_seq <- seq(-1, 1, length.out = 100)  # Ensure scale covers full correlation range


# Custom softer rainbow-like palette
library(scico)
library(pheatmap)

pheatmap(correlation_matrix, 
         color = scico(100, palette = "berlin"),
         display_numbers = FALSE, 
         breaks = breaks_seq,number_color = "gray80",fontsize_number=6,
         main = "Signature correlation in mouse TBO Allografts", cluster_cols = TRUE, cluster_rows=TRUE)


######## Assess Therapeutic Targets (Fig. 6h) #########
VlnPlot(TBO_seurat,features = c("Dll3", "Ncam1", "Sez6", "Tacstd2"),group.by=c("Pheno"),cols=pheno_col,alpha=0.05, ncol=4)

# Save object with all signatures before proceeding to CellTag analyses for Fig 4 ###
saveRDS(TBO_seurat,"05_2025_RPM_RPMA_Allos_Seurat.rds")
        

##############################################################################################################
##############################################################################################################
# Begin CellTag/clonal analysis in FA space for Fig. 4 
##############################################################################################################
##############################################################################################################

# Read in previously generated seurat object
setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat")
TBO_seurat<-readRDS("05_2025_RPM_RPMA_Allos_Seurat.rds")

####### Add FA projections ######################################################################################
library(zellkonverter)
library(rhdf5)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SeuratDisk)


# Read in anndata object that includes FA projections and pseudotime values (dpt) as SCE
# Adata object previously generated with Scanpy/scVI in python as in: https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/Fig3-4_ExtFig5-6_RPM_RPMA_Allografts_CellTag_CellRank_Final_Clean.ipynb
adata2 <- readH5AD("050125_RPM_RPMA_TBOAllo_CellTagAnalysis_New_1.2_fate_FAprojection_DPT_final.h5ad")


table(adata2$UnID)
# RPMA_Allo    RPM_Allo3    RPM_Allo4 RPM_Allo_New 
# 10160         6718         6860         2880 

table(adata2$cell_type)
# ATOH1                  Basal                     NE            NE_Neuronal 
# 1343                    277                   6426                   4332 
# Neuronal SL-Epith_basal_hillock         SL-Mesenchymal         SL-stem_prolif 
# 4006                    990                   1783                   7039 
# Tuft 
# 422 

dim(assay(adata2,"counts"))
counts(adata2)<-assay(adata2,"counts")

assay(adata2,"norm")
rowData(adata2)
length(rowData(adata2)$gene_ids)

# Add dpt pseudotime coordinates to original TBO_seurat processed object
dpt<-adata2@colData$dpt_pseudotime
TBO_seurat@meta.data$dpt_pseudotime<-dpt


# Add fa embeddings from the new anndata object, adata2, to the original TBO_seurat processed object
reducedDim(adata2, "X_draw_graph_fa")
test<-reducedDim(adata2, "X_draw_graph_fa")
colnames(test)<-c("FA_1","FA_2")
rownames(test)<-colnames(TBO_seurat)
head(colnames(TBO_seurat))

TBO_seurat[['fa']] <- CreateDimReducObject(test, key="FA_", assay = "RNA")
TBO_seurat[['fa']]@cell.embeddings

# Fig. 4c UMAPs #
DimPlot(TBO_seurat,group.by='leiden_scVI_1.2',cols=my_colors, reduction='fa',label=FALSE,label.size=6)&NoAxes()
DimPlot(TBO_seurat,group.by='Pheno',cols=pheno_col, reduction='fa',label=FALSE,label.size=6)&NoAxes()
DimPlot(TBO_seurat,group.by='GenoCT',cols=my_colors, shuffle=TRUE,reduction='fa',label=FALSE,label.size=6)&NoAxes()

# Fig. 4c UMAPs split by genotype #
a<-DimPlot(TBO_seurat,group.by="leiden_scVI_1.2",cols=my_colors,reduction='fa',shuffle=TRUE,split.by="Genotype")& NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                                                 legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                 panel.background = element_rect(fill = "transparent"),
                                                                                                                                 panel.grid.major = element_blank(),
                                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                                 plot.background = element_rect(fill = "transparent",
                                                                                                                                                                color = NA))


ggsave("FAbyGeno.png", plot= a, width=9, height=5, dpi=300, bg = "transparent")

######### Annotate CellTag Data #############

md<-read.csv("TBO_metadata_withallclonedata.csv")
# Unique clone ID (remember if number is same in Allo3/4, same clone)
md$Clone_ID
#indicates all cells in clones with >5 cells for new RPM Allo3/4 or >10 for old RPMA and original RPM
md$CellTagged
ct_clone<-md$Clone_ID
ct_bin<-md$CellTagged


TBO_seurat@meta.data$CellTag_Clone<-ct_clone
table(TBO_seurat@meta.data$CellTag_Clone)
# None                  RPM_Clone_13 
# 20307                            23 
# RPM_Clone_14                  RPM_Clone_16 
# 7                             4 
# RPM_Clone_19                   RPM_Clone_2 
# 1                            23 
# RPM_Clone_22                  RPM_Clone_23 
# 3                           542 
# RPM_Clone_33                  RPM_Clone_36 
# 6                             9 
# RPM_Clone_37                  RPM_Clone_40 
# 6                             3 
# RPM_Clone_5                  RPM_Clone_53 
# 2                            10 
# RPM_Clone_55                   RPM_Clone_6 
# 18                             9 
# RPM_Clone_9 RPM_CTpostCre_Allo3_Clone_142 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_149 RPM_CTpostCre_Allo3_Clone_150 
# 2                           216 
# RPM_CTpostCre_Allo3_Clone_151 RPM_CTpostCre_Allo3_Clone_152 
# 962                            60 
# RPM_CTpostCre_Allo3_Clone_153 RPM_CTpostCre_Allo3_Clone_154 
# 3                            11 
# RPM_CTpostCre_Allo3_Clone_155 RPM_CTpostCre_Allo3_Clone_156 
# 20                             4 
# RPM_CTpostCre_Allo3_Clone_157 RPM_CTpostCre_Allo3_Clone_158 
# 8                             3 
# RPM_CTpostCre_Allo3_Clone_159 RPM_CTpostCre_Allo3_Clone_160 
# 2                            33 
# RPM_CTpostCre_Allo3_Clone_161 RPM_CTpostCre_Allo3_Clone_162 
# 22                             4 
# RPM_CTpostCre_Allo3_Clone_163 RPM_CTpostCre_Allo3_Clone_164 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_165 RPM_CTpostCre_Allo3_Clone_166 
# 2                            14 
# RPM_CTpostCre_Allo3_Clone_167 RPM_CTpostCre_Allo3_Clone_168 
# 1                             1 
# RPM_CTpostCre_Allo3_Clone_169 RPM_CTpostCre_Allo3_Clone_170 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_171 RPM_CTpostCre_Allo3_Clone_172 
# 3                            11 
# RPM_CTpostCre_Allo3_Clone_173 RPM_CTpostCre_Allo3_Clone_174 
# 6                             2 
# RPM_CTpostCre_Allo3_Clone_175 RPM_CTpostCre_Allo3_Clone_177 
# 3                             2 
# RPM_CTpostCre_Allo3_Clone_178 RPM_CTpostCre_Allo3_Clone_179 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_18 RPM_CTpostCre_Allo3_Clone_180 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_181 RPM_CTpostCre_Allo3_Clone_182 
# 5                             1 
# RPM_CTpostCre_Allo3_Clone_183 RPM_CTpostCre_Allo3_Clone_184 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_185 RPM_CTpostCre_Allo3_Clone_186 
# 6                             3 
# RPM_CTpostCre_Allo3_Clone_187 RPM_CTpostCre_Allo3_Clone_188 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_189 RPM_CTpostCre_Allo3_Clone_190 
# 6                             7 
# RPM_CTpostCre_Allo3_Clone_191 RPM_CTpostCre_Allo3_Clone_192 
# 4                             6 
# RPM_CTpostCre_Allo3_Clone_194 RPM_CTpostCre_Allo3_Clone_195 
# 3                             6 
# RPM_CTpostCre_Allo3_Clone_196 RPM_CTpostCre_Allo3_Clone_197 
# 4                             2 
# RPM_CTpostCre_Allo3_Clone_198 RPM_CTpostCre_Allo3_Clone_199 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_201 RPM_CTpostCre_Allo3_Clone_202 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_203 RPM_CTpostCre_Allo3_Clone_204 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_205 RPM_CTpostCre_Allo3_Clone_206 
# 1                             2 
# RPM_CTpostCre_Allo3_Clone_207 RPM_CTpostCre_Allo3_Clone_208 
# 4                             1 
# RPM_CTpostCre_Allo3_Clone_209 RPM_CTpostCre_Allo3_Clone_210 
# 1                             1 
# RPM_CTpostCre_Allo3_Clone_211 RPM_CTpostCre_Allo3_Clone_212 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_213 RPM_CTpostCre_Allo3_Clone_214 
# 14                             2 
# RPM_CTpostCre_Allo3_Clone_215 RPM_CTpostCre_Allo3_Clone_216 
# 3                             1 
# RPM_CTpostCre_Allo3_Clone_217 RPM_CTpostCre_Allo3_Clone_218 
# 5                             2 
# RPM_CTpostCre_Allo3_Clone_219 RPM_CTpostCre_Allo3_Clone_220 
# 1                             3 
# RPM_CTpostCre_Allo3_Clone_221 RPM_CTpostCre_Allo3_Clone_222 
# 1                             8 
# RPM_CTpostCre_Allo3_Clone_223 RPM_CTpostCre_Allo3_Clone_224 
# 3                             1 
# RPM_CTpostCre_Allo3_Clone_225 RPM_CTpostCre_Allo3_Clone_226 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_227 RPM_CTpostCre_Allo3_Clone_228 
# 10                             4 
# RPM_CTpostCre_Allo3_Clone_229 RPM_CTpostCre_Allo3_Clone_230 
# 5                             1 
# RPM_CTpostCre_Allo3_Clone_232 RPM_CTpostCre_Allo3_Clone_233 
# 3                             3 
# RPM_CTpostCre_Allo3_Clone_234 RPM_CTpostCre_Allo3_Clone_235 
# 3                             3 
# RPM_CTpostCre_Allo3_Clone_236 RPM_CTpostCre_Allo3_Clone_237 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_238 RPM_CTpostCre_Allo3_Clone_239 
# 3                             3 
# RPM_CTpostCre_Allo3_Clone_240 RPM_CTpostCre_Allo3_Clone_241 
# 4                             3 
# RPM_CTpostCre_Allo3_Clone_242 RPM_CTpostCre_Allo3_Clone_243 
# 2                             4 
# RPM_CTpostCre_Allo3_Clone_244 RPM_CTpostCre_Allo3_Clone_245 
# 4                             1 
# RPM_CTpostCre_Allo3_Clone_246 RPM_CTpostCre_Allo3_Clone_247 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_248 RPM_CTpostCre_Allo3_Clone_249 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_250 RPM_CTpostCre_Allo3_Clone_251 
# 1                             3 
# RPM_CTpostCre_Allo3_Clone_252 RPM_CTpostCre_Allo3_Clone_253 
# 2                             4 
# RPM_CTpostCre_Allo3_Clone_254 RPM_CTpostCre_Allo3_Clone_255 
# 1                             4 
# RPM_CTpostCre_Allo3_Clone_256 RPM_CTpostCre_Allo3_Clone_257 
# 6                             7 
# RPM_CTpostCre_Allo3_Clone_258 RPM_CTpostCre_Allo3_Clone_259 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_260 RPM_CTpostCre_Allo3_Clone_261 
# 4                             2 
# RPM_CTpostCre_Allo3_Clone_262 RPM_CTpostCre_Allo3_Clone_263 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_264 RPM_CTpostCre_Allo3_Clone_265 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_266 RPM_CTpostCre_Allo3_Clone_267 
# 1                             3 
# RPM_CTpostCre_Allo3_Clone_268 RPM_CTpostCre_Allo3_Clone_269 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_270 RPM_CTpostCre_Allo3_Clone_271 
# 1                             2 
# RPM_CTpostCre_Allo3_Clone_272 RPM_CTpostCre_Allo3_Clone_273 
# 2                             6 
# RPM_CTpostCre_Allo3_Clone_274 RPM_CTpostCre_Allo3_Clone_275 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_276 RPM_CTpostCre_Allo3_Clone_277 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_279  RPM_CTpostCre_Allo3_Clone_28 
# 2                            14 
# RPM_CTpostCre_Allo3_Clone_280 RPM_CTpostCre_Allo3_Clone_281 
# 1                             1 
# RPM_CTpostCre_Allo3_Clone_282 RPM_CTpostCre_Allo3_Clone_283 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_284 RPM_CTpostCre_Allo3_Clone_285 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_286 RPM_CTpostCre_Allo3_Clone_288 
# 1                             1 
# RPM_CTpostCre_Allo3_Clone_289 RPM_CTpostCre_Allo3_Clone_292 
# 4                             1 
# RPM_CTpostCre_Allo3_Clone_293 RPM_CTpostCre_Allo3_Clone_294 
# 3                             1 
# RPM_CTpostCre_Allo3_Clone_295 RPM_CTpostCre_Allo3_Clone_296 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_297 RPM_CTpostCre_Allo3_Clone_298 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_299 RPM_CTpostCre_Allo3_Clone_300 
# 1                             2 
# RPM_CTpostCre_Allo3_Clone_301 RPM_CTpostCre_Allo3_Clone_302 
# 2                             3 
# RPM_CTpostCre_Allo3_Clone_303 RPM_CTpostCre_Allo3_Clone_305 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_306 RPM_CTpostCre_Allo3_Clone_307 
# 1                             2 
# RPM_CTpostCre_Allo3_Clone_308 RPM_CTpostCre_Allo3_Clone_309 
# 2                             1 
# RPM_CTpostCre_Allo3_Clone_310 RPM_CTpostCre_Allo3_Clone_311 
# 2                             2 
# RPM_CTpostCre_Allo3_Clone_312  RPM_CTpostCre_Allo3_Clone_57 
# 2                             4 
# RPM_CTpostCre_Allo3_Clone_85 RPM_CTpostCre_Allo4_Clone_150 
# 1                             1 
# RPM_CTpostCre_Allo4_Clone_151 RPM_CTpostCre_Allo4_Clone_152 
# 1                             1 
# RPM_CTpostCre_Allo4_Clone_287 RPM_CTpostCre_Allo4_Clone_307 
# 36                             1 
# RPM_CTpostCre_Allo4_Clone_313 RPM_CTpostCre_Allo4_Clone_314 
# 60                           133 
# RPM_CTpostCre_Allo4_Clone_315 RPM_CTpostCre_Allo4_Clone_316 
# 9                            27 
# RPM_CTpostCre_Allo4_Clone_317 RPM_CTpostCre_Allo4_Clone_320 
# 30                             1 
# RPM_CTpostCre_Allo4_Clone_321 RPM_CTpostCre_Allo4_Clone_322 
# 1                            20 
# RPM_CTpostCre_Allo4_Clone_323 RPM_CTpostCre_Allo4_Clone_324 
# 10                             7 
# RPM_CTpostCre_Allo4_Clone_325 RPM_CTpostCre_Allo4_Clone_326 
# 3                             5 
# RPM_CTpostCre_Allo4_Clone_327 RPM_CTpostCre_Allo4_Clone_328 
# 3                             1 
# RPM_CTpostCre_Allo4_Clone_329 RPM_CTpostCre_Allo4_Clone_330 
# 2                             1 
# RPM_CTpostCre_Allo4_Clone_331 RPM_CTpostCre_Allo4_Clone_333 
# 1                             2 
# RPM_CTpostCre_Allo4_Clone_334 RPM_CTpostCre_Allo4_Clone_335 
# 2                             3 
# RPM_CTpostCre_Allo4_Clone_336 RPM_CTpostCre_Allo4_Clone_337 
# 2                             2 
# RPM_CTpostCre_Allo4_Clone_338 RPM_CTpostCre_Allo4_Clone_339 
# 1                             2 
# RPM_CTpostCre_Allo4_Clone_340 RPM_CTpostCre_Allo4_Clone_341 
# 2                             2 
# RPMA_Clone_1                RPMA_Clone_100 
# 2                            11 
# RPMA_Clone_11                 RPMA_Clone_12 
# 227                            69 
# RPMA_Clone_13                 RPMA_Clone_17 
# 96                            17 
# RPMA_Clone_2                 RPMA_Clone_23 
# 85                            67 
# RPMA_Clone_24                 RPMA_Clone_26 
# 31                            79 
# RPMA_Clone_28                  RPMA_Clone_3 
# 213                           105 
# RPMA_Clone_32                 RPMA_Clone_36 
# 155                            70 
# RPMA_Clone_4                 RPMA_Clone_42 
# 2                            21 
# RPMA_Clone_44                  RPMA_Clone_5 
# 13                             4 
# RPMA_Clone_53                 RPMA_Clone_55 
# 18                            55 
# RPMA_Clone_58                 RPMA_Clone_63 
# 234                           216 
# RPMA_Clone_64                 RPMA_Clone_65 
# 249                            39 
# RPMA_Clone_66                 RPMA_Clone_67 
# 130                            53 
# RPMA_Clone_68                 RPMA_Clone_69 
# 725                           171 
# RPMA_Clone_7                 RPMA_Clone_72 
# 2                            14 
# RPMA_Clone_74                 RPMA_Clone_75 
# 27                            15 
# RPMA_Clone_76                 RPMA_Clone_78 
# 35                            48 
# RPMA_Clone_8                 RPMA_Clone_81 
# 29                            18 
# RPMA_Clone_82                 RPMA_Clone_84 
# 24                            36 
# RPMA_Clone_86                 RPMA_Clone_87 
# 11                            21 
# RPMA_Clone_90                 RPMA_Clone_93 
# 17                            21 
# RPMA_Clone_95                 RPMA_Clone_96 
# 11                            25 


TBO_seurat@meta.data$CellTag_Binary<-ct_bin
table(TBO_seurat@meta.data$CellTag_Binary, TBO_seurat@meta.data$UnID)

#             RPMA_Allo RPM_Allo3 RPM_Allo4 RPM_Allo_New
# CellTagged      3511      1504       343          668
# No                 0         0        29            0
# None            6649      5214      6488         2212


########################################################################
######### Visualize CellTag Data ##################
########################################################################
cs<-table(TBO_seurat@meta.data$CellTag_Clone, TBO_seurat@meta.data$UnID)
cs
write.csv(cs,"042925_postQC_clonesize.csv")

########################################################################
######## Visualizing Celltag/clone data #########
########################################################################
#From table above, added whether clones were robust or not, now uploading that metadata
# Robust defined as >5 cells per clone post-QC
# CellTag metadata, clone information can also be found in Supplementary Table 4 of Ireland et al, 2025

md<-read.csv("TBO_metadata_withallclonedata.csv")
# Unique clone ID (remember if number is same in Allo3/4, same clone)
robust<-md$Robust


TBO_seurat@meta.data$Robust<-robust
Idents(TBO_seurat)<-'Robust'


table(TBO_seurat@meta.data$Robust, TBO_seurat@meta.data$UnID)
#          RPMA_Allo RPM_Allo3 RPM_Allo4 RPM_Allo_New
# No          6659      5244      6523         2227
# Robust      3501      1474       337          653

table(TBO_seurat@meta.data$Robust, TBO_seurat@meta.data$GenoCT)
#             RPMA_CTpreCre RPM_CTpostCre RPM_CTpreCre
# No              6659         11767         2227
# Robust          3501          1811          653

########################################################################
# Save seurat object that has FA, UMAP, and diffusion pseudotime coordinates 
saveRDS(TBO_seurat,"05_2025_RPM_RPMA_TBO_CellTag_Seurat_wSigs_FA_dpt_final.rds")
########################################################################

# Read back in later
setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag/RPMvRPMA_Seurat")
TBO_seurat<-readRDS("050225_RPM_RPMA_TBO_CellTag_Seurat_wSigs_FA_dpt_final.rds")


# Subset just clones with >5 cells to analyze clonal dynamis #
Idents(TBO_seurat)<-'Robust'
clones<-subset(TBO_seurat,idents=c("Robust"))
table(clones@meta.data$Robust, clones@meta.data$Genotype)
#       RPM   RPMA
# Robust 2464 3501


# Check to ensure only clones with >5 cells left
test<-table(clones@meta.data$CellTag_Clone)
test
# RPM_Clone_13                  RPM_Clone_14                   RPM_Clone_2 
# 23                             7                            23 
# RPM_Clone_23                  RPM_Clone_33                  RPM_Clone_36 
# 542                             6                             9 
# RPM_Clone_37                  RPM_Clone_53                  RPM_Clone_55 
# 6                            10                            18 
# RPM_Clone_6 RPM_CTpostCre_Allo3_Clone_150 RPM_CTpostCre_Allo3_Clone_151 
# 9                           216                           962 
# RPM_CTpostCre_Allo3_Clone_152 RPM_CTpostCre_Allo3_Clone_154 RPM_CTpostCre_Allo3_Clone_155 
# 60                            11                            20 
# RPM_CTpostCre_Allo3_Clone_157 RPM_CTpostCre_Allo3_Clone_160 RPM_CTpostCre_Allo3_Clone_161 
# 8                            33                            22 
# RPM_CTpostCre_Allo3_Clone_166 RPM_CTpostCre_Allo3_Clone_172 RPM_CTpostCre_Allo3_Clone_173 
# 14                            11                             6 
# RPM_CTpostCre_Allo3_Clone_181 RPM_CTpostCre_Allo3_Clone_185 RPM_CTpostCre_Allo3_Clone_189 
# 5                             6                             6 
# RPM_CTpostCre_Allo3_Clone_190 RPM_CTpostCre_Allo3_Clone_192 RPM_CTpostCre_Allo3_Clone_195 
# 7                             6                             6 
# RPM_CTpostCre_Allo3_Clone_213 RPM_CTpostCre_Allo3_Clone_217 RPM_CTpostCre_Allo3_Clone_222 
# 14                             5                             8 
# RPM_CTpostCre_Allo3_Clone_227 RPM_CTpostCre_Allo3_Clone_229 RPM_CTpostCre_Allo3_Clone_256 
# 10                             5                             6 
# RPM_CTpostCre_Allo3_Clone_257 RPM_CTpostCre_Allo3_Clone_273  RPM_CTpostCre_Allo3_Clone_28 
# 7                             6                            14 
# RPM_CTpostCre_Allo4_Clone_287 RPM_CTpostCre_Allo4_Clone_313 RPM_CTpostCre_Allo4_Clone_314 
# 36                            60                           133 
# RPM_CTpostCre_Allo4_Clone_315 RPM_CTpostCre_Allo4_Clone_316 RPM_CTpostCre_Allo4_Clone_317 
# 9                            27                            30 
# RPM_CTpostCre_Allo4_Clone_322 RPM_CTpostCre_Allo4_Clone_323 RPM_CTpostCre_Allo4_Clone_324 
# 20                            10                             7 
# RPM_CTpostCre_Allo4_Clone_326                RPMA_Clone_100                 RPMA_Clone_11 
# 5                            11                           227 
# RPMA_Clone_12                 RPMA_Clone_13                 RPMA_Clone_17 
# 69                            96                            17 
# RPMA_Clone_2                 RPMA_Clone_23                 RPMA_Clone_24 
# 85                            67                            31 
# RPMA_Clone_26                 RPMA_Clone_28                  RPMA_Clone_3 
# 79                           213                           105 
# RPMA_Clone_32                 RPMA_Clone_36                 RPMA_Clone_42 
# 155                            70                            21 
# RPMA_Clone_44                 RPMA_Clone_53                 RPMA_Clone_55 
# 13                            18                            55 
# RPMA_Clone_58                 RPMA_Clone_63                 RPMA_Clone_64 
# 234                           216                           249 
# RPMA_Clone_65                 RPMA_Clone_66                 RPMA_Clone_67 
# 39                           130                            53 
# RPMA_Clone_68                 RPMA_Clone_69                 RPMA_Clone_72 
# 725                           171                            14 
# RPMA_Clone_74                 RPMA_Clone_75                 RPMA_Clone_76 
# 27                            15                            35 
# RPMA_Clone_78                  RPMA_Clone_8                 RPMA_Clone_81 
# 48                            29                            18 
# RPMA_Clone_82                 RPMA_Clone_84                 RPMA_Clone_86 
# 24                            36                            11 
# RPMA_Clone_87                 RPMA_Clone_90                 RPMA_Clone_93 
# 21                            17                            21 
# RPMA_Clone_95                 RPMA_Clone_96 
# 11                            25 



# Assess clones by leiden cluster for Fig. 4d #

# First, just look at bar graph, no particular order #
library(ggpubr)
Idents(clones)<-'leiden_scVI_1.2'
Idents(clones)

x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))

# proportions$Cluster
colnames(proportions)<-c("Cluster", "Sample", "Frequency")
proportions
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")


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
  "#66C2A5", # teal
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
  "#66A61E"  # olive green (NOT same green as before)
)

p + scale_fill_manual(values=my_colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)


# Cluster distribution in unbiased way. 

# Uset this file to get x axis labels by removing duplicates
write.csv(p$data,"test_stacked_cluster_042925.csv")
df<-p$data
head(df)

# Convert to matrix
library(tidyr)
library(dplyr)
library(tibble)

# Pivot to wide format: Cluster × Sample
mat <- df %>%
  pivot_wider(names_from = Sample, values_from = Frequency, values_fill = 0) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

# Normalize to proportions if needed (e.g., % per bar)
mat_norm <- prop.table(mat, margin = 1)  # normalize each row to sum to 1
str(mat_norm)
dim(mat_norm)

write.csv(mat_norm,"mat__norm_test.csv")

mat_trans<-t(mat)
mat_trans

library(pheatmap)
t<-pheatmap(mat_trans,cutree_rows = 1,cutree_cols = 8,cellwidth = 5, cellheight = 5,fontsize = 10, cluster_rows=FALSE,border_color=NA,color = colorRampPalette(c("darkturquoise","black","red2"))(30))
clone_order<-colnames(mat_trans)[t$tree_col$order]
clone_order
t<-pheatmap(mat_trans,cutree_rows = 1,cutree_cols = 1,cellwidth = 5, cellheight = 5,fontsize = 10, cluster_rows=FALSE,border_color=NA,color = colorRampPalette(c("darkturquoise","black","red2"))(30))
length(clone_order)


# Now, plot clones by Leiden, ordered, final # (Fig. 4d)
Idents(clones)<-'leiden_scVI_1.2'

x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")
proportions$Cluster<-factor(proportions$Cluster, clone_order)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=my_colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)



## Annotate by group and look at dynamics ##
write.csv(clone_order,"order_clones_for_group_anno.csv")
write.csv(clones@meta.data, "new_meta_042925_patterns.csv")

# Add pattern manually to clone_order file
# Vlookup to add pattern data to metadata and read back in
md<-read.csv("new_meta_042925_patterns_rev.csv")
table(md$Pattern)
# Pattern_1 Pattern_2 Pattern_3 Pattern_4 Pattern_5 Unknown_1 Unknown_2 
# 553       675       693      2797      1032       167        48 

clones@meta.data$Clone_Dynamics<-md$Pattern
Idents(clones)<-'Clone_Dynamics'
Idents(clones)
table(clones@meta.data$Clone_Dynamics)


# Visualize clonal patterns in FA 
# Patterns 1-5, Unknown 1 = Pattern 6, Unknown 2= Pattern 7

table(clones@meta.data$Clone_Dynamics)

p1<-subset(clones,idents=c('Pattern_1'))
p2<-subset(clones,idents=c('Pattern_2'))
p3<-subset(clones,idents=c('Pattern_3'))
p4<-subset(clones,idents=c('Pattern_4'))
p5<-subset(clones,idents=c('Pattern_5'))
u1<-subset(clones,idents=c('Unknown_1'))
u2<-subset(clones,idents=c('Unknown_2'))


p1<-rownames(p1@meta.data)
p2<-rownames(p2@meta.data)
p3<-rownames(p3@meta.data)
p4<-rownames(p4@meta.data)
p5<-rownames(p5@meta.data)
u1<-rownames(u1@meta.data)
u2<-rownames(u2@meta.data)

# Patterns in FA space for Fig. 4e #

DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p1,sizes.highlight=1, cols.highlight=c("orange"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p2,sizes.highlight=1, cols.highlight=c("green2"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p3,sizes.highlight=1, cols.highlight=c("red"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p4,sizes.highlight=1, cols.highlight=c("royalblue2"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p5,sizes.highlight=1, cols.highlight=c("purple"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=u1,sizes.highlight=1, cols.highlight=c("gray40"))+ggtitle("") & NoLegend() & NoAxes()
DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=u2,sizes.highlight=1, cols.highlight=c("black"))+ggtitle("") & NoLegend() & NoAxes()



# For Fig. 4f overlap cell fate by pattern #
DimPlot(TBO_seurat,group.by="Genotype",cols=c("grey","grey"),reduction='fa')+ggtitle("") & NoLegend() & NoAxes()

# Export to overlap with plain gray fa map
p1<-subset(clones,idents=c('Pattern_1'))
p2<-subset(clones,idents=c('Pattern_2'))
p3<-subset(clones,idents=c('Pattern_3'))
p4<-subset(clones,idents=c('Pattern_4'))
p5<-subset(clones,idents=c('Pattern_5'))
u1<-subset(clones,idents=c('Unknown_1'))
u2<-subset(clones,idents=c('Unknown_2'))


# Example code, change pheno_rev color vector to match clones in each pattern
# Export transparent FA map with colored cells per pattern
# Overlay on gray FA map background in Illustrator 

pheno_rev<-c("brown2","darkorchid4","dodgerblue","turquoise4")
clusterumap<-DimPlot(u1,group.by="Pheno",reduction='fa',order=TRUE, cols=pheno_rev, pt.size=8, shuffle=FALSE)+ggtitle("") & NoAxes()&theme(legend.position="none",
                                                                                                                                                 legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                                 panel.background = element_rect(fill = "transparent"),
                                                                                                                                                 panel.grid.major = element_blank(),
                                                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                                                 plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                color = NA))


clusterumap
ggsave("Pheno_Pattern6.png", plot= clusterumap, width=5, height=4, dpi=300, bg = "transparent")

# Repeat above for each pattern, double check color vector each time


# For Fig. 4g dpt psuedotime in FA space 
library(viridis)
FeaturePlot(TBO_seurat, features = c("dpt_pseudotime"), pt.size=0.01, reduction='fa',)+scale_color_viridis(option="viridis",direction=-1)& NoAxes()

# For Fig. 4h overlaps, example for Pattern 6, repeat for each pattern and export image, overlay on gray FA in illustrator
clusterumap<-
  FeaturePlot(u1, features = c("dpt_pseudotime"), pt.size=6, reduction='fa')+scale_color_viridis(option="turbo",direction=-1)+ggtitle("")& NoAxes()&theme(legend.position="none",
                                                                                                                                              panel.background = element_rect(fill = "transparent"),
                                                                                                                                              panel.grid.major = element_blank(),
                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                              plot.background = element_rect(fill = "transparent",
                                                                                                                                                                             color = NA))

clusterumap
ggsave("p6_pseudo.png", plot= clusterumap, width=5, height=4, dpi=300, bg = "transparent")


## Visualize individual clones for Ext Data Fig. 6b,c ##

# Group them by sample..

# First, original RPM clones, pattern 1
Idents(clones)<-'Clone_Dynamics'
table(clones$Clone_Dynamics)
p1<-subset(clones,idents=c('Pattern_1'))
library(dplyr)    
test<-as.data.frame(p1$CellTag_Clone)
test$Barcodes<-rownames(test)



test <- test %>% group_by(p1$CellTag_Clone)
test<-group_split(test)
test

# Plot in for loop all RPM clones in Pattern 1
plot_lst <- vector("list", length = 26)
for (i in 1:26) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("orange"))+
    ggtitle(paste0(test[[i]]$`p1$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=5)


###### Pattern 2 clones ##########

p2<-subset(clones,idents=c('Pattern_2'))
library(dplyr)    
test<-as.data.frame(p2$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p2$CellTag_Clone)
test
test<-group_split(test)
test


DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
        cells.highlight=test[[1]]$Barcodes,sizes.highlight=2, 
        cols.highlight=c("green"))+
  ggtitle(paste0(test[[1]]$`p2$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()


plot_lst <- vector("list", length = 7)
for (i in 1:7) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("green"))+
    ggtitle(paste0(test[[i]]$`p2$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=5)

###### Pattern 3 clones ##########
p3<-subset(clones,idents=c('Pattern_3'))
library(dplyr)    
test<-as.data.frame(p3$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p3$CellTag_Clone)
test
test<-group_split(test)
test



plot_lst <- vector("list", length = 16)
for (i in 1:16) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("red"))+
    ggtitle(paste0(test[[i]]$`p3$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=8)


###### Pattern 4 clones ##########
p4<-subset(clones,idents=c('Pattern_4'))
library(dplyr)    
test<-as.data.frame(p4$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p4$CellTag_Clone)
test
test<-group_split(test)
test



plot_lst <- vector("list", length = 23)
for (i in 1:23) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("dodgerblue2"))+
    ggtitle(paste0(test[[i]]$`p4$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=8)


###### Pattern 5 clones ##########
p5<-subset(clones,idents=c('Pattern_5'))
library(dplyr)    
test<-as.data.frame(p5$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p5$CellTag_Clone)
test
test<-group_split(test)
test

plot_lst <- vector("list", length = 4)
for (i in 1:4) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("purple"))+
    ggtitle(paste0(test[[i]]$`p5$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=4)


###### Pattern 6 clones ##########
Idents(clones)
p6<-subset(clones,idents=c('Unknown_1'))
library(dplyr)    
test<-as.data.frame(p6$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p6$CellTag_Clone)
test
test<-group_split(test)
test

plot_lst <- vector("list", length = 7)
for (i in 1:7) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("gray30"))+
    ggtitle(paste0(test[[i]]$`p6$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=7)

###### Pattern 7 clones ##########
Idents(clones)
p7<-subset(clones,idents=c('Unknown_2'))
library(dplyr)    
test<-as.data.frame(p7$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(p7$CellTag_Clone)
test
test<-group_split(test)
test

plot_lst <- vector("list", length = 3)
for (i in 1:3) {
  g<-DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, 
             cells.highlight=test[[i]]$Barcodes,sizes.highlight=2, 
             cols.highlight=c("yellow"))+
    ggtitle(paste0(test[[i]]$`p7$CellTag_Clone`[1]))+theme(plot.title = element_text(size = 8,face = "plain")) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst, ncol=3)

# DONE #



## Visualize only clones matching in vivo in FA projection, for Ext Data Fig. 7f
Idents(clones)<-"CellTag_Clone"
# Subset just the clones that match the in vitro, pattern 1
invitromatch<-subset(clones,idents=c("RPM_Clone_14","RPM_Clone_2","RPM_Clone_33","RPM_Clone_36","RPM_Clone_6"))
table(invitromatch$CellTag_Clone)
table(invitromatch$Clone_Dynamics)
p1_cells<-colnames(invitromatch)

DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=p1_cells,sizes.highlight=2, cols.highlight=c("orange"))+ggtitle("") & NoLegend() & NoAxes()


# Save clones as seurat object
saveRDS(clones,"05_2025_RPM_RPMA_TBOAllo_CellTagClones_Onlyclones.rds")
# Read back in later
clones<-readRDS("05_2025_RPM_RPMA_TBOAllo_CellTagClones_Onlyclones.rds")
