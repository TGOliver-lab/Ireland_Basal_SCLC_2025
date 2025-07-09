# RPM and RPMA +/- Cre In Vitro to In Vivo CellTag/Clonal Analyses
# Ireland et al, 2025
# Ext. Data Fig. 7a-f


# First, set working directory
setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/030525_CellTag_InVitroInVivo_final")

# If using Isilon share
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/ExtFig7_RPM_RPMA_InVitroClones")


# Load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratData)
  library(SeuratDisk)
  library(SeuratObject)
  library(CellTagR)
  library(zellkonverter)
})


############################## CELLTAG EXPERIMENT KEY ##########################
# Matches metadata on GEO deposit GSE279200

# CellTag Experiment "A"= RPM matching organoids and resulting allograft, CellTagged pre-Cre

# CellTag Experiment "B"= RPM matching organoids and resulting allografts (n=2), CellTagged post-Cre; 

# CellTag Experiment "C"= RPMA matching organoids and resulting allograft, CellTagged pre-Cre.
############################################################################################################## First, analyze CellTag data from RPM organoids and allografts CellTaged post-Cre (Experiment B)  "#############################################

########################################################################################################
## Now, bring new WT+RPM+RPMA basal organoid object into R as seurat object to label CellTag info 
# All data for Ext Fig 7 are from "CellTag pre-Cre" experiments (CellTag experiments 'A' and 'C')
########################################################################################################
# Read in adata object as SCE, anndata was previously generated using Scanpy/scVI in python as in: https://github.com/TGOliver-lab/Ireland_Basal_SCLC_2025/blob/main/Python_Code/ExtFig7_RPM_RPMA_Organoid_Allo_Final_Clean.ipynb

adata<-readH5AD("030525_RPM_RPMA_WT_Organoids_forCellTagwStates.h5ad")

# Check object to ensure proper dataset
table(adata$Cre)
# Cre No_Cre 
# 7322   1767 

table(adata$leiden_scVI_1.2)
# 0    1    2    3    4    5    6    7 
# 1705 1624 1574 1572 1427  577  450  160 

table(adata$UnID)
# WT_Org_NoCre  RPM_Org_Cre RPMA_Org_Cre 
# 1767         1191         6131 

dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)


#Convert SCE to seurat
TBO_seurat_invitro <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
TBO_seurat_invitro
# An object of class Seurat 
# 55800 features across 9089 samples within 1 assay 
# Active assay: RNA (55800 features, 0 variable features)
# 1 layer present: counts

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

TBO_seurat_invitro
TBO_seurat_invitro[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(TBO_seurat_invitro)<-'norm'
table(TBO_seurat_invitro@meta.data$UnID)
# WT_Org_NoCre  RPM_Org_Cre RPMA_Org_Cre 
# 1767         1191         6131 

# Add embeddings of umap, and X_scVI_1.2 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(TBO_seurat_invitro)
dim(test)
TBO_seurat_invitro[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
TBO_seurat_invitro[['umap']]@cell.embeddings


# Define color vector to plot UMAP by leiden cluster for Ext Data Fig. 7b (Ext. Data Fig. 7a from Python analyses)
colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#4ba4d3', '#b5d7e4','#ef9171','#a6eb9a')
DimPlot(TBO_seurat_invitro,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=TRUE,label.size=8)&NoAxes()

# View assigned basal states (in Python) in UMAP as in Ext. Data Fig. 7d
table(TBO_seurat_invitro@meta.data$CellType_Sig_Based)
# Proteosomal_Basal     Differentiating_Basal Luminal_Hillock/Secretory       Krt8_Interm/Hillock 
# 2155                      2201                      1574                      1732 
# Proliferative_Basal 
# 1427 

DimPlot(TBO_seurat_invitro,group.by='CellType_Sig_Based',cols="Dark2", reduction='umap',label=FALSE,label.size=6)&NoAxes()


# Now, add CellTag annotations
write.csv(TBO_seurat_invitro@meta.data, "seurat_metadata_RPM_RPMA_TBO_InVitro.csv")

# To do this, in above file, add CellTag Clone ID from celltag metadata stored in Supplementary Table 4 in Ireland et al, 2025 
# Or if using lab drive...
# For RPM: "./collapse_0.8_clonecomp_RPM_Final.csv"
# For RPMA: "./collapse_0.8_clonecomp_RPMA_Final.csv"

# Also annotate corresponding in vivo patterns (as defined in Fig. 4e) using this csv file and vlookup or similar
# ./cloneswmatch.csv

md<-read.csv("051525_seurat_metadata_RPM_RPMA_TBO_InVitro_wCellTagAnno_Unlinked.csv")
ct_clone<-md$CellTag_Clone
ct_clone
pattern<-md$Patterns
ct_bin<-md$CellTag_Binary

TBO_seurat_invitro@meta.data$CellTag_Clone<-ct_clone
table(TBO_seurat_invitro@meta.data$CellTag_Clone)
# None  RPM_Clone_13  RPM_Clone_14  RPM_Clone_16  RPM_Clone_19   RPM_Clone_2  RPM_Clone_22  RPM_Clone_23 
# 6619            53             5            11            38           178            16            11 
# RPM_Clone_33  RPM_Clone_36  RPM_Clone_37  RPM_Clone_40   RPM_Clone_5  RPM_Clone_53  RPM_Clone_55   RPM_Clone_6 
# 5             6             2             5             7             1             1            11 
# RPM_Clone_7   RPM_Clone_9  RPMA_Clone_1 RPMA_Clone_11 RPMA_Clone_12 RPMA_Clone_13 RPMA_Clone_14 RPMA_Clone_17 
# 14            10            94            11            18            19            23            47 
# RPMA_Clone_18  RPMA_Clone_2 RPMA_Clone_22 RPMA_Clone_23 RPMA_Clone_24 RPMA_Clone_25 RPMA_Clone_26 RPMA_Clone_28 
# 18          1241            17             6            14             9             6             2 
# RPMA_Clone_29  RPMA_Clone_3 RPMA_Clone_30 RPMA_Clone_32 RPMA_Clone_36  RPMA_Clone_4 RPMA_Clone_42 RPMA_Clone_44 
# 31           230            13             1             1           156             1             3 
# RPMA_Clone_5 RPMA_Clone_53 RPMA_Clone_55 RPMA_Clone_58  RPMA_Clone_6  RPMA_Clone_7  RPMA_Clone_8 
# 20             1             1             2            12            47            52 

TBO_seurat_invitro@meta.data$CellTag_Pattern<-pattern
table(TBO_seurat_invitro@meta.data$CellTag_Pattern, TBO_seurat_invitro@meta.data$UnID)
#               WT_Org_NoCre RPM_Org_Cre RPMA_Org_Cre
# None              1767         918         4475
# Pattern_1            0         207            0
# Pattern_2            0          13            0
# Pattern_4            0           0         1656
# Pattern_5            0          53            0

table(TBO_seurat_invitro@meta.data$CellTag_Clone,TBO_seurat_invitro@meta.data$CellTag_Pattern)

# Save and read in later
saveRDS(TBO_seurat_invitro,"05.15.25_RPM_RPMA_WT_Organoids_wCT_Clones.rds")
TBO_seurat_invitro<-readRDS("05.15.25_RPM_RPMA_WT_Organoids_wCT_Clones.rds")


##############################################################################################################
# Begin clonal analysis in UMAP space for in vitro clones only, to start 
########################################################################################
TBO_seurat_invitro@meta.data$CellTag_Binary<- ct_bin
table(TBO_seurat_invitro@meta.data$CellTag_Binary)  
# Celltagged       None 
# 2470       6619 

Idents(TBO_seurat_invitro)<-'CellTag_Binary'

# Subset seurat object to only those with CellTags
clones<-subset(TBO_seurat_invitro,idents=c("Celltagged"))

table(clones@meta.data$UnID)
# RPM_Org_Cre RPMA_Org_Cre 
# 374         2096 

# Generate and write a table of clone size
test<-table(clones@meta.data$CellTag_Clone)
test
write.csv(test,"celltag_invitro_clones.csv")
Idents(clones)<-'CellTag_Clone'

# Label clones with > or = 5 cells as 'Robust', otherwise 'Not_Robust'
keep<-read.csv("invitro_5ormore_clones.csv")
keep<-keep$clones.5

table(clones@meta.data$CellTag_Clone %in% keep)

clones@meta.data$Robust<-ifelse(clones@meta.data$CellTag_Clone %in% keep, "Robust","Not_Robust")
table(clones@meta.data$Robust)
# Not_Robust     Robust 
# 16       2454 

# Subset seurat object to just Robust clones
Idents(clones)<-'Robust'
clones<-subset(clones,idents=c("Robust"))

# Here are the clones 
table(clones@meta.data$CellTag_Clone)
# RPM_Clone_13  RPM_Clone_14  RPM_Clone_16  RPM_Clone_19   RPM_Clone_2  RPM_Clone_22  RPM_Clone_23  RPM_Clone_33 
# 53             5            11            38           178            16            11             5 
# RPM_Clone_36  RPM_Clone_40   RPM_Clone_5   RPM_Clone_6   RPM_Clone_7   RPM_Clone_9  RPMA_Clone_1 RPMA_Clone_11 
# 6             5             7            11            14            10            94            11 
# RPMA_Clone_12 RPMA_Clone_13 RPMA_Clone_14 RPMA_Clone_17 RPMA_Clone_18  RPMA_Clone_2 RPMA_Clone_22 RPMA_Clone_23 
# 18            19            23            47            18          1241            17             6 
# RPMA_Clone_24 RPMA_Clone_25 RPMA_Clone_26 RPMA_Clone_29  RPMA_Clone_3 RPMA_Clone_30  RPMA_Clone_4  RPMA_Clone_5 
# 14             9             6            31           230            13           156            20 
# RPMA_Clone_6  RPMA_Clone_7  RPMA_Clone_8 
# 12            47            52 


# Here are their resulting in vivo 'Patterns'
table(clones@meta.data$CellTag_Pattern, clones@meta.data$UnID)

#             RPM_Org_Cre  RPMA_Org_Cre
# None              101          440
# Pattern_1         205            0
# Pattern_2          11            0
# Pattern_4           0         1644
# Pattern_5          53            0


table(clones@meta.data$CellTag_Pattern, clones@meta.data$CellTag_Clone)

# For plotting barplot in Ext. Fig. 7e
#Extract current order of clones
test<-table(clones@meta.data$CellTag_Pattern, clones@meta.data$CellTag_Clone)
colnames(test)
write.csv(colnames(test),"old_order.csv")

#Desired order of clones (obtain from manual sort based on group patterns)
order<-read.csv("new_order.csv")
order<-order$Clone
order

clones@meta.data$CellTag_Clone<-factor(clones@meta.data$CellTag_Clone, order)
table(clones@meta.data$CellTag_Clone)

# Subset clones to just those with matching in vivo clones
match<-read.csv("cloneswmatch.csv")
match$Pattern
clone_match<-match$Clone
clone_match
Idents(clones)<-'CellTag_Clone'
Idents(clones)
clones_subset<-subset(clones,idents=clone_match)
table(clones_subset$CellTag_Pattern)
# Pattern_1 Pattern_2 Pattern_4 Pattern_5 
# 205        11      1644        53 

order<-match$Clone
order

#Plot clones by Leiden cluster, ordered by in vivo pattern as above#

Idents(clones_subset)<-'CellType_Sig_Based'


x<-table(clones_subset@meta.data$CellTag_Clone,Idents(clones_subset))
proportions <- as.data.frame(100*prop.table(x, margin = 1))

colnames(proportions)<-c("Cluster", "Sample", "Frequency")
proportions
proportions$Cluster<-factor(proportions$Cluster, order)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")


# Ext. Data Fig. 7e #
colors=c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666')
p + scale_fill_manual(values=colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)


# Add pattern annotation 
x<-match$Pattern
x
# Create annotation data
annotation_df <- data.frame(
  Sample = levels(proportions$Cluster),  # match x-axis labels
  Group = x,  # Define annotation groups
  y = -0.5  # Place below x-axis, adjust as needed
)
annotation_df
# Optional: Define colors for annotation bar
annotation_colors <- c("Pattern_1" = "orange", "Pattern_2" = "green", "Pattern_5"="purple","Pattern_4"="blue")

# Add to original plot, p, for final Ext. Data Fig. 7e. Add genotype information manually in Illustrator.
library(ggnewscale)
p +
  # Original plot
  scale_fill_manual(values = colors, name="Basal cell state") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18)
  ) +
  rotate_x_text(size = 7, angle = 90) +
  
  # Reset fill scale so annotation can have its own
  new_scale_fill() +
  
  # Add annotation bar
  geom_tile(
    data = annotation_df,
    aes(x = Sample, y = y, fill = Group),
    height = 3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(name = "In vivo clonal pattern", values = annotation_colors,
                    guide = guide_legend(override.aes = list(size = 5)))



# Visualize clonal patterns in UMAP

Idents(clones_subset)<-'CellTag_Pattern'
Idents(clones_subset)

#Plot RPM patterns 1, 2, and 5 that have matching in vivo clones
p1<-subset(clones_subset,idents=c('Pattern_1'))
p2<-subset(clones_subset,idents=c('Pattern_2'))
p5<-subset(clones_subset,idents=c('Pattern_5'))

table(Idents(clones_subset))

# For Ext. Data Fig. 7f, overlap
# Plot and export clear cells, colored by basal cell state, to overlap with plain gray fa map in Illustrator
# May have to adjust width and height parameters for best overlay fit.

#####  Pattern 1 ##### 
colors=c('#1B9E77', '#D95F02', '#7570B3', '#66A61E')
clusterumap<-DimPlot(p1,group.by="CellType_Sig_Based",reduction='umap',order=TRUE, cols=colors, pt.size=6, shuffle=TRUE)+ggtitle("Pattern 1") & NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                                                                 legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                                 panel.background = element_rect(fill = "transparent"),
                                                                                                                                                 panel.grid.major = element_blank(),
                                                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                                                 plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                color = NA))


clusterumap
ggsave("UMAPp1_finalv2.png", plot= clusterumap, width=10, height=5, dpi=300, bg = "transparent")

##### Pattern 2 ##### 
colors=c('#D95F02', '#7570B3','#66A61E')
clusterumap<-DimPlot(p2,group.by="CellType_Sig_Based",reduction='umap',order=TRUE, cols=colors, pt.size=6, shuffle=TRUE)+ggtitle("Pattern 2") & NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                                                                               legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                                               panel.background = element_rect(fill = "transparent"),
                                                                                                                                                               panel.grid.major = element_blank(),
                                                                                                                                                               panel.grid.minor = element_blank(),
                                                                                                                                                               plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                              color = NA))


clusterumap
ggsave("UMAPp2_finalv2.png", plot= clusterumap, width=10, height=5, dpi=300, bg = "transparent")


#####  Pattern 5 ##### 
colors=c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666')

clusterumap<-DimPlot(p5,group.by="CellType_Sig_Based",reduction='umap',order=TRUE, cols=colors, pt.size=6, shuffle=TRUE)+ggtitle("Pattern 5") & NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                                                                               legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                                               panel.background = element_rect(fill = "transparent"),
                                                                                                                                                               panel.grid.major = element_blank(),
                                                                                                                                                               panel.grid.minor = element_blank(),
                                                                                                                                                               plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                              color = NA))


clusterumap
ggsave("UMAPp5.png", plot= clusterumap, width=10, height=5, dpi=300, bg = "transparent")



# Use the following UMAPs as templates to overlay cells colored by fate, in Illustrator
p1<-rownames(p1@meta.data)
p2<-rownames(p2@meta.data)
p5<-rownames(p5@meta.data)

DimPlot(TBO_seurat_invitro,group.by="CellTag_Clone",reduction='umap',order=TRUE, cells.highlight=p2,sizes.highlight=4, cols.highlight=c("green2"))+ggtitle("Pattern 2") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_invitro,group.by="CellTag_Clone",reduction='umap',order=TRUE, cells.highlight=p5,sizes.highlight=3, cols.highlight=c("darkorchid3"))+ggtitle("Pattern 5") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_invitro,group.by="CellTag_Clone",reduction='umap',order=TRUE, cells.highlight=p1,sizes.highlight=3, cols.highlight=c("orange"))+ggtitle("Pattern 1") & NoLegend() & NoAxes()

# NOTE: To generate corresponding FA maps in Ext. Data Fig. 7f that just include the clones also captured in vitro, see the end of R script 'Fig3_4_6_ExtFig5-6-10_RPM_RPMA_Allos_CellTag.R'

# Save objects
saveRDS(clones,"05.18.25_RPM_RPMA_invitro_CTclones_all.rds")
saveRDS(clones_subset,"05.18.25_RPM_RPMA_invitro_CTclones_matchinginvivo.rds")

# End of script, any additional formatting done in Illustrator

