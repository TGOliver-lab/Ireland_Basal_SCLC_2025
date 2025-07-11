# CellTag and Clone Calling for RPM (CellTagged pre-Cre and CellTagged post-Cre) and RPMA basal-organoid-derived allograft tumours
# To generate CellTag metadata/annotations
# Ireland et al, 2025
# Related to *Old Fig. 4c-j, and Extended Data Fig. 6b-c*
# Related to *Final Fig. 3c-j and Extended Data Fig. 6b-c*

### Load packages##

suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(CellTagR)
})

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

library("devtools")
devtools::install_github("morris-lab/CellTagR")

############################## CELLTAG EXPERIMENT KEY ##########################
# Matches metadata on GEO deposit GSE279200

# CellTag Experiment "A"= RPM matching organoids and resulting allograft, CellTagged pre-Cre

# CellTag Experiment "B"= RPM matching organoids and resulting allografts (n=2), CellTagged post-Cre; 

# CellTag Experiment "C"= RPMA matching organoids and resulting allograft, CellTagged pre-Cre.
############################################################################################################## First, analyze CellTag data from RPM organoids and allografts CellTaged post-Cre (Experiment B)  "#############################################

########################################################################################
# CellTag extraction/Clone calling performed following published workflow on Github at: https://github.com/morris-lab/CellTagR  

########################################################################################
# First, ID clones in RPM organoids and allografts CellTagged post-Cre (Experiment "B")
########################################################################################

setwd("/Users/abbieireland/Desktop/scRNAseq/042525_RPMTBO_CTpostCre_CellTag")

# If using Isilon share
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig3-4-6_ExtFig5-6-10_RPM_RPMA_Allos_UMAP_FA_CellTag/CellTag_Clone_Calling/042525_RPMTBO_CTpostCre_CellTag")

#Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "./bams")
bam.test.obj
#Sample-1 is In Vitro Organoids + CMV Cre post CellTag
#Sample-2 is RPM "Allo3" (called Allo1 in bam) In Vivo (CT post-Cre)
#Sample-3 is "RPM Allo4" (called Allo2 in bam) In Vivo (CT post-Cre)

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
bam.test.obj

# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
Barcode.Aggregate(list("./barcodes/1-invitro-barcodes.tsv", "./barcodes/2-Allo1-barcodes.tsv","./barcodes/3-Allo2-barcodes.tsv"), "./barcodes/all_barcodes.tsv")


# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file ="./barcodes/all_barcodes.tsv")

bam.test.obj
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  7682 
# Raw Number of Cells with CellTag =  31302 
# Collapsed CellTag Counts =  0 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)
tail(colnames(bam.test.obj@raw.count))
head(bam.test.obj@raw.count)
head(bam.test.obj@bam.parse.rslt)


# Generating the collapsing file
bam.test.obj.collapsed <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "./collapsing.txt")


#Run starcode in terminal per sample
# starcode/starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

collapsed.rslt.dir <- "./collapsing_results"

# Recount and generate collapsed matrix
bam.test.obj.collapsed <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj.collapsed, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
bam.test.obj.collapsed
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  7682 
# Raw Number of Cells with CellTag =  31302 
# Collapsed CellTag Counts =  6279 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)
bam.test.obj.collapsed@collapsed.count
bam.test.obj.collapsed@raw.count

# Calling binarization
bam.test.obj.collapsed <- SingleCellDataBinarization(bam.test.obj.collapsed, 2)

# Read the RDS file and get the object
dt.mtx.whitelist.path <- "./v1_whitelist_ours_082023_75percentile.csv"
bam.test.obj.collapsed <- SingleCellDataWhitelist(bam.test.obj.collapsed, dt.mtx.whitelist.path)
bam.test.obj.collapsed

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  7682 
# Raw Number of Cells with CellTag =  31302 
# Collapsed CellTag Counts =  6279 
# Whitelisted CellTag Counts =  6265 
# Whitelisted Number of Cells with CellTag =  20685 

MetricPlots(bam.test.obj.collapsed)
# Average:  2.020546 
# Frequency:  6.671189 

bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 20, comparison = "less")
bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 2, comparison = "greater")

MetricPlots(bam.test.obj.collapsed)
# Average:  5.277042 
# Frequency:  5.022666 

bam.test.obj.collapsed
bam.test.obj.collapsed <- JaccardAnalysis(bam.test.obj.collapsed, fast = T)


# Call clones
bam.test.obj.collapsed <- CloneCalling(celltag.obj = bam.test.obj.collapsed, correlation.cutoff=0.8)

# Check them out!!
bam.test.obj.collapsed@clone.composition[["v1"]]
bam.test.obj.collapsed@clone.size.info[["v1"]]

comp<-bam.test.obj.collapsed@clone.composition[["v1"]]
size<-bam.test.obj.collapsed@clone.size.info[["v1"]]
tail(bam.test.obj.collapsed@whitelisted.count)

write.csv(comp,"./collapse_0.8_clonecomp.csv")
write.csv(size,"./collapse_0.8_clonesize.csv")


## Add a column called Sample_ID and then read in comp data to calculate quickly comp/clone 

comp_clone<-read.csv("collapse_0.8_clonecomp.csv")
tail(comp_clone)
table(comp_clone$clone.id,comp_clone$Sample_ID)

table(comp_clone$Sample_ID)
# RPM_Allo1       RPM_Allo2     RPM_InVitro_TBO 
# 2038             456             384 

df<-table(comp_clone$clone.id, comp_clone$Sample_ID)
df
write.csv(df, "compperclonebySample.csv")


## Make an alluvial plot to track bottlenecks/changes in celltag populations over each sample
library(ggalluvial)
library(tidyverse)

### Showing All clones here 
db3<-read.csv("percent_clone_persample.csv")
db3<-data.frame(db3)
data_bar3<-db3
data_bar3
table(data_bar3$Sample)
# RPM_Allo1       RPM_Allo2 RPM_InVitro_TBO 
# 341             341             341 
data_bar3$Sample<-factor(data_bar3$Sample, c("RPM_InVitro_TBO","RPM_Allo1","RPM_Allo2"))


## first need to calculate the actual percentage by group/Sample
data_bar3 %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .2) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="none")


### Showing All clones, no ZEROS
db3<-read.csv("percent_clone_persample_NOZERO.csv")
db3<-data.frame(db3)
data_bar3<-db3
data_bar3
table(data_bar3$Sample)
# RPM_InVitro_TBO       RPM_Allo1       RPM_Allo2 
# 148             169              34 
data_bar3$Sample<-factor(data_bar3$Sample, c("RPM_InVitro_TBO","RPM_Allo2","RPM_Allo1"))


## Plot alluvial plot
data_bar3 %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .2) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="none")


########################################################################################
# Next, ID CellTags in  RPM organoids and allograft CellTagged Pre-Cre (Experiment "A") 
########################################################################################
setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPM_TBO_Allo2_CellTagAnalysis/0531_InVitroCretoInVivo")

# If using Isilon share...
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig3-4-6_ExtFig5-6-10_RPM_RPMA_Allos_UMAP_FA_CellTag/CellTag_Clone_Calling/0531_RPM_InVitroCretoInVivo_CTpreCre")


#Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "./bams")
bam.test.obj
#Sample-1 is In Vitro + CT pre-CMV Cre
#Sample-2 is RPM Allo In Vivo, CT pre-CMV-Cre

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
bam.test.obj

# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
Barcode.Aggregate(list("./barcodes/0531_RPM_InVitro_Cre_barcodes.tsv", "./barcodes/0531_RPM_TBO_barcodes.tsv"), "./barcodes/all_barcodes.tsv")


# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file ="./barcodes/all_barcodes.tsv")

bam.test.obj
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  0 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)
tail(colnames(bam.test.obj@raw.count))
head(bam.test.obj@raw.count)
head(bam.test.obj@bam.parse.rslt)



##The collapsing steps are actually optional.. try with first ###
# Generating the collapsing file
bam.test.obj.collapsed <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "./collapsing.txt")


#Run starcode in terminal per sample
# starcode/starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

collapsed.rslt.dir <- "./collapsing_results/"

# Recount and generate collapsed matrix
bam.test.obj.collapsed <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj.collapsed, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
bam.test.obj.collapsed
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  1447 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)
bam.test.obj.collapsed@collapsed.count
bam.test.obj.collapsed@raw.count


# Calling binarization
bam.test.obj.collapsed <- SingleCellDataBinarization(bam.test.obj.collapsed, 2)


# Read the RDS file and get the object
dt.mtx.whitelist.path <- "./v1_whitelist_ours_082023_75percentile.csv"
bam.test.obj.collapsed <- SingleCellDataWhitelist(bam.test.obj.collapsed, dt.mtx.whitelist.path)
bam.test.obj.collapsed

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  1447 
# Whitelisted CellTag Counts =  1427 
# Whitelisted Number of Cells with CellTag =  4840 

MetricPlots(bam.test.obj.collapsed)
# Average:  1.542149 
# Frequency:  5.230554 

bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 20, comparison = "less")
bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 2, comparison = "greater")

MetricPlots(bam.test.obj.collapsed)

# Average:  2.705283 
# Frequency:  4.0911
bam.test.obj.collapsed
bam.test.obj.collapsed <- JaccardAnalysis(bam.test.obj.collapsed, fast = T)


# Call clones
bam.test.obj.collapsed <- CloneCalling(celltag.obj = bam.test.obj.collapsed, correlation.cutoff=0.8)

# Check them out!!
bam.test.obj.collapsed@clone.composition[["v1"]]
bam.test.obj.collapsed@clone.size.info[["v1"]]

comp<-bam.test.obj.collapsed@clone.composition[["v1"]]
size<-bam.test.obj.collapsed@clone.size.info[["v1"]]
tail(bam.test.obj.collapsed@whitelisted.count)

write.csv(comp,"./collapse_0.8_clonecomp.csv")
write.csv(size,"./collapse_0.8_clonesize.csv")


## Read in comp data to calculate quickly comp/clone 

comp_clone<-read.csv("collapse_0.8_clonecomp.csv")
head(comp_clone)
table(comp_clone$clone.id,comp_clone$Sample_ID)

table(comp_clone$Sample_ID)
# Allograft InVitro_Cre 
# 1488         516 

df<-table(comp_clone$clone.id, comp_clone$Sample_ID)
df
write.csv(df, "compperclonebySample.csv")


###############################################################################################################
######## Last, ID CellTags/clones in RPMA organoids and resulting allograft (Experiment "C")
########################################################################################################################
########################################################################################################################

setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPMA_TBO_scRNAseq_031824/CellTagAnalysis/InVitroCretoAllo")
# If using Isilon share...
# setwd("All_Staff/Abbie/Basal_Lung_Manuscript/scRNAseq/Final_R/Fig3-4-6_ExtFig5-6-10_RPM_RPMA_Allos_UMAP_FA_CellTag/CellTag_Clone_Calling/RPMA_TBO_scRNAseq_031824")

library("devtools")
devtools::install_github("morris-lab/CellTagR")
library("CellTagR")


#Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "./bams")
bam.test.obj
#Sample-1 is In Vitro Organoids, CellTagged pre-Cre 
#Sample-2 is RPMA Allo In Vivo

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
bam.test.obj

# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
Barcode.Aggregate(list("./barcodes/1-RPMA_TBO_CMV_barcodes.tsv", "./barcodes/2-RPMA_Allo_barcodes.tsv"), "./barcodes/all_barcodes.tsv")


# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file ="./barcodes/all_barcodes.tsv")
bam.test.obj

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  0 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)
tail(colnames(bam.test.obj@raw.count))
head(bam.test.obj@raw.count)
head(bam.test.obj@bam.parse.rslt)


# Generating the collapsing file
bam.test.obj.collapsed <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "collapsing.txt")


#Run starcode in terminal per sample
# starcode/starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

collapsed.rslt.dir <- "collapsing_results"

# Recount and generate collapsed matrix
bam.test.obj.collapsed <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj.collapsed, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
bam.test.obj.collapsed
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  4188 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)
bam.test.obj.collapsed@collapsed.count
bam.test.obj.collapsed@raw.count


# Calling binarization
bam.test.obj.collapsed <- SingleCellDataBinarization(bam.test.obj.collapsed, 2)
bam.test.obj.collapsed

# Read the RDS file and get the object
dt.mtx.whitelist.path <- "./v1_whitelist_ours_082023_75percentile.csv"
bam.test.obj.collapsed <- SingleCellDataWhitelist(bam.test.obj.collapsed, dt.mtx.whitelist.path)
bam.test.obj.collapsed

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  4188 
# Whitelisted CellTag Counts =  4090 
# Whitelisted Number of Cells with CellTag =  15167 

MetricPlots(bam.test.obj.collapsed)
# Average:  1.786115 
# Frequency:  6.623472 

bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 20, comparison = "less")
bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 2, comparison = "greater")

MetricPlots(bam.test.obj.collapsed)

# W collapse 75 percentile of celltag whitelist
# Average:  3.480637 
# Frequency:  5.559658 

bam.test.obj.collapsed
bam.test.obj.collapsed <- JaccardAnalysis(bam.test.obj.collapsed, fast = T)


# Call clones
bam.test.obj.collapsed <- CloneCalling(celltag.obj = bam.test.obj.collapsed, correlation.cutoff=0.8)

# Check them out!!
bam.test.obj.collapsed@clone.composition[["v1"]]
bam.test.obj.collapsed@clone.size.info[["v1"]]

comp<-bam.test.obj.collapsed@clone.composition[["v1"]]
size<-bam.test.obj.collapsed@clone.size.info[["v1"]]
tail(bam.test.obj.collapsed@whitelisted.count)

write.csv(comp,"./collapse_0.8_clonecomp.csv")
write.csv(size,"./collapse_0.8_clonesize.csv")


### Read in comp data to calculate quickly comp/clone ###

getwd()
comp_clone<-read.csv("collapse_0.8_clonecomp.csv")

table(comp_clone$Sample)
# Allograft InVitro_Cre 
# 3916        2321 

df<-table(comp_clone$clone.id, comp_clone$Sample)
df
write.csv(df, "RPMA_compperclonebySample.csv")


## End of script, now use the CellTag metadata generated in this script to annotate the RPM/RPMA basal allograft Seurat object and perform additional clonal analyses ##
## To do this, see script, "Fig3_4_6_ExtFig5-6-10_RPM_RPMA_Allos_CellTag.R" ##
