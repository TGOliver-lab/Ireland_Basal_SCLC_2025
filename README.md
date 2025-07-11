# Ireland_Basal_SCLC_2025
#### Repository for code used in Ireland et al, _In Press at Nature_, 2025. 

## **Basal cell of origin resolves neuroendocrine-tuft lineage plasticity in cancer**

#### Authors:
Abbie S. Ireland<sup>1</sup>, Daniel Xie<sup>1</sup>, Sarah B. Hawgood<sup>1</sup>,  Margaret W. Barbier<sup>1</sup>, Scarlett Lucas-Randolph<sup>1</sup>, Darren R. Tyson<sup>1</sup>, Lisa Y. Zuo<sup>1</sup>, Benjamin E. Hanna<sup>1</sup>, Benjamin L. Witt<sup>2</sup>, Ramaswamy Govindan<sup>3</sup>, Afshin Dowlati<sup>4</sup>, Justin C. Moser<sup>5</sup>, Anish Thomas<sup>6</sup>, Sonam Puri<sup>7</sup>, Charles M. Rudin<sup>8</sup>, Joseph M. Chan<sup>9</sup>, Andrew Elliott<sup>10</sup>, Trudy G. Oliver<sup>1*</sup>

#### Affiliations:	
<sup>1</sup>Department of Pharmacology and Cancer Biology, Duke University, Durham, NC, 27710, USA;
<sup>2</sup>Department of Pathology, University of Utah, Salt Lake City, UT 84112, USA; 
<sup>3</sup>Division of Oncology, Department of Medicine, Alvin J. Siteman Cancer Center, Washington University School of Medicine, St Louis, MO, 63110, USA;
<sup>4</sup>Division of Hematology and Oncology, Department of Medicine, University Hospitals Seidman Cancer Center and Case Western Reserve University, Cleveland, OH, 44106, USA;
<sup>5</sup>HonorHealth Research Institute, Scottsdale, AZ 85254, USA;
<sup>6</sup>Developmental Therapeutics Branch, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, 20892, USA;
<sup>7</sup>Department of Thoracic Oncology, Moffitt Cancer Center, Tampa, FL 33612, USA;
<sup>8</sup>Department of Medicine, Memorial Sloan Kettering Cancer Center, New York, NY 10065, USA;
<sup>9</sup>Human Oncology and Pathogenesis Program, Memorial Sloan Kettering Cancer Center, New York, NY, 10065, USA;
<sup>10</sup>Caris Lifes Sciences, Tempe, AZ, 85040, USA.
<sup>*</sup>Corresponding author. Email: tgo@duke.edu

## Summary paragraph:
Neuroendocrine and tuft cells are rare, chemosensory epithelial lineages defined by expression of ASCL1 and POU2F3 transcription factors, respectively. Neuroendocrine cancers, including small cell lung cancer (SCLC), frequently display tuft-like subsets, a feature linked to poor patient outcomes. The mechanisms driving neuroendocrine–tuft tumour heterogeneity, and the origins of tuft-like cancers are unknown. Using multiple genetically-engineered animal models of SCLC, we demonstrate that a basal cell of origin (but not the accepted neuroendocrine origin) generates neuroendocrine–tuft-like tumours that highly recapitulate human SCLC. Single-cell clonal analyses of basal-derived SCLC further uncovers unexpected transcriptional states and lineage trajectories underlying neuroendocrine–tuft plasticity. Uniquely in basal cells, introduction of genetic alterations enriched in human tuft-like SCLC, including high MYC, PTEN loss, and ASCL1 suppression, cooperate to promote tuft-like tumours. Transcriptomics of 944 human SCLCs reveal a basal-like subset and a tuft-ionocyte-like state that altogether demonstrate remarkable conservation between cancer states and normal basal cell injury response mechanisms. Together, these data suggest that the basal cell is a plausible origin for SCLC and other neuroendocrine-tuft cancers that can explain neuroendocrine–tuft heterogeneity—offering new insights for targeting lineage plasticity.

## Manuscript
Manuscript pre-print is available online [here](https://www.biorxiv.org/content/10.1101/2024.11.13.623500v1).

## Data
* Data are accessible on GEO, under Superseries [GSE279200](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279200).  
* Raw data are provided as FASTQ and processed CellRanger count filtered barcodes, features, and matrices.  
* Processed scRNA-seq data have been uploaded to [Zenodo in record 15857303](https://zenodo.org/records/15857303)  

#### NOTE
Data to reproduce the figures can be downloaded automatically using [code/00_download_data.ipynb](code/00_download_data.ipynb)  

## Code
Code to replicate analyses performed in [Ireland et al. *Nature* 2025]

1. Code to reproduce the figures are provided in [code](code/) and preprocessing code is provided in [code/Preprocessing_for_reference]([code/Preprocessing_for_reference)  
1. Preprocessed data needed to reproduce the figures can automatically be downloaded using [code/00_download_data.ipynb](code/00_download_data.ipynb)  
1. Raw or processed data is also available from NCBI GEO [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279200).
1. See Jupyter Notebooks and R notebooks in [code](code) for manuscript specific dataset integration (preprocessing) methods and analyses. 
	* All samples were first subject to Scanpy/scVI-based QC and clustering in Python, and processed anndata objects were saved. 
	* R scripts call in the resulting anndata objects and convert them to Seurat objects for additional analyses. 
	* Many files called in R scripts are included in the [code/Preprocessing_for_reference/R_Code/metadata_files](code/Preprocessing_for_reference/R_Code/metadata_files).  
1. Please refer to [Scanpy](https://scanpy.readthedocs.io/en/stable/), [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/index.html), [Single-cell best practices](https://www.sc-best-practices.org/), [CellTagR](https://github.com/morris-lab/CellTagR) and [CellRank](https://cellrank.readthedocs.io/en/stable/notebooks/tutorials/general/100_getting_started.html) tutorials for system requirements and package-specific code that was adapted for datasets in this manuscript. 

#### Note 
To install necessary packages for executing Python code we recommend using the `scrnaseq312` conda environment that can be set up using the following (assuming you are in the top directory of this repository):

```
conda env create -f scrnaseq312.yaml
```

## Contact
Please consult methods described in our manuscript for more details or contact corresponding author for specific requests.

