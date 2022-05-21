# STATegra
Analysis of multi-omics dataset of B-cell differentiation in mouse

The recent breakthroughs in molecular technologies allow high throughput genomics, transcriptomics, metabolomics, and proteomics data generation to address
complex biological questions. The current task focuses on two high throughput sequencing data: transcriptomics and metabolomics. Using a dataset generated 
by Gomez-Cabrero et al. 2019 (https://www.nature.com/articles/s41597-019-0202-7; PMID: 31672995), following analyses are performed:

Experimental design: The mouse B3 cell line models the pre-BI (or Hardy fraction C’) stage. Upon nuclear translocation of the Ikaros transcription factor these 
cells progress to the pre-BII (or Hardy fraction D) stage, where B cell progenitors undergo growth arrest and differentiation. The B3 cell line was
retrovirally transduced with a vector encoding an Ikaros-REt2 fusion protein, which allows control of nuclear levels of Ikaros upon exposure to the drug
Tamoxifen. In parallel, cells were transfected with an empty vector to serve as control for the Tamoxifen effect.
After drug treatment, cultures were harvested at 0 h, 2 h, 6 h, 12 h, 18 h and 24 hs and profiled by several omics technologies:
long messenger RNA-seq (mRNA-seq) and micro RNA-seq (miRNA-seq) to measure gene expression; reduced representation by bisulfite sequencing (RRBS) to
measure DNA methylation; DNase-seq to measure chromatin accessibility as DNaseI Hypersensitive Sites (DHS) and transcription factor footprints, 
shotgun proteomics and targeted metabolomics of primary carbon and amino-acid metabolism. Moreover, single-cell RNA-seq (scRNA-seq) data for the entire 
time-series, while bulk ATAC-seq (ATAC-seq) and single-cell ATAC-seq (scATAC-seq) were obtained in a later round of experiments for 0 h and 24 h-time points 
of Ikaros induction only (no control series were run for these datasets). The dataset is complemented by existing ChIP-seq data on the same system equivalent
to our 0 h and 24 h time points34. In total, 793 different samples across the different omics datasets define the STATegra data collection.

1. Data exploration, visualization, and preprocessing strategies: A descriptive summary of the transcriptomics and metabolomics data using plots to explore the data 
characteristics and detect potential quality issues. Preprocessing strategies are applied to generate high quality gene and metabolite matrices for downstream analyses. 

2. Data prediction strategies: Multi-omics integration techniques on the two data types, i.e. transcriptomics and metabolomics ia applied to train an effective model
to predict sample conditions (i.e. exposed vs control). The contributin of  each data type to the prediction is also explored.

3. Data interpretation strategies: Genes, metabolites, functional modules or pathways that associated with control and exposed experimental conditions are identified.
How do they change over time? 
