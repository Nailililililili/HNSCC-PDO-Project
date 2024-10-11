Genomic and single-cell characterization of patient-derived tumor organoid models of head and neck squamous cell carcinoma


Introduction
This project contains code, data, and output results for data analysis and visualization. Below is a detailed description of the folder structure.

1. code Folder

--01_TCGA_subtype_analysis.R          # Generates Fig4. This script determines the molecular subtype enrichment score of each single cell using the published genesets for each molecular subtype.

--02_gene_modlue_identification.R     # Generates Fig5. This script employs non-negative matrix factorization (NMF) on each of the 8 HNSCC PDOs to capture tumor ITH.

--03_RHP_scoring.R                    # Generates Fig6 and Fig7. This script scores the gene modlues in each single cell using a random gene set approach.

--04_HPV_HNSCC_analysis.R             # Generates Fig8. This script uses datasets from TCGA, CCLE and the Human Protein Atlas for HPV-negative samples to compare the correlation of the gene AREG with other genes.

--control_geneset.R                   # Functions used for gene modlue identification

--custom_magma.R                      # Functions used for gene modlue identification

--nmf_cell_class.R                    # Functions used for gene modlue identification

--nmf_programs.R                      # Functions used for gene modlue identification

--robust_nmf_programs.R               # Functions used for gene modlue identification

--seurat_functions_public.R           # Functions used for gene modlue scoring

control_geneset.R, custom_magma.R, nmf_cell_class.R, nmf_programs.R: These functions are used for gene modlue identification and are sourced from relevant literature: DOI: 10.1038/s41588-020-00726-6

seurat_functions_public.R: This function is used for gene modlue scoring and are sourced from relevant literature: DOI: 10.1038/s41588-022-01141-9


2. data Folder
   
This folder contains all the data files used in the project. 
The files srat_harmony_dims50_res0.5.RDS and HNSCC_expr.RDS in the data folder can be found on Zenodo：https://zenodo.org/records/13917309

--srat_harmony_dims50_res0.5.RDS      # merged seurat objects of 8 HNSCC PDOs 

--HNSCC_expr.RDS                      # the matrix of single cell data after log-normalization

--CCLE_heterogeneity_Rfiles           # HNSCC primary tumors from: DOI: 10.1038/s41588-020-00726-6

--CSCC_MPlist.rds                     # CSCC gene modules from: 

--epithelial_program.RDS              # ESCC gene modules from: DOI: 10.1038/s41467-021-25539-x

--MPgenelist.RDS                      # pan-cancer gene modules from: DOI: 10.1038/s41588-023-01570-0

--RHP_reference.RDS                   # HNSCC gene modules from: DOI: 10.1038/s41467-023-36691-x

--TCGA_HPV.RDS                        # published genesets for each molecular subtype

--bulk_expr.csv                       # bulk RNA-seq expression matrix of 13 HNSCC PDOs

--HNSCC_samples.all.txt               # metadata of 13 HNSCC PDOs bulk RNA-seq samples

--GSE9349_series_matrix.txt           # published data from GEO: GSE9349

--GPL201-30390.txt                    # Probe information for the GSE9349_series_matrix.txt file

--HPA.areg.correlation.csv            # datasets from TCGA for HPV-negative samples to compare the correlation of the gene AREG with other genes

--tcga.areg.correlation.csv           # datasets from CCLE and the Human Protein Atlas for HPV-negative samples to compare the correlation of the gene AREG with other genes


3. output Folder		
This folder is used to store the output results and related files for each figure. All generated figures and reports will be saved here.

 
