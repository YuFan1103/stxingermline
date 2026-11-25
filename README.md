# Single-cell RNA sequencing of Drosophila testis during aging
Analysis scripts for scRNA-seq of young (1wo) and aged (9wo) Drosophila testes, including cell clustering, gene expression analysis, subclustering of early germ cells and cyst cells, identification of GSCs.

# Custom packages and installation:
These scripts rely on packages developed for this study. Install the following packages using BiocManager::install() function or install.packages().
1. Package: Seurat R-package (v4.4.0); download: https://satijalab.org/seurat/   description: single-cell RNA tool kit. We used standard workflows in Seurat package for quality control, counts normalization, cell clustering, and marker gene expression analysis.
2. Package: SCpubr R-package (v2.0.2) download:(https://enblacar.github.io/SCpubr-book-v1/). description: We used this package to visualize UMAP and dot plots.
3. Package: Monocle3 download: (https://cole-trapnell lab.github.io/monocle3/). description: we used this package to perform pseudotime analysis.

# Scripts:
1. UMAP visualization of single cell clusters of young and aged Drosophila testes, related to Fig.4&6 and &Extended Data Fig.5&6.R
   description: single-cell processing, clustering, annotation, and marker gene expression.
2. Monocle3, related to Fig.4.R
   description: pseudotime analysis of GSC differentiation.
3. Identification of niche-retained and displaced GSCs, related to Fig.4.R
   description: Identification of GSCs via marker stg and zpg, and marker gene expression in spermatogonia, spermatocyte in those GSCs.

# Data source and Usage
Input data: GEO dataset under the accession number GSE299425 (released upon publication).
Data source: https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299425 (reviewer token: gzydqeemflcdfot).

# Overview
Step 1: load the file of T1_T60_filtered_feature_bc_matrix and integrate young (T1) and aged (T60) dataset into a combined T1.combined.
Step 2: Run Harmony to correct batch effects
Step 3: Perform unsupervised subclustering of cyst cells and earely germ cells.
Step 4: Analysis of marker gene expression.
Expected output file: T1.combined.rds and Harmony_markers.csv
