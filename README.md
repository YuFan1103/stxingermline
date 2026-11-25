# Single-cell RNA sequencing of Drosophila testis during aging
Analysis scripts for scRNA-seq of young (1wo) and aged (9wo) Drosophila testes, including cell clustering, gene expression analysis, subclustering of early germ cells and cyst cells, identification of GSCs.

# Custom packages and installation:
These scripts rely on following packages installed in this study. Install the following packages using BiocManager::install() function or install.packages().
1. Seurat (v4.4.0)
   Source: https://satijalab.org/seurat/
   Installation: install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.4.0.tar.gz", repos = NULL, type = "source")
   Purpose: single-cell RNA-seq toolkit. We used the standard Seurat workflow for quality control, count normalization, cell clustering and marker-gene expression analysis.

2. SCpubr (v2.0.2)
   Source: https://github.com/enblacar/SCpubr (documentation: https://enblacar.github.io/SCpubr-book-v1/)
   Installation: remotes::install_github("enblacar/SCpubr@v2.0.2")
   Purpose: visualization of UMAP embeddings and dot plots.

3. Monocle3 (v1.3.1)
   Source: https://cole-trapnell-lab.github.io/monocle3/
   Installation: remotes::install_github("cole-trapnell-lab/monocle3@v1.3.1")
   Note: Monocle3 is distributed through GitHub. Its Bioconductor dependencies must be installed beforehand.
   Purpose: pseudotime/trajectory analysis of GSC differentiation.

## Installation script
Run the following commands in a fresh R session, in this order:

install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("batchelor", "BiocGenerics"))
install.packages("remotes")
remotes::install_github("enblacar/SCpubr@v2.0.2")
remotes::install_github("cole-trapnell-lab/monocle3@v1.3.1")
install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.4.0.tar.gz", repos = NULL, type = "source")

# Scripts and execution order

The three scripts should be run in numerical order. Each script reads the output of
the previous step and writes to the `output/` directory.

1. `UMAP visualization of single cell clusters of young and aged Drosophila testes.R`
   Description: single-cell processing, clustering, annotation, and marker gene expression.
   Input: 10x Genomics filtered feature-barcode matrices of T1 and T60.
   Output: `T1.combined.rds`, `Harmony_markers.csv`.

2. `Monocle3, related to Fig.4.R`
   Description: pseudotime analysis of GSC differentiation.
   Input: `T1.combined.rds`.
   Output: Monocle3 trajectory and pseudotime plots.

3. `Identification of niche-retained and displaced GSCs, related to Fig.4.R`
   Description: identification of GSCs via the markers stg and zpg, and marker-gene
   expression in spermatogonia and spermatocytes within these GSCs.
   Input: `T1.combined.rds` and the annotation from script 1.
   Output: GSC classification table and summary figures.

# Data source and usage

Input data: NCBI GEO accession GSE299425.
Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299425
The dataset is publicly released upon publication. Access for editors and reviewers before publication is provided through the private reviewer link and token issued
with the submission; the token is available in main manuscript of this study.

# Overview

All scripts were executed in RStudio on a [Windows 11] workstation. Full pipeline (scripts 1–3) wall time: approximately [6] h; script 1 ≈ [4] h, script 2 ≈ [1] h, script 3 ≈ [30] min.

Step 1: Load the filtered feature-barcode matrices of T1 and T60 and combine the
young (T1) and aged (T60) datasets into a single object, `T1.combined`.

Step 2: Run Harmony to correct batch effects between the two age groups.

Step 3: Perform unsupervised subclustering of cyst cells and early germ cells.

Step 4: Analyse marker-gene expression, including the GSC markers stg and zpg.

Expected output files: `T1.combined.rds` and `Harmony_markers.csv`.

