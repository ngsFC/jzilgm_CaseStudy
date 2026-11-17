# JZILGM case study: tumour-associated macrophages in pediatric high-grade glioma

Real-data application of **JZILGM** [github.com/FranMigl/JZILGM](https://github.com/FranMigl/JZILGM)) to single-cell RNA-seq data of macrophages from six pediatric high-grade glioma samples.

- **Data:** Single-cell Pediatric Cancer Atlas (ScPCA), project **SCPCP000001**, samples SCPCS000002, 014, 016,
  019, 020, 022 (10x Chromium v3; ScPCA processing workflow v0.9.3). Macrophages = ScPCA c nsensus cell-type
  annotation "macrophage".

The raw ScPCA data are not included; `R/01_download_data.R` downloads them into `scpca_data_hgg6/`.

## How to run

```bash
export SCPCA_EMAIL="you@example.org"    
Rscript R/01_download_data.R             # download raw data (ScPCA account required)
Rscript R/02_preprocess.R                # preprocessing step 
Rscript R/03_fit_single_patient.R        # JZILGM per patient
Rscript R/04_fit_joint.R                 # JZILGM joint networks
Rscript R/05_network_analysis.R          # Jaccard, communities, ARI, hubs, degree
Rscript R/06_enrichment_msigdb.R         # MSigDB enrichment
Rscript R/07_figures.R                   # all figures
Rscript R/08_supplementary_tables.R      # results/Supplementary_tables.xlsx
```

## Requirements

R ≥ 4.5 and the following packages (versions used in `sessionInfo.txt`):

- **Modelling:** JZILGM (GitHub, FranMigl/JZILGM), Matrix
- **Preprocessing:** ScPCAr, SingleCellExperiment, Seurat, org.Hs.eg.db, AnnotationDbi
- **Analysis:** igraph, RobustRankAggreg, clusterProfiler, msigdbr (MSigDB 2026.1.Hs), dplyr, tidyr
- **Figures and tables:** ggplot2, ggraph, ggrepel, patchwork, ggVennDiagram, openxlsx

```r
install.packages("pak")
pak::pak(c("FranMigl/JZILGM", "igraph", "RobustRankAggreg", "dplyr", "tidyr", "ggplot2", "ggraph",
           "ggrepel", "patchwork", "ggVennDiagram", "openxlsx", "msigdbr", "Seurat",
           "bioc::clusterProfiler", "bioc::org.Hs.eg.db", "bioc::SingleCellExperiment",
           "AlexsLemonade/ScPCAr"))
```

## License

Code: MIT (see `LICENSE`). The ScPCA data are subject to the ScPCA Portal terms of use.
