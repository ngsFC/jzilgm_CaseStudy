# R scripts — real-data application (JZILGM, pediatric HGG macrophages)

Run from the repository root, in order (`Rscript R/05_network_analysis.R`), or set `REALDATA_DIR`.
Every script sources `00_config.R` (paths, sample IDs, gene-set sizes, seed, shared helpers).

| Script | Does | Main outputs |
|---|---|---|
| `00_config.R` | settings + helpers (Louvain, Jaccard, module labels, colours) | — |
| `01_download_data.R` | download the 6 ScPCA samples (needs ScPCA token) | `scpca_data_hgg6/` |
| `02_preprocess.R` | macrophages, gene filter, M1/M2 scores, top 300/500/700 genes, cell/gene counts | `data/macrophage_raw_top*_*.rds`, `data/polarization_df.rds`, `results/cell_gene_counts_per_patient.csv` |
| `03_fit_single_patient.R` | JZILGM per patient (slow) | `data/res_real_<patient>_*`, `data/adj_<patient>_*` |
| `04_fit_joint.R` | patient triplets + joint JZILGM (slow) | `results/joint_triplets.csv`, `data/res_real_joint*`, `data/adj_joint*` |
| `05_network_analysis.R` | edge Jaccard, Louvain communities, ARI, consensus hubs (RRA), degree vs expression | `results/*_matrix_top*.rds`, `communities_*`, `community_ari_*`, `hub_*`, `degree_expression_top300.csv` |
| `06_enrichment_msigdb.R` | MSigDB enrichment of 700-gene communities | `results/msigdb_enrichment_{full,short}_top700.csv` |
| `07_figures.R` | all figures | `results/macrophage_M1M2_composition.*`, `results/networks/*` |
| `08_supplementary_tables.R` | Supplementary Table S1 | `results/Supplementary_tables.xlsx` (overwrites; `SUPP_OUT` to redirect) |

