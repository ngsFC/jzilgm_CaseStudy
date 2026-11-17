# Preprocessing of the 6-patient pediatric HGG cohort (downloaded by 01_download_data.R).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(Matrix)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(Seurat)
})

source("R/00_config.R")

# ScPCA sample metadata (cell counts after each ScPCA filtering step)
meta <- read.delim(file.path(data_dir, "single-cell_metadata.tsv"))

# --- Sex-chromosome genes (excluded) --------------------------------------
sce_ref <- readRDS(file.path(data_dir, sample_ids[1],
                              paste0(library_ids[[sample_ids[1]]], "_processed.rds")))
all_ensembl_ids <- rowData(sce_ref)$gene_ids

gene_chr <- mapIds(org.Hs.eg.db, keys = all_ensembl_ids,
                    column = "CHR", keytype = "ENSEMBL", multiVals = "first")

xy_ensembl_ids <- all_ensembl_ids[gene_chr %in% c("X", "Y")]
rm(sce_ref)

# --- Load each sample and select macrophages --------------
macrophage_data <- list()
count_rows <- list()

for (sample_id in sample_ids) {

  lib_id <- library_ids[[sample_id]]
  sce_path <- file.path(data_dir, sample_id, paste0(lib_id, "_processed.rds"))
  sce <- readRDS(sce_path)

  is_macrophage <- colData(sce)$consensus_celltype_annotation == "macrophage"
  stopifnot(sum(is_macrophage) > 0)

  sce_mac <- sce[, is_macrophage]

  gene_symbol <- rowData(sce_mac)$gene_symbol
  gene_ids    <- rowData(sce_mac)$gene_ids

  keep_gene <- !is.na(gene_symbol) &
    !grepl("MALAT1", gene_symbol) &
    !grepl("^MT-", gene_symbol) &
    !grepl("^RP[SL]", gene_symbol) &
    !grepl("^HB[^(P|E|S)]", gene_symbol) &
    !(gene_ids %in% xy_ensembl_ids)

  raw_mat  <- as.matrix(counts(sce_mac)[keep_gene, , drop = FALSE])
  norm_mat <- as.matrix(logcounts(sce_mac)[keep_gene, , drop = FALSE])
  rownames(raw_mat)  <- gene_symbol[keep_gene]
  rownames(norm_mat) <- gene_symbol[keep_gene]

  macrophage_data[[sample_id]] <- list(raw = raw_mat, norm = norm_mat)

  # cell and gene counts at each step
  m <- meta[meta$scpca_sample_id == sample_id, ]
  det <- Matrix::rowSums(counts(sce_mac) > 0) > 0
  count_rows[[sample_id]] <- data.frame(
    sample = sample_id, diagnosis = m$diagnosis, age = m$age,
    cells_unfiltered_droplets = m$unfiltered_cells,
    cells_after_emptyDrops = m$filtered_cell_count,
    cells_after_QC_miQC_minGenes200 = m$processed_cells,
    macrophages_consensus = ncol(sce_mac),
    pct_macrophages = round(100 * ncol(sce_mac) / m$processed_cells, 1),
    median_UMI_per_macrophage = median(colSums(counts(sce_mac))),
    median_genes_per_macrophage = median(colSums(counts(sce_mac) > 0)),
    genes_reference = nrow(sce),
    genes_detected_in_macrophages = sum(det),
    genes_detected_after_gene_filter = sum(det & keep_gene),
    macrophages_called_doublet = sum(sce_mac$scDblFinder_class == "doublet", na.rm = TRUE))
}

# --- M1 / M2 polarization scoring -----------------------------------------
# Markers restricted to the UPREGULATED side of each state as reported in:
#   Mantovani A, Sozzani S, Locati M, Allavena P, Sica A. Macrophage
#   polarization: tumor-associated macrophages as a paradigm for polarized
#   M2 mononuclear phagocytes. Trends Immunol. 2002;23(11):549-555.
#   DOI: 10.1016/S1471-4906(02)02302-5
m1_markers <- c("CD80", "CD86", "TNF", "IL1B")
m2_markers <- c("CD163", "MRC1", "MSR1", "TGFB1")

polarization_list <- list()

for (sid in sample_ids) {

  raw_mat  <- macrophage_data[[sid]]$raw
  norm_mat <- macrophage_data[[sid]]$norm
  dedup <- !duplicated(rownames(raw_mat))

  seu <- CreateSeuratObject(counts = raw_mat[dedup, , drop = FALSE])
  seu <- SetAssayData(seu, layer = "data", new.data = norm_mat[dedup, , drop = FALSE])

  m1_present <- intersect(m1_markers, rownames(seu))
  m2_present <- intersect(m2_markers, rownames(seu))

  seu <- AddModuleScore(
    seu,
    features = list(M1 = m1_present, M2 = m2_present),
    name = "Pol",
    seed = 123
  )

  polarization_list[[sid]] <- data.frame(
    sample_id = sid,
    cell_id   = colnames(norm_mat),
    m1_score  = seu@meta.data$Pol1,
    m2_score  = seu@meta.data$Pol2,
    diff      = seu@meta.data$Pol1 - seu@meta.data$Pol2
  )
}

polarization_df <- bind_rows(polarization_list)

polarization_summary <- polarization_df %>%
  group_by(sample_id) %>%
  summarise(
    n_cells          = n(),
    mean_m1_score    = mean(m1_score),
    mean_m2_score    = mean(m2_score),
    mean_diff        = mean(diff),
    pct_M1_dominant  = mean(diff > 0) * 100,
    pct_M2_dominant  = mean(diff < 0) * 100,
    .groups = "drop"
  )

print(polarization_summary)
write.csv(polarization_summary, "results/macrophage_polarization_summary.csv", row.names = FALSE)

saveRDS(polarization_df, "data/polarization_df.rds")

composition_table <- polarization_df %>%
  group_by(sample_id) %>%
  summarise(M1 = sum(diff > 0), M2 = sum(diff < 0), .groups = "drop")

composition_mat <- as.matrix(composition_table[, c("M1", "M2")])
rownames(composition_mat) <- composition_table$sample_id

composition_chisq <- chisq.test(composition_mat)
print(composition_chisq)

saveRDS(composition_chisq, "results/macrophage_M1M2_composition_chisq.rds")

# --- Gene ranking and top-N selection ---------------------------------------
gene_means <- lapply(macrophage_data, function(x) rowMeans(x$norm))
gene_ranks <- lapply(gene_means, function(gm) rank(-gm, ties.method = "average"))

common_genes <- Reduce(intersect, lapply(gene_ranks, names))

mean_ranks <- sapply(common_genes, function(gene) {
  mean(sapply(gene_ranks, function(gr) gr[gene]))
})

cutoffs <- c(300, 500, 700)

top_genes_list <- list()

for (n_top in cutoffs) {

  top_genes <- names(sort(mean_ranks))[1:n_top]
  top_genes_list[[paste0("top", n_top)]] <- top_genes

  topN_list <- lapply(macrophage_data, function(x) {
    mat_sub <- x$raw[top_genes, , drop = FALSE]
    t(mat_sub)
  })

  topN_list_dense <- lapply(topN_list, function(mat) {
    if (!inherits(mat, "matrix")) as.matrix(mat) else mat
  })

  out_path <- paste0("data/macrophage_raw_top", n_top, "_list_dense_6pts.rds")
  saveRDS(topN_list_dense, file = out_path)
}

# --- Cell and gene counts table ---------------------------------------------
counts_tab <- do.call(rbind, count_rows)
for (n_top in cutoffs) counts_tab[[paste0("genes_in_panel_", n_top)]] <- n_top
write.csv(counts_tab, "results/cell_gene_counts_per_patient.csv", row.names = FALSE)
