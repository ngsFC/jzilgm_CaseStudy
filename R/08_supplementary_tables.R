# Supplementary Table

source("R/00_config.R")
suppressPackageStartupMessages(library(openxlsx))
LAB <- c(`1` = "Translation/OXPHOS", `2` = "Macrophage identity")

# --- genes, expression ------------------------------------------------------
X <- readRDS("data/macrophage_raw_top700_list_dense_6pts.rds")
genes <- colnames(X[[1]])                       # ordered by mean expression rank
mean_umi <- rowMeans(sapply(X, function(m) colMeans(m)[genes]))
pct_det  <- rowMeans(sapply(X, function(m) colMeans(m > 0)[genes])) * 100

ca_all <- read.csv("results/communities_joint_networks_all_sizes.csv")
comm_all <- list(); deg_all <- list()
for (j in joint_ids) for (n in cutoffs) {
  x <- ca_all[ca_all$joint == j & ca_all$genes_in_network == n, ]
  comm_all[[paste0(j, "_", n)]] <- x[, c("gene", "community", "module")]
  deg_all[[paste0(j, "_", n)]]  <- setNames(x$degree, x$gene)
}

# --- hubs --------------------------------------------------------------------
hub <- read.csv("results/hub_degree_by_cutoff.csv")

# --- Sheet: Genes (700) -------------------------------------------------------
G <- data.frame(gene = genes, expression_rank = seq_along(genes),
                mean_UMI_per_cell = round(mean_umi, 2), pct_cells_detected = round(pct_det, 1))
for (j in joint_ids) {
  ca <- comm_all[[paste0(j, "_700")]]; ca <- ca[match(genes, ca$gene), ]
  G[[paste0("community_", j)]] <- ca$community
  G[[paste0("module_", j)]]    <- ca$module
}
mods <- G[, paste0("module_", joint_ids)]
G$consensus_module <- apply(mods, 1, function(x) if (length(unique(x)) == 1) x[1] else "Mixed")
for (j in joint_ids) G[[paste0("degree_", j)]] <- unname(deg_all[[paste0(j, "_700")]][genes])
G$mean_degree_700 <- round(rowMeans(G[, paste0("degree_", joint_ids)]), 1)
for (n in cutoffs) G[[paste0("hub_top", n)]] <- ifelse(genes %in% hub$gene[hub[[paste0("is_hub_top", n)]]], "yes", "")
G$hub_all_sizes <- ifelse(G$hub_top300 == "yes" & G$hub_top500 == "yes" & G$hub_top700 == "yes", "yes", "")
G <- G[order(-G$mean_degree_700), ]

# --- Sheet: Hubs --------------------------------------------------------------
H <- hub
names(H) <- sub("mean_degree_top", "mean_degree_", sub("is_hub_top", "hub_", names(H)))
for (n in cutoffs) H[[paste0("hub_", n)]] <- ifelse(H[[paste0("hub_", n)]], "yes", "")
H$hub_all_sizes <- ifelse(H$hub_300 == "yes" & H$hub_500 == "yes" & H$hub_700 == "yes", "yes", "")
H <- merge(H, G[, c("gene", "expression_rank", "consensus_module",
                    paste0("degree_", joint_ids))], by = "gene", all.x = TRUE)
H$max_possible_degree_700 <- 699
H <- H[order(-H$mean_degree_700), c("gene", "consensus_module", "expression_rank",
        paste0("degree_", joint_ids), "mean_degree_700", "max_possible_degree_700",
        "mean_degree_300", "mean_degree_500", "hub_300", "hub_500", "hub_700", "hub_all_sizes")]
H[, c("mean_degree_300", "mean_degree_500", "mean_degree_700")] <-
  round(H[, c("mean_degree_300", "mean_degree_500", "mean_degree_700")], 1)

# --- Sheet: Communities (all sizes, summary) ----------------------------------
enr_full <- read.csv("results/msigdb_enrichment_full_top700.csv")
sig <- enr_full[enr_full$significant, ]
C <- do.call(rbind, lapply(names(comm_all), function(key) {
  j <- sub("_.*", "", key); n <- as.integer(sub(".*_", "", key))
  ca <- comm_all[[key]]; d <- deg_all[[key]]
  hubs_n <- hub$gene[hub[[paste0("is_hub_top", n)]]]
  ca$grp <- ifelse(ca$module == "Unconnected/singleton", "singletons", ca$community)
  do.call(rbind, lapply(split(ca, ca$grp), function(x) {
    cm <- if (x$module[1] == "Unconnected/singleton") "-" else as.character(x$community[1])
    s <- if (n == 700) sig[sig$joint == j & sig$community == cm, ] else sig[0, ]
    s <- s[order(s$p.adjust), ]
    data.frame(joint = j, patients = paste(sub("SCPCS000", "", unlist(members[[j]])), collapse = "+"),
      genes_in_network = n, community = cm, module = x$module[1], n_genes = nrow(x),
      n_hubs = sum(x$gene %in% hubs_n),
      hub_genes = paste(intersect(names(sort(d, decreasing = TRUE)), intersect(x$gene, hubs_n)), collapse = ", "),
      top10_genes_by_degree = paste(head(names(sort(d[x$gene], decreasing = TRUE)), 10), collapse = ", "),
      n_sig_Hallmark = if (n == 700 && cm != "-") sum(s$collection == "Hallmark") else NA,
      n_sig_C2_CP = if (n == 700 && cm != "-") sum(s$collection == "C2_CP") else NA,
      n_sig_C8 = if (n == 700 && cm != "-") sum(s$collection == "C8_celltype") else NA,
      top_sig_gene_sets_FDR = if (n == 700 && nrow(s)) paste(sprintf("%s (%.2g)", head(s$gene_set, 5), head(s$p.adjust, 5)), collapse = "; ")
                              else if (cm == "-") "not tested (no edges)" else if (n == 700) "none" else "not tested")
  }))
}))
C <- C[order(C$genes_in_network, C$joint, -C$n_genes), ]

# --- Sheet: Membership all sizes -----------------------------------------------
M <- data.frame(gene = genes)
for (n in cutoffs) for (j in joint_ids) {
  ca <- comm_all[[paste0(j, "_", n)]]
  M[[paste0(j, "_", n)]] <- ca$module[match(genes, ca$gene)]
}
M[is.na(M)] <- "not in panel"
M <- cbind(M[, 1, drop = FALSE], expression_rank = seq_along(genes), M[, -1])

# --- Sheets: enrichment ----------------------------------------------------------
E_short <- read.csv("results/msigdb_enrichment_short_top700.csv")
E_short$module <- LAB[as.character(E_short$community)]
E_short <- E_short[, c("collection", "community", "module", setdiff(names(E_short), c("collection", "community", "module")))]

# --- README ------------------------------------------------------------------------
readme <- data.frame(Sheet = c("README", "Cell_gene_counts", "Degree_expression_300", "Hubs", "Communities", "Genes_700", "Membership_all_sizes",
                               "Enrichment_summary", "Enrichment_significant", "Enrichment_all_tested"),
  Content = c(
    "Summary of the joint JZILGM networks of tumour-associated macrophages (6 pediatric HGG patients, ScPCA SCPCP000001). Joint networks: joint1 = SCPCS000014+019+022; joint2 = 016+020+022; joint3 = 002+014+020.",
    "Cell and gene counts per patient at each filtering step. Cells: unfiltered droplets -> emptyDropsCellRanger -> miQC (prob. compromised > 0.75 removed) and >= 200 detected genes (ScPCA processed) -> ScPCA consensus 'macrophage' annotation. Genes: reference -> detected in >= 1 macrophage -> after removing genes without symbol, MALAT1, MT-, RPS/RPL, haemoglobin and X/Y-chromosome genes -> 300/500/700 panel (shared by all patients). M1/M2: each macrophage scored with Seurat AddModuleScore on genes upregulated in each state (Mantovani et al. 2002; M1: CD80, CD86, TNF, IL1B; M2: CD163, MRC1, MSR1, TGFB1); M1-dominant if M1 score > M2 score, M2-dominant otherwise. Composition differs between patients (chi-squared test on the 6 x 2 table, given below the table).",
    "Data behind the hub-degree figure (top 300 genes): for every gene, mean log-normalised expression (log1p of counts per 10,000 over the 300 genes, averaged over cells and the 6 patients), expression rank, degree in each single-patient (P002 = SCPCS000002, ...) and joint network (max 299), min/max degree rank across the 9 networks, and whether the gene is among the 15 shown in the figure. Spearman correlation between degree and expression per network is given below the table.",
    "Consensus hub genes: per gene-set size, genes ranked by degree in each joint network, rankings combined with Robust Rank Aggregation (RRA); top 30 = hubs. Lists all 33 genes that are hubs at any size, with degree in each joint network (700 genes) and mean degree at 300/500/700 genes.",
    "One row per Louvain community (igraph cluster_louvain, seed 123) for every joint network and gene-set size, including small communities and singletons. Module labels assigned by Jaccard overlap with the two consensus modules of the 700-gene networks. Enrichment counts refer to the 700-gene networks (FDR < 0.05).",
    "All 700 genes: expression rank, mean UMI per cell and % cells detected (averaged over 6 patients), community and module in each joint network, degree, hub status at each size.",
    "Module assignment of every gene in each joint network at 300, 500 and 700 genes ('not in panel' = gene not among the top 300/500).",
    "Significant MSigDB gene sets aggregated over the three joint networks (700 genes): joints where significant, overlap, fold enrichment, FDR range, contributing genes.",
    "All significant gene sets (FDR < 0.05) per joint network and community, with contributing genes and adjusted p-values.",
    "All gene sets tested. Method: one-sided hypergeometric test (clusterProfiler::enricher); MSigDB 2026.1.Hs Hallmark, C2 canonical pathways (CP*), C8; universe = 700 analysed genes; gene-set size 5-500; Benjamini-Hochberg within each collection x community; significant = p.adjust < 0.05."))

cg <- read.csv("results/cell_gene_counts_per_patient.csv", check.names = FALSE)
cg <- cg[, c("sample", "age", "cells_unfiltered_droplets", "cells_after_emptyDrops",
             "cells_after_QC_miQC_minGenes200", "macrophages_consensus", "median_UMI_per_macrophage",
             "median_genes_per_macrophage", "genes_reference", "genes_detected_in_macrophages",
             "genes_detected_after_gene_filter")]
pol <- readRDS("data/polarization_df.rds")
pol_n <- t(sapply(cg$sample, function(s) { d <- pol$diff[pol$sample_id == s]; c(sum(d > 0), sum(d < 0)) }))
cg$macrophages_M1_dominant     <- pol_n[, 1]
cg$macrophages_M2_dominant     <- pol_n[, 2]
cg$pct_M1_dominant             <- round(100 * pol_n[, 1] / rowSums(pol_n), 1)
cg$pct_M2_dominant             <- round(100 * pol_n[, 2] / rowSums(pol_n), 1)
chi <- readRDS("results/macrophage_M1M2_composition_chisq.rds")
tot <- cg[1, ]; tot[] <- NA; tot$sample <- "Total"
for (k in c("cells_unfiltered_droplets", "cells_after_emptyDrops", "cells_after_QC_miQC_minGenes200",
            "macrophages_consensus", "macrophages_M1_dominant", "macrophages_M2_dominant")) tot[[k]] <- sum(cg[[k]])
tot$pct_M1_dominant <- round(100 * tot$macrophages_M1_dominant / (tot$macrophages_M1_dominant + tot$macrophages_M2_dominant), 1)
tot$pct_M2_dominant <- round(100 - tot$pct_M1_dominant, 1)
cg <- rbind(cg, tot)

dexp <- read.csv("results/degree_expression_top300.csv", check.names = FALSE)
dexp$in_top15_figure <- NULL
degm <- dexp[, grep("^degree_(P|joint)", names(dexp))]
readme <- readme[readme$Sheet %in% c("README", "Cell_gene_counts", "Degree_expression_300", "Hubs", "Communities",
                                     "Genes_700", "Membership_all_sizes", "Enrichment_summary"), ]

# --- write -------------------------------------------------------------------------------
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom", wrapText = TRUE, valign = "top")
add <- function(name, df, widths = NULL) {
  df[] <- lapply(df, function(x) if (is.character(x)) replace(x, x == "", NA) else x)
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs, withFilter = name != "README")
  freezePane(wb, name, firstRow = TRUE, firstCol = name != "README")
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = if (is.null(widths)) "auto" else widths)
}
add("README", readme, c(24, 140))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(readme) + 1), cols = 2)
add("Cell_gene_counts", cg)
writeData(wb, "Cell_gene_counts", data.frame(test = "Chi-squared, M1/M2 composition across patients",
          X_squared = round(unname(chi$statistic), 2), df = unname(chi$parameter),
          p_value = signif(chi$p.value, 3)), startRow = nrow(cg) + 3, headerStyle = hs)
add("Degree_expression_300", dexp)
sp <- data.frame(network = colnames(degm), spearman_rho_degree_vs_expression =
                 round(sapply(colnames(degm), function(k) cor(degm[[k]], dexp$mean_lognorm_expr, method = "spearman")), 3))
writeData(wb, "Degree_expression_300", sp, startRow = nrow(dexp) + 3, headerStyle = hs)
add("Hubs", H)
add("Communities", C, c(8, 10, 9, 10, 22, 8, 7, 60, 60, 9, 9, 9, 90))
add("Membership_all_sizes", M)
add("Genes_700", G)
add("Enrichment_summary", E_short, c(12, 10, 20, 50, 10, 8, 14, 12, 10, 10, 10, 90))

# highlight module colours
for (sh in c("Genes_700", "Membership_all_sizes", "Communities", "Hubs")) {
  df <- switch(sh, Genes_700 = G, Membership_all_sizes = M, Communities = C, Hubs = H)
  for (col in which(sapply(df, is.character))) {
    for (lv in c("Macrophage identity", "Translation/OXPHOS")) {
      r <- which(df[[col]] == lv)
      if (length(r)) addStyle(wb, sh, createStyle(fgFill = if (lv == "Macrophage identity") "#FBE3CC" else "#D6E4F2"),
                              rows = r + 1, cols = col, stack = TRUE)
    }
  }
}
out <- Sys.getenv("SUPP_OUT", "results/Supplementary_tables.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
