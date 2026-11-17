# Network comparison and summary statistics (single-patient + joint networks).

source("R/00_config.R")
suppressPackageStartupMessages(library(RobustRankAggreg))

adj  <- list()   # adj[[paste(id, n)]]
memb <- list()   # Louvain membership, same keys
for (n in cutoffs) for (id in all_ids) {
  k <- paste(id, n); adj[[k]] <- read_adj(id, n); memb[[k]] <- louvain(adj[[k]])
}

# --- 1. Edge Jaccard -------------------------------------------------------------------
for (n in cutoffs) {
  J <- outer(all_ids, all_ids, Vectorize(function(i, j) edge_jaccard(adj[[paste(i, n)]], adj[[paste(j, n)]])))
  dimnames(J) <- list(all_ids, all_ids)
  saveRDS(J, sprintf("results/jaccard_matrix_top%d.rds", n))
}

# --- 2. Communities --------------------------------------------------------------------
for (j in joint_ids) {
  m <- memb[[paste(j, 700)]]; d <- degree_of(adj[[paste(j, 700)]])
  tab <- data.frame(gene = names(m), community = as.integer(m), degree = as.integer(d[names(m)]))
  tab <- tab[order(tab$community, -tab$degree), ]
  write.csv(tab, sprintf("results/communities_%s_top700.csv", j), row.names = FALSE)
}
core_mac <- Reduce(intersect, lapply(joint_ids, function(j) { m <- memb[[paste(j, 700)]]; names(m)[m == 2] }))
core_trl <- Reduce(intersect, lapply(joint_ids, function(j) { m <- memb[[paste(j, 700)]]; names(m)[m == 1] }))

comm_all <- do.call(rbind, lapply(cutoffs, function(n) do.call(rbind, lapply(joint_ids, function(j) {
  m <- memb[[paste(j, n)]]; d <- degree_of(adj[[paste(j, n)]])
  data.frame(joint = j, genes_in_network = n, gene = names(m), community = as.integer(m),
             module = label_modules(m, core_mac, core_trl), degree = as.integer(d[names(m)]))
}))))
write.csv(comm_all, "results/communities_joint_networks_all_sizes.csv", row.names = FALSE)
print(table(comm_all$genes_in_network, comm_all$module))

# --- 3. Community stability (ARI) ---------------------------------------------------------
ari <- function(a, b) compare(a, b, method = "adjusted.rand")
by_size <- NULL; held <- NULL
for (n in cutoffs) {
  m <- lapply(all_ids, function(id) memb[[paste(id, n)]]); names(m) <- all_ids
  g <- names(m[[1]]); m <- lapply(m, function(x) x[g])
  A <- outer(all_ids, all_ids, Vectorize(function(i, j) ari(m[[i]], m[[j]]))); dimnames(A) <- list(all_ids, all_ids)
  saveRDS(A, sprintf("results/ari_matrix_top%d.rds", n))

  ss <- A[sample_ids, sample_ids][upper.tri(A[sample_ids, sample_ids])]
  jj <- A[joint_ids, joint_ids][upper.tri(A[joint_ids, joint_ids])]
  for (j in joint_ids) for (h in setdiff(sample_ids, members[[j]]))
    held <- rbind(held, data.frame(genes = n, joint = j, heldout = h,
                                   ARI_joint_vs_heldout = A[j, h],
                                   ARI_member_vs_heldout = mean(A[members[[j]], h])))
  hn <- held[held$genes == n, ]
  by_size <- rbind(by_size, data.frame(
    genes = n,
    n_comm_single = paste(sapply(m[sample_ids], function(x) sum(table(x) >= 10)), collapse = ","),
    n_comm_joint  = paste(sapply(m[joint_ids],  function(x) sum(table(x) >= 10)), collapse = ","),
    ARI_single_mean = mean(ss), ARI_single_min = min(ss), ARI_single_max = max(ss),
    ARI_joint_min = min(jj), ARI_joint_max = max(jj),
    heldout_joint_mean = mean(hn$ARI_joint_vs_heldout), heldout_member_mean = mean(hn$ARI_member_vs_heldout),
    heldout_joint_wins = sum(hn$ARI_joint_vs_heldout > hn$ARI_member_vs_heldout)))
}

g300 <- rownames(read_adj("joint1", 300))
across <- do.call(rbind, lapply(joint_ids, function(j) {
  m <- lapply(cutoffs, function(n) memb[[paste(j, n)]][g300]); names(m) <- cutoffs
  mac <- lapply(cutoffs, function(n) { x <- comm_all[comm_all$joint == j & comm_all$genes_in_network == n, ]
                                       intersect(x$gene[x$module == "Macrophage identity"], g300) }); names(mac) <- cutoffs
  macJ <- sapply(cutoffs, function(n) { x <- comm_all[comm_all$joint == j & comm_all$genes_in_network == n, ]
                                        jaccard(x$gene[x$module == "Macrophage identity"], core_mac) })
  data.frame(joint = j, ARI_300v500 = ari(m[["300"]], m[["500"]]), ARI_500v700 = ari(m[["500"]], m[["700"]]),
             ARI_300v700 = ari(m[["300"]], m[["700"]]),
             mac_jaccard_core_300 = macJ[1], mac_jaccard_core_500 = macJ[2], mac_jaccard_core_700 = macJ[3],
             mac300_retained_at_700 = length(intersect(mac[["300"]], mac[["700"]])) / length(mac[["300"]]))
}))
rnd <- function(d) { d[] <- lapply(d, function(x) if (is.numeric(x) && any(x %% 1 != 0)) round(x, 3) else x); d }
write.csv(rnd(by_size), "results/community_ari_by_size.csv", row.names = FALSE)
write.csv(rnd(held),    "results/community_ari_heldout.csv", row.names = FALSE)
write.csv(rnd(across),  "results/community_ari_across_sizes.csv", row.names = FALSE)
print(rnd(by_size), row.names = FALSE); print(rnd(across), row.names = FALSE)

# --- 4. Consensus hubs (RRA) and edge counts -------------------------------------------------
deg_joint <- lapply(cutoffs, function(n) lapply(setNames(joint_ids, joint_ids), function(j) degree_of(adj[[paste(j, n)]])))
names(deg_joint) <- cutoffs
hub_sets <- lapply(cutoffs, function(n) {
  lists <- lapply(deg_joint[[as.character(n)]], function(d) names(sort(d, decreasing = TRUE)))
  agg <- aggregateRanks(glist = lists, N = n, method = "RRA")
  head(agg$Name[order(agg$Score)], N_HUBS)
}); names(hub_sets) <- paste0("top", cutoffs)
saveRDS(hub_sets, "results/hub_sets.rds")

hub_union <- Reduce(union, hub_sets)
hub_tab <- data.frame(gene = hub_union)
for (n in cutoffs) {
  hub_tab[[paste0("mean_degree_top", n)]] <- sapply(hub_union, function(gn) {
    mean(sapply(deg_joint[[as.character(n)]], function(d) if (gn %in% names(d)) d[[gn]] else NA_real_), na.rm = TRUE) })
  hub_tab[[paste0("is_hub_top", n)]] <- hub_union %in% hub_sets[[paste0("top", n)]]
}
hub_tab <- hub_tab[order(-hub_tab$mean_degree_top700), ]
write.csv(hub_tab, "results/hub_degree_by_cutoff.csv", row.names = FALSE)
cat(sprintf("Hubs: %d distinct genes, %d hubs at all sizes\n", nrow(hub_tab),
            sum(hub_tab$is_hub_top300 & hub_tab$is_hub_top500 & hub_tab$is_hub_top700)))

for (j in joint_ids) write.csv(data.frame(cutoff = cutoffs, n_edges = sapply(cutoffs, function(n) {
  a <- adj[[paste(j, n)]]; sum(a[upper.tri(a)] > 0) })), sprintf("results/edge_counts_%s.csv", j), row.names = FALSE)

# --- 5. Degree vs expression (300 genes) ------------------------------------------------------
N <- 300
deg <- sapply(all_ids, function(id) degree_of(adj[[paste(id, N)]]))
X <- readRDS(sprintf("data/macrophage_raw_top%d_list_dense_6pts.rds", N))
mean_expr <- rowMeans(sapply(X, function(m) colMeans(lognorm(m))))[rownames(deg)]
rk  <- apply(-deg, 2, rank, ties.method = "min")
top <- names(sort(rowMeans(deg[, joint_ids]), decreasing = TRUE))[1:15]
dex <- data.frame(gene = rownames(deg), mean_lognorm_expr = round(mean_expr, 3),
                  expression_rank = rank(-mean_expr, ties.method = "min"), deg, check.names = FALSE)
names(dex)[names(dex) %in% all_ids] <- paste0("degree_", ifelse(grepl("^joint", names(dex)[names(dex) %in% all_ids]),
                                              names(dex)[names(dex) %in% all_ids], short(names(dex)[names(dex) %in% all_ids])))
dex$degree_rank_min <- apply(rk, 1, min); dex$degree_rank_max <- apply(rk, 1, max)
dex$in_top15_figure <- ifelse(dex$gene %in% top, "yes", "")
dex <- dex[order(-rowMeans(deg[dex$gene, joint_ids])), ]
write.csv(dex, sprintf("results/degree_expression_top%d.csv", N), row.names = FALSE)
print(sapply(joint_ids, function(j) cor(deg[, j], mean_expr, method = "spearman")))
