# Joint (3-patient) JZILGM networks, top 300/500/700 genes.

library(JZILGM)
library(Matrix)
library(igraph)

source("R/00_config.R")

FIT_SEED   <- SEED
COMBO_SEED <- 123
N_COMBOS   <- 3

available <- parallel::detectCores(logical = FALSE)
ncores <- available - 1


# --- Pick N_COMBOS of the 20 possible 3-of-6 triplets-------
all_triplets <- combn(sample_ids, 3)
stopifnot(ncol(all_triplets) == 20)

set.seed(COMBO_SEED)
chosen_idx <- sample(seq_len(ncol(all_triplets)), N_COMBOS)
chosen_triplets <- lapply(chosen_idx, function(i) all_triplets[, i])

# Manual override: the seed-selected joint3 shared 2 of its 3 patients with
# joint1 (SCPCS000014, SCPCS000019). Swap SCPCS000019 -> SCPCS000002 so each
# patient is reused at most twice across the 3 triplets instead of 3 times.
chosen_triplets[[3]][chosen_triplets[[3]] == "SCPCS000019"] <- "SCPCS000002"

triplet_df <- do.call(rbind, lapply(seq_along(chosen_triplets), function(k) {
  data.frame(joint_id = paste0("joint", k),
             patient_1 = chosen_triplets[[k]][1],
             patient_2 = chosen_triplets[[k]][2],
             patient_3 = chosen_triplets[[k]][3])
}))

dir.create("results", showWarnings = FALSE)
write.csv(triplet_df, "results/joint_triplets.csv", row.names = FALSE)
print(triplet_df)

for (n_top in cutoffs) {

  in_path <- paste0("data/macrophage_raw_top", n_top, "_list_dense_6pts.rds")
  ListMat <- readRDS(in_path)
  stopifnot(length(ListMat) == 6)

  for (k in seq_along(chosen_triplets)) {

    trio <- chosen_triplets[[k]]
    joint_id <- paste0("joint", k)

    trio_list <- ListMat[trio]
    gene_names <- colnames(trio_list[[1]])

    running_time <- system.time({
      result <- zinb_LGM_grp(trio_list, sym = "AND", nCores = ncores,
                              seed = FIT_SEED, verbose = 1)
    })[3]
    rt_h <- running_time / 3600

    out_path  <- paste0("data/res_real_", joint_id, "_top", n_top, ".RDS")
    time_path <- paste0("data/time_", joint_id, "_top", n_top, ".RDS")
    saveRDS(result, file = out_path)
    saveRDS(rt_h, file = time_path)

    adj <- as.matrix(result$quantile_network)
    rownames(adj) <- colnames(adj) <- gene_names
    saveRDS(adj, file = paste0("data/adj_", joint_id, "_top", n_top, ".rds"))

    g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    png(paste0("results/networks/net_", joint_id, "_top", n_top, ".png"), width = 1600, height = 1600, res = 200)
    set.seed(FIT_SEED)
    plot(g, layout = layout_with_fr(g), vertex.size = 3, vertex.label = NA,
         edge.width = 0.6, main = paste0(joint_id, " (", paste(trio, collapse = "+"), ") - top ", n_top))
    dev.off()
  }
}
