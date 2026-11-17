# Single-patient JZILGM networks (zinb_LGM_grp on one patient), top 300/500/700 genes.

library(JZILGM)
library(Matrix)
library(igraph)

source("R/00_config.R")
FIT_SEED <- SEED

available <- parallel::detectCores(logical = FALSE)
ncores <- available - 1


for (n_top in cutoffs) {

  in_path <- paste0("data/macrophage_raw_top", n_top, "_list_dense_6pts.rds")
  ListMat <- readRDS(in_path)
  stopifnot(length(ListMat) == 6)
  stopifnot(all(sapply(ListMat, ncol) == n_top))

  for (sid in sample_ids) {

    mat <- ListMat[[sid]]
    gene_names <- colnames(mat)

    single_list <- list(mat)
    names(single_list) <- sid

    running_time <- system.time({
      result <- zinb_LGM_grp(single_list, sym = "AND", nCores = ncores,
                              seed = FIT_SEED, verbose = 1)
    })[3]
    rt_h <- running_time / 3600

    out_path  <- paste0("data/res_real_", sid, "_top", n_top, ".RDS")
    time_path <- paste0("data/time_", sid, "_top", n_top, ".RDS")
    saveRDS(result, file = out_path)
    saveRDS(rt_h, file = time_path)

    adj <- as.matrix(result$quantile_network)
    rownames(adj) <- colnames(adj) <- gene_names
    saveRDS(adj, file = paste0("data/adj_", sid, "_top", n_top, ".rds"))

    g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    png(paste0("results/networks/net_", sid, "_top", n_top, ".png"), width = 1600, height = 1600, res = 200)
    set.seed(FIT_SEED)
    plot(g, layout = layout_with_fr(g), vertex.size = 3, vertex.label = NA,
         edge.width = 0.6, main = paste0(sid, " - top ", n_top, " genes"))
    dev.off()
  }
}
