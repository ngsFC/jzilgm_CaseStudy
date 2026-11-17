suppressPackageStartupMessages(library(igraph))

proj_dir <- Sys.getenv("REALDATA_DIR", ".")
if (dir.exists(proj_dir)) setwd(proj_dir)
dir.create("data", showWarnings = FALSE)
dir.create("results/networks", showWarnings = FALSE, recursive = TRUE)

# --- samples, gene-set sizes, seeds ---------------------------------------------
sample_ids  <- c("SCPCS000002", "SCPCS000014", "SCPCS000016",
                 "SCPCS000019", "SCPCS000020", "SCPCS000022")
library_ids <- setNames(sub("SCPCS", "SCPCL", sample_ids), sample_ids)
data_dir    <- "scpca_data_hgg6/SCPCP000001_single-cell"
cutoffs     <- c(300, 500, 700)
joint_ids   <- c("joint1", "joint2", "joint3")
all_ids     <- c(sample_ids, joint_ids)
SEED        <- 123
N_HUBS      <- 30 
PADJ_CUTOFF <- 0.05

trip_file <- "results/joint_triplets.csv"
trip    <- if (file.exists(trip_file)) read.csv(trip_file) else NULL
members <- if (!is.null(trip)) lapply(split(trip[, 2:4], trip$joint_id), unlist) else NULL

# --- helpers -----------------------------------------------------------------------
short    <- function(x) sub("SCPCS000", "P", x)                       
read_adj <- function(id, n) readRDS(sprintf("data/adj_%s_top%d.rds", id, n))
as_graph <- function(a) graph_from_adjacency_matrix(a, mode = "undirected", diag = FALSE)
degree_of <- function(a) { a <- as.matrix(a) > 0; diag(a) <- FALSE; rowSums(a) }
lognorm  <- function(m) log1p(m / rowSums(m) * 1e4)                   
jaccard  <- function(a, b) length(intersect(a, b)) / length(union(a, b))

louvain <- function(a) { set.seed(SEED); membership(cluster_louvain(as_graph(a))) }

edge_jaccard <- function(a1, a2) {
  e1 <- a1[upper.tri(a1)] > 0; e2 <- a2[upper.tri(a2)] > 0
  u <- sum(e1 | e2); if (u == 0) NA_real_ else sum(e1 & e2) / u
}

label_modules <- function(memb, core_mac, core_trl) {
  sz <- table(memb); comms <- names(sz)[sz > 1]
  lab <- setNames(rep("Unconnected/singleton", length(sz)), names(sz))
  jm <- sapply(comms, function(c) jaccard(names(memb)[memb == c], core_mac))
  lab[comms[which.max(jm)]] <- "Macrophage identity"
  rest <- setdiff(comms, comms[which.max(jm)])
  if (length(rest)) {
    jt <- sapply(rest, function(c) jaccard(names(memb)[memb == c], core_trl))
    lab[rest[which.max(jt)]] <- "Translation/OXPHOS"
    oth <- setdiff(rest, rest[which.max(jt)])
    lab[oth[order(-sz[oth])]] <- paste("Other", seq_along(oth))
  }
  unname(lab[as.character(memb)])
}

# --- figure colours ------------------------------------------------------------------
col_modules <- c("Translation/OXPHOS" = "#3b75af", "Macrophage identity" = "#e3812b")
col_extra   <- c("#4f9a6b", "#8e6bb8", "#c4577a", "#8c8c3a", "#5bb3c4")
col_heat    <- c(low = "#eef4fa", high = "#1a4f7a")
