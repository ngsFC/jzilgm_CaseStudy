# All figures of the real-data application
#
#  1. results/macrophage_M1M2_composition.{png,pdf}      M1/M2 composition per patient
#  2. results/networks/jaccard_heatmap_panel.{png,pdf}    edge Jaccard, single + joint, 300/500/700
#  3. results/networks/ari_heatmap_panel.{png,pdf}        community ARI, single + joint, 300/500/700
#  4. results/networks/joint_networks_panel_3x3.{png,pdf} joint networks x gene-set size, communities
#  5. results/networks/hub_venn_cutoffs.png               consensus hubs across gene-set sizes
#  6. results/networks/hub_degree_robustness{,_top500,_top700,_panel}.{png,pdf}
#                                                        top-15 degrees + degree vs expression (300/500/700)
#  7. results/networks/edge_venn_joint{k}.png             edges shared across gene-set sizes
#  8. results/networks/community_network_joint{k}_top700.png  communities of each 700-gene joint network

source("R/00_config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(ggraph); library(patchwork); library(ggrepel)
  library(ggVennDiagram); library(dplyr); library(tidyr)
})

save_fig <- function(p, file, w, h, pdf = TRUE) {
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  if (pdf) ggsave(paste0(file, ".pdf"), p, width = w, height = h, bg = "white", device = cairo_pdf)
}

# --- 1. M1/M2 composition ----------------------------------------------------------------
pol <- readRDS("data/polarization_df.rds") %>%
  group_by(sample_id) %>%
  summarise(pct_M1 = mean(diff > 0) * 100, pct_M2 = mean(diff < 0) * 100, .groups = "drop") %>%
  pivot_longer(c(pct_M1, pct_M2), names_to = "state", values_to = "pct")
p1 <- ggplot(pol, aes(sample_id, pct, fill = state)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(pct_M1 = "#2980b9", pct_M2 = "#e74c3c"), labels = c("M1-dominant", "M2-dominant")) +
  labs(title = "M1/M2 cell-state composition per patient", x = "Sample", y = "% of macrophage cells", fill = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave("results/macrophage_M1M2_composition.png", p1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/macrophage_M1M2_composition.pdf", p1, width = 8, height = 6)

# --- 2-3. Jaccard / ARI heatmap panels ----------------------------------------------------
heat_panel <- function(M, n, limits, fill_name, white_above, show_y = TRUE) {
  ids <- rownames(M)
  lab <- ifelse(grepl("^joint", ids), ids, short(ids))
  grp <- factor(ifelse(grepl("^joint", ids), "Joint", "Single patient"), c("Single patient", "Joint"))
  d <- expand.grid(i = seq_along(ids), j = seq_along(ids))
  d$x <- factor(lab[d$j], lab); d$y <- factor(lab[d$i], rev(lab))
  d$gx <- grp[d$j]; d$gy <- grp[d$i]
  d$val <- M[cbind(d$i, d$j)]; d$val[d$i == d$j] <- NA
  d$member <- mapply(function(a, b) (grepl("^joint", a) && b %in% members[[a]]) ||
                                    (grepl("^joint", b) && a %in% members[[b]]), ids[d$i], ids[d$j])
  ggplot(d, aes(x, y)) +
    geom_tile(aes(fill = val), colour = "white", linewidth = 0.6) +
    geom_tile(data = d[d$member, ], fill = NA, colour = "black", linewidth = 0.7, width = 0.94, height = 0.94) +
    geom_text(aes(label = ifelse(is.na(val), "", sprintf("%.2f", val)), colour = !is.na(val) & val > white_above), size = 2.9) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_fill_gradient(low = col_heat[["low"]], high = col_heat[["high"]], limits = limits,
                        oob = scales::squish, na.value = "grey92", name = fill_name) +
    facet_grid(gy ~ gx, scales = "free", space = "free", switch = "y") +
    labs(title = sprintf("Top %d genes", n), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), panel.spacing = unit(7, "pt"),
          strip.text = element_text(face = "bold", size = 9, colour = "grey25"),
          strip.text.y.left = element_text(angle = 90), strip.placement = "outside",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
          axis.text.y = if (show_y) element_text(size = 8.5) else element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
}
heat_figure <- function(prefix, limits, fill_name, white_above) {
  ps <- lapply(seq_along(cutoffs), function(k) {
    n <- cutoffs[k]
    heat_panel(readRDS(sprintf("results/%s_matrix_top%d.rds", prefix, n)), n, limits, fill_name, white_above, show_y = k == 1)
  })
  wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect") &
    theme(legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))
}
save_fig(heat_figure("jaccard", c(0.25, 0.70), "Edge Jaccard", 0.5), "results/networks/jaccard_heatmap_panel", 16, 5.6)
save_fig(heat_figure("ari", c(0, 0.75), "Community ARI", 0.45), "results/networks/ari_heatmap_panel", 16, 5.6)

# --- 4. Joint networks, 3 x 3 panel ---------------------------------------------------------
memb700  <- lapply(joint_ids, function(j) { d <- read.csv(sprintf("results/communities_%s_top700.csv", j)); setNames(d$community, d$gene) })
core_mac <- Reduce(intersect, lapply(memb700, function(m) names(m)[m == 2]))
core_trl <- Reduce(intersect, lapply(memb700, function(m) names(m)[m == 1]))
col_main <- c(Translation = col_modules[["Translation/OXPHOS"]], Macrophage = col_modules[["Macrophage identity"]])

net_panel <- function(jid, n, n_labels = 6) {
  g <- as_graph(read_adj(jid, n))
  set.seed(SEED); V(g)$comm <- as.integer(membership(cluster_louvain(g)))
  V(g)$degree <- degree(g)
  n_isolated <- sum(V(g)$degree == 0)
  g <- induced_subgraph(g, which(V(g)$degree > 0))

  comms <- sort(unique(V(g)$comm)); sizes <- sapply(comms, function(c) sum(V(g)$comm == c))
  ov_mac <- sapply(comms, function(c) jaccard(V(g)$name[V(g)$comm == c], core_mac))
  ov_trl <- sapply(comms, function(c) jaccard(V(g)$name[V(g)$comm == c], core_trl))
  lab <- setNames(rep(NA_character_, length(comms)), comms)
  lab[as.character(comms[which.max(ov_mac)])] <- "Macrophage"
  rest <- comms[is.na(lab)]
  if (length(rest)) lab[as.character(rest[which.max(ov_trl[match(rest, comms)])])] <- "Translation"
  others <- comms[is.na(lab)]; others <- others[order(-sizes[match(others, comms)])]
  other_names <- if (length(others)) paste0("Other ", seq_along(others)) else character(0)
  lab[as.character(others)] <- other_names
  V(g)$group <- unname(lab[as.character(V(g)$comm)])
  pal <- c(col_main, setNames(col_extra[seq_along(others)], other_names))

  el <- as_edgelist(g, names = FALSE)
  same <- V(g)$comm[el[, 1]] == V(g)$comm[el[, 2]]
  E(g)$within <- same; E(g)$group <- ifelse(same, V(g)$group[el[, 1]], "between")
  set.seed(SEED); xy <- layout_with_fr(g, weights = ifelse(same, 1, 0.05), niter = 2000)
  top <- order(V(g)$degree, decreasing = TRUE)[seq_len(n_labels)]
  V(g)$label <- NA_character_; V(g)$label[top] <- V(g)$name[top]
  ne <- ecount(g); a_in <- min(0.35, 60 / ne * 4); a_bt <- a_in / 4
  mem <- members[[jid]]
  size_txt <- paste(sort(sizes[sizes > 1], decreasing = TRUE), collapse = " / ")

  ggraph(g, layout = "manual", x = xy[, 1], y = xy[, 2]) +
    geom_edge_link(aes(filter = !within), colour = "grey55", alpha = a_bt, width = 0.15) +
    geom_edge_link(aes(filter = within, colour = group), alpha = a_in, width = 0.2, show.legend = FALSE) +
    geom_node_point(aes(fill = group, size = degree), shape = 21, colour = "white", stroke = 0.25) +
    geom_node_text(aes(label = label), size = 2.6, fontface = "bold", colour = "grey15",
                   repel = TRUE, bg.colour = "white", bg.r = 0.15, min.segment.length = 0,
                   segment.colour = "grey40", segment.size = 0.25, max.overlaps = Inf) +
    scale_fill_manual(values = pal, breaks = names(pal), name = NULL) +
    scale_edge_colour_manual(values = pal, guide = "none") +
    scale_size_continuous(range = c(0.6, 5), guide = "none") +
    coord_equal(clip = "off") +
    labs(title = sprintf("%s  ·  top %d genes", jid, n),
         subtitle = sprintf("%s\n%s edges · communities: %s%s",
                            paste(sub("SCPCS000", "", mem), collapse = " + "), format(ne, big.mark = ","), size_txt,
                            ifelse(n_isolated > 0, sprintf(" · %d unconnected", n_isolated), ""))) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, colour = "grey35", hjust = 0.5, lineheight = 1.1),
          legend.position = "none", plot.margin = margin(6, 6, 6, 6))
}
plots <- list(); for (jid in joint_ids) for (n in cutoffs) plots[[paste(jid, n)]] <- net_panel(jid, n)

legend_plot <- ggplot(data.frame(x = 1:3, y = 1, group = c("Translation", "Macrophage", "Other")), aes(x, y, fill = group)) +
  geom_point(shape = 21, size = 4, colour = "white") +
  scale_fill_manual(values = c(col_main, Other = col_extra[1]), breaks = c("Translation", "Macrophage", "Other"),
                    labels = c("Translation / OXPHOS module (Hallmark MYC targets, OXPHOS)",
                               "Macrophage-identity module (MSigDB C8 myeloid)",
                               "Additional network-specific communities (green, purple, ...)"), name = NULL) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 5))) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10, margin = margin(r = 28)),
        legend.key.spacing.x = unit(4, "pt"))
lg <- ggplotGrob(legend_plot)
shared_legend <- wrap_elements(full = lg$grobs[[grep("guide-box", sapply(lg$grobs, function(x) x$name))[1]]])
p4 <- wrap_plots(plots, ncol = 3) / shared_legend + plot_layout(heights = c(1, 0.04)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
ggsave("results/networks/joint_networks_panel_3x3.png", p4, width = 15, height = 16, dpi = 300, bg = "white")
ggsave("results/networks/joint_networks_panel_3x3.pdf", p4, width = 15, height = 16, bg = "white")

# --- 5. Hub Venn diagram ----------------------------------------------------------------------
hub_sets <- readRDS("results/hub_sets.rds")
p5 <- ggVennDiagram(hub_sets, label = "count") +
  scale_fill_gradient(low = "#eaf4fd", high = "#1a5276") +
  labs(title = "Top-30 aggregated hub genes in common across gene-count cutoffs") +
  theme(legend.position = "none")
ggsave("results/networks/hub_venn_cutoffs.png", p5, width = 7, height = 6, dpi = 300)

# --- 6. Hub degree robustness (300 / 500 / 700 genes) --------------------------------------------
hub_degree_fig <- function(N, title = NULL) {
  dex <- read.csv(sprintf("results/degree_expression_top%d.csv", N), check.names = FALSE)
  net_cols <- c(paste0("degree_", short(sample_ids)), paste0("degree_", joint_ids))
  net_lab  <- c(short(sample_ids), joint_ids)
  top <- dex$gene[dex$in_top15_figure == "yes"]
  dA <- dex[dex$gene %in% top, c("gene", net_cols)] %>%
    pivot_longer(-gene, names_to = "net", values_to = "degree") %>%
    mutate(x = factor(net_lab[match(net, net_cols)], net_lab), y = factor(gene, rev(top)),
           frac = degree / (N - 1),
           grp = factor(ifelse(grepl("joint", net), "Joint", "Single patient"), c("Single patient", "Joint")))
  pA <- ggplot(dA, aes(x, y)) +
    geom_tile(aes(fill = frac), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = degree, colour = frac > 0.55), size = 2.8) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_fill_gradient(low = col_heat[["low"]], high = col_heat[["high"]], limits = c(0, 1), labels = scales::percent,
                        name = sprintf("Degree / max\n(%d)", N - 1)) +
    facet_grid(. ~ grp, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), panel.spacing = unit(7, "pt"),
          strip.text = element_text(face = "bold", size = 9, colour = "grey25"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
          axis.text.y = element_text(size = 8.5, face = "italic"),
          legend.position = "bottom", legend.title = element_text(size = 8.5),
          legend.key.width = unit(22, "pt"), legend.key.height = unit(8, "pt"))

  dB <- dex[, c("gene", "mean_lognorm_expr", paste0("degree_", joint_ids))] %>%
    pivot_longer(starts_with("degree_"), names_to = "joint", values_to = "degree") %>%
    mutate(joint = sub("degree_", "", joint))
  rho <- dB %>% group_by(joint) %>% summarise(rho = cor(degree, mean_lognorm_expr, method = "spearman"))
  dB <- dB %>% mutate(panel = sprintf("%s  (Spearman ρ = %.2f)", joint, rho$rho[match(joint, rho$joint)]),
                      hub = gene %in% top)
  pB <- ggplot(dB, aes(mean_lognorm_expr, degree)) +
    geom_point(data = dB[!dB$hub, ], colour = "grey65", size = 1, alpha = 0.7) +
    geom_point(data = dB[dB$hub, ], colour = "#c0392b", size = 1.6) +
    geom_text_repel(data = dB[dB$gene %in% top[1:8], ], aes(label = gene), size = 2.6, fontface = "italic",
                    min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf, seed = 1) +
    facet_wrap(~ panel, nrow = 1) +
    scale_y_continuous(limits = c(0, N - 1)) +
    labs(x = "Mean log-normalised expression (6 patients)", y = sprintf("Degree (top %d genes)", N)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(size = 9))
  if (!is.null(title)) pA <- pA + labs(title = title) + theme(plot.title = element_text(face = "bold", size = 12))
  pA + pB + plot_layout(widths = c(1, 2.3))
}

hub_figs <- lapply(setNames(cutoffs, cutoffs), hub_degree_fig)
save_fig(hub_figs[["300"]], "results/networks/hub_degree_robustness", 15, 5.6)   # original 300-gene figure
for (n in c(500, 700)) save_fig(hub_figs[[as.character(n)]], sprintf("results/networks/hub_degree_robustness_top%d", n), 15, 5.6)
hub_panel <- wrap_plots(lapply(cutoffs, function(n) hub_degree_fig(n, sprintf("%d genes", n))), ncol = 1)
save_fig(hub_panel, "results/networks/hub_degree_robustness_panel", 15, 17)


# --- 7. Edge Venn diagrams per joint network --------------------------------------------------
for (j in joint_ids) {
  edge_sets <- lapply(cutoffs, function(n) {
    a <- read_adj(j, n); g <- rownames(a); idx <- which(upper.tri(a) & a > 0, arr.ind = TRUE)
    paste(pmin(g[idx[, 1]], g[idx[, 2]]), pmax(g[idx[, 1]], g[idx[, 2]]), sep = "__")
  }); names(edge_sets) <- paste0("top", cutoffs)
  p7 <- ggVennDiagram(edge_sets, label = "count") +
    scale_fill_gradient(low = "#eaf4fd", high = "#1a5276") +
    labs(title = paste0(j, ": edges in common across gene-count cutoffs")) +
    theme(legend.position = "none")
  ggsave(sprintf("results/networks/edge_venn_%s.png", j), p7, width = 7, height = 6, dpi = 300)
}

# --- 8. Community network per 700-gene joint network --------------------------------------------
for (j in joint_ids) {
  g <- as_graph(read_adj(j, 700))
  cm <- read.csv(sprintf("results/communities_%s_top700.csv", j))
  V(g)$community <- factor(cm$community[match(V(g)$name, cm$gene)])
  V(g)$degree <- degree(g)
  n_comm <- sum(table(cm$community) >= 10)
  g <- induced_subgraph(g, which(V(g)$degree > 0))
  set.seed(SEED)
  p8 <- ggraph(g, layout = "fr") +
    geom_edge_link(color = "grey80", alpha = 0.5, width = 0.3) +
    geom_node_point(aes(color = community, size = degree), alpha = 0.9) +
    scale_color_viridis_d(option = "turbo") + scale_size_continuous(range = c(1, 6)) +
    labs(title = paste0(j, " community structure - top 700 genes"),
         subtitle = paste0(n_comm, " communities with >= 10 genes")) +
    theme_void() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave(sprintf("results/networks/community_network_%s_top700.png", j), p8, width = 9, height = 9, dpi = 300)
}
