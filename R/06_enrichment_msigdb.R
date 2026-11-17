# MSigDB over-representation analysis 

suppressPackageStartupMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(msigdbr)
})

source("R/00_config.R")
MIN_COMM_SIZE <- 10

collections <- list(
  Hallmark   = msigdbr(species = "Homo sapiens", collection = "H"),
  C2_CP      = msigdbr(species = "Homo sapiens", collection = "C2") %>%
                 filter(grepl("^CP", gs_subcollection)),
  C8_celltype = msigdbr(species = "Homo sapiens", collection = "C8")
)
collections <- lapply(collections, function(m) {
  m %>% transmute(gs_name, ncbi_gene = as.character(ncbi_gene)) %>% distinct()
})

all_res <- list()

for (jid in joint_ids) {

  comm_tab <- read.csv(paste0("results/communities_", jid, "_top700.csv"))
  map <- bitr(comm_tab$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  universe_entrez <- unique(map$ENTREZID)

  comm_list <- split(comm_tab$gene, comm_tab$community)
  comm_list <- comm_list[lengths(comm_list) >= MIN_COMM_SIZE]

  for (cid in names(comm_list)) {
    entrez <- unique(map$ENTREZID[map$SYMBOL %in% comm_list[[cid]]])

    for (coll in names(collections)) {
      res <- enricher(gene = entrez, universe = universe_entrez,
                      TERM2GENE = collections[[coll]],
                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                      minGSSize = 5, maxGSSize = 500)
      if (is.null(res) || nrow(as.data.frame(res)) == 0) next
      res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      df <- as.data.frame(res)
      all_res[[paste(jid, cid, coll)]] <- data.frame(
        joint = jid, community = cid, n_genes_community = length(comm_list[[cid]]),
        collection = coll, gene_set = df$ID,
        GeneRatio = df$GeneRatio, BgRatio = df$BgRatio,
        fold_enrichment = sapply(df$GeneRatio, function(x) eval(parse(text = x))) /
                          sapply(df$BgRatio, function(x) eval(parse(text = x))),
        pvalue = df$pvalue, p.adjust = df$p.adjust, Count = df$Count,
        genes = df$geneID, row.names = NULL)
    }
  }
}

full <- bind_rows(all_res) %>%
  mutate(significant = p.adjust < PADJ_CUTOFF) %>%
  arrange(joint, community, collection, p.adjust)
write.csv(full, "results/msigdb_enrichment_full_top700.csv", row.names = FALSE)

short <- full %>%
  filter(significant) %>%
  group_by(collection, community, gene_set) %>%
  summarise(
    joints_significant = paste(sub("joint", "", sort(joint)), collapse = ","),
    n_joints = n(),
    genes_in_community = paste(sub("/.*", "", GeneRatio), collapse = "; "),
    genes_in_background = sub("/.*", "", BgRatio[1]),
    fold_enrichment = sprintf("%.2f", mean(fold_enrichment)),
    padj_min = min(p.adjust), padj_max = max(p.adjust),
    shared_genes = paste(Reduce(intersect, strsplit(genes, "/")), collapse = ", "),
    .groups = "drop") %>%
  arrange(factor(collection, c("Hallmark", "C2_CP", "C8_celltype")), community,
          desc(n_joints), padj_min)
write.csv(short, "results/msigdb_enrichment_short_top700.csv", row.names = FALSE)
