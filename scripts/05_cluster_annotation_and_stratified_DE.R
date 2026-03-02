suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

merged <- qs::qread("results/seurat/merged_H1H4_S1S4.qs")
merged <- JoinLayers(merged)
DefaultAssay(merged) <- "RNA"

# Ensure cluster ids
merged$cluster <- Idents(merged)

# A) Markers per cluster (top 20)
Idents(merged) <- "cluster"
all_markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

readr::write_csv(all_markers, "results/tables/cluster_top_markers_top20.csv")

# B) DotPlot marker panel for annotation
panel <- c(
  "MS4A1","CD79A","CD74",
  "TCL1A","IGHD","IGHM",
  "CD27","TNFRSF13B","BANK1",
  "CD69","JUN","FOS","NR4A2",
  "PRDM1","XBP1","MZB1","JCHAIN","TNFRSF17","SDC1",
  "MKI67","TOP2A"
)

panel <- panel[panel %in% rownames(merged)]

p <- DotPlot(merged, features = panel, group.by = "cluster") +
  RotatedAxis() +
  ggtitle("Marker panel per cluster (annotation aid)")

ggsave("results/figures/Cluster_marker_panel_dotplot.png", p, width = 11, height = 5)

# C) Stratified DE within high-expansion clusters
# Identify clusters with >=1% expanded cells
meta <- merged@meta.data
meta$cluster <- as.character(meta$cluster)
meta$is_expanded <- ifelse(meta$is_expanded == TRUE, "Expanded", "Non-expanded")

cluster_exp <- meta %>%
  group_by(cluster) %>%
  summarise(frac_expanded = mean(is_expanded == "Expanded"), .groups="drop") %>%
  filter(frac_expanded >= 0.01) %>%
  arrange(desc(frac_expanded))

readr::write_csv(cluster_exp, "results/tables/high_expansion_clusters.csv")

high_clusters <- cluster_exp$cluster

# For each cluster: DE Expanded vs Non-expanded within cluster
strat_DE <- lapply(high_clusters, function(cl) {
  sub <- subset(merged, subset = cluster == cl)
  sub$is_expanded2 <- ifelse(sub$is_expanded == TRUE, "Expanded", "Non-expanded")
  Idents(sub) <- "is_expanded2"

  # If no expanded or too few, skip safely
  n_exp <- sum(Idents(sub) == "Expanded")
  n_non <- sum(Idents(sub) == "Non-expanded")
  if (n_exp < 10 || n_non < 50) return(NULL)

  de <- FindMarkers(sub, ident.1="Expanded", ident.2="Non-expanded",
                    logfc.threshold = 0.25, min.pct = 0.1) %>%
    rownames_to_column("gene") %>%
    mutate(cluster = cl) %>%
    arrange(p_val_adj)

  de
})

strat_DE <- bind_rows(strat_DE)
readr::write_csv(strat_DE, "results/tables/DE_expanded_vs_nonexpanded_stratified.csv")
