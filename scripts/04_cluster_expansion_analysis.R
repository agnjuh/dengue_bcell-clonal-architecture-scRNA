suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

merged <- qs::qread("results/seurat/merged_H1H4_S1S4.qs")

meta <- merged@meta.data
if (!"is_expanded" %in% colnames(meta)) stop("is_expanded column missing")

# Add cluster + expansion status
meta$cluster <- as.character(Idents(merged))
meta$is_expanded <- ifelse(meta$is_expanded == TRUE, "Expanded", "Non-expanded")


# 1) Expansion fraction per cluster
# 
cluster_expansion <- meta %>%
  group_by(cluster) %>%
  summarise(
    total = n(),
    expanded = sum(is_expanded == "Expanded"),
    frac_expanded = expanded / total,
    .groups = "drop"
  )

readr::write_csv(cluster_expansion, "results/tables/cluster_expansion_fraction.csv")

p1 <- ggplot(cluster_expansion,
             aes(x = reorder(cluster, frac_expanded),
                 y = frac_expanded)) +
  geom_col() +
  coord_flip() +
  labs(title = "Fraction of expanded cells per cluster",
       x = "Cluster",
       y = "Fraction expanded")

ggsave("results/figures/Cluster_expansion_fraction.png", p1, width = 6, height = 5)


# 2) Cluster × group interaction
# 
cluster_group_expansion <- meta %>%
  group_by(cluster, group) %>%
  summarise(
    total = n(),
    expanded = sum(is_expanded == "Expanded"),
    frac_expanded = expanded / total,
    .groups = "drop"
  )

readr::write_csv(cluster_group_expansion, "results/tables/cluster_group_expansion.csv")

p2 <- ggplot(cluster_group_expansion,
             aes(x = reorder(cluster, frac_expanded),
                 y = frac_expanded,
                 fill = group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Expanded fraction per cluster (H vs S)",
       x = "Cluster",
       y = "Fraction expanded")

ggsave("results/figures/Cluster_group_expansion.png", p2, width = 7, height = 5)


# 3) Differential expression: Expanded vs Non-expanded
#    Seurat v5: need JoinLayers before DE
# 
merged <- JoinLayers(merged)

# Set identities to expansion state
merged$is_expanded <- meta$is_expanded
Idents(merged) <- "is_expanded"

markers_expanded <- FindMarkers(
  merged,
  ident.1 = "Expanded",
  ident.2 = "Non-expanded",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

markers_expanded <- markers_expanded %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)

readr::write_csv(markers_expanded, "results/tables/DE_expanded_vs_nonexpanded.csv")
