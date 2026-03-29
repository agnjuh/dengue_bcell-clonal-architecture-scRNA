suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(patchwork)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/seurat", recursive = TRUE, showWarnings = FALSE)

donors <- c("H1","H4","S1","S4")

objs <- lapply(donors, function(d) {
  qs::qread(file.path("results/seurat", paste0(d, "_gex_vdj.qs")))
})
names(objs) <- donors

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = donors)

# Recompute a joint embedding (robust & interpretable)
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
merged <- RunUMAP(merged, dims = 1:30)

qs::qsave(merged, "results/seurat/merged_H1H4_S1S4.qs")

# UMAP by group
p1 <- DimPlot(merged, group.by = "group", reduction = "umap") +
  ggtitle("B cells – group (H vs S)")
ggsave("results/figures/UMAP_group_H_vs_S.png", p1, width=6, height=5)

# UMAP by expansion
if (!"is_expanded" %in% colnames(merged@meta.data)) {
  stop("is_expanded not found in merged metadata. Check VDJ integration step.")
}

merged$is_expanded <- factor(
  ifelse(
    is.na(merged$is_expanded),
    "Non-expanded",
    ifelse(merged$is_expanded, "Expanded", "Non-expanded")
  ),
  levels = c("Non-expanded", "Expanded")
)

print(table(merged$is_expanded, useNA = "ifany"))

p2 <- DimPlot(merged, group.by = "is_expanded", reduction = "umap") +
  ggtitle("B cells – expanded vs non-expanded clones")

ggsave("results/figures/UMAP_expanded.png", p2, width = 6, height = 5)


# V gene usage
meta <- merged@meta.data

v_usage <- meta %>%
  filter(has_igh) %>%
  count(group, v_gene) %>%
  group_by(group) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

readr::write_csv(v_usage, "results/tables/v_gene_usage.csv")

top_v <- v_usage %>%
  group_by(v_gene) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:15) %>%
  pull(v_gene)

p3 <- v_usage %>%
  filter(v_gene %in% top_v) %>%
  ggplot(aes(x=reorder(v_gene, freq), y=freq, fill=group)) +
  geom_col(position="dodge") +
  coord_flip() +
  labs(title="Top IGHV usage (H vs S)", x="V gene", y="Frequency")

ggsave("results/figures/V_gene_usage_top15.png", p3, width=7, height=6)

# Shannon diversity per donor
shannon <- meta %>%
  filter(has_igh) %>%
  count(donor, clone_id, name = "n") %>%
  group_by(donor) %>%
  summarise(
    shannon = -sum((n/sum(n)) * log(n/sum(n))),
    .groups = "drop"
  )

readr::write_csv(shannon, "results/tables/shannon_diversity.csv")
