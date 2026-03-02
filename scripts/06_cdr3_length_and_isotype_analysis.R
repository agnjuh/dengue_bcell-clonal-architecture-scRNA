suppressPackageStartupMessages({
  library(tidyverse)
  library(qs)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load per-cell VDJ mapping table if saved, otherwise load merged seurat meta
merged <- qs::qread("results/seurat/merged_H1H4_S1S4.qs")
meta <- merged@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(
    cluster = as.character(seurat_clusters),
    is_expanded = ifelse(is_expanded == TRUE, "Expanded", "Non-expanded")
  )

# these columns from VDJ integration step:
required <- c("chain","cdr3","cdr3_nt","v_gene","j_gene","c_gene","has_igh")
missing <- setdiff(required, colnames(meta))
if (length(missing) > 0) {
  stop("Missing columns in merged meta: ", paste(missing, collapse = ", "),
       "\nThis script expects VDJ fields to be stored in Seurat meta.data.")
}

# CDR3 AA length
meta <- meta %>%
  mutate(
    cdr3_aa_len = ifelse(!is.na(cdr3) & cdr3 != "", nchar(cdr3), NA_integer_)
  )

# Focus clusters with expansion signal
focus_clusters <- c("11","10","8")


# 1) Heavy-chain CDR3 length
# 
heavy <- meta %>%
  filter(chain == "IGH", !is.na(cdr3_aa_len))

readr::write_csv(heavy, "results/tables/cdr3_lengths_IGH_per_cell.csv")

p1 <- heavy %>%
  filter(cluster %in% focus_clusters) %>%
  ggplot(aes(x = cdr3_aa_len, fill = is_expanded)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~cluster, nrow = 1) +
  labs(title = "IGH CDR3 amino-acid length (clusters 11/10/8)",
       x = "CDR3 AA length", y = "Density")

ggsave("results/figures/CDR3_len_IGH_density_by_cluster.png", p1, width = 11, height = 4)

p2 <- heavy %>%
  ggplot(aes(x = is_expanded, y = cdr3_aa_len)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  labs(title = "IGH CDR3 AA length: Expanded vs Non-expanded",
       x = NULL, y = "CDR3 AA length")

ggsave("results/figures/CDR3_len_IGH_violin.png", p2, width = 6, height = 5)


# 2) Light-chain CDR3 length
# 
light <- meta %>%
  filter(chain %in% c("IGK","IGL"), !is.na(cdr3_aa_len)) %>%
  mutate(chain = factor(chain, levels = c("IGK","IGL")))

readr::write_csv(light, "results/tables/cdr3_lengths_IGK_IGL_per_cell.csv")

p3 <- light %>%
  filter(cluster %in% focus_clusters) %>%
  ggplot(aes(x = cdr3_aa_len, fill = is_expanded)) +
  geom_density(alpha = 0.4) +
  facet_grid(chain ~ cluster) +
  labs(title = "Light-chain CDR3 AA length (clusters 11/10/8)",
       x = "CDR3 AA length", y = "Density")

ggsave("results/figures/CDR3_len_light_density_by_cluster.png", p3, width = 11, height = 6)


# 3) Isotype usage (heavy only)
#
# c_gene examples: IGHM, IGHD, IGHA1, IGHG1...
isotype <- meta %>%
  filter(chain == "IGH", !is.na(c_gene), c_gene != "") %>%
  mutate(isotype = str_replace(c_gene, "[0-9]+$", "")) %>%  # IGHA1 -> IGHA
  count(is_expanded, isotype) %>%
  group_by(is_expanded) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

readr::write_csv(isotype, "results/tables/isotype_usage_expanded_vs_nonexpanded.csv")

p4 <- ggplot(isotype, aes(x = isotype, y = freq, fill = is_expanded)) +
  geom_col(position = "dodge") +
  labs(title = "Isotype usage (IGH) – Expanded vs Non-expanded",
       x = NULL, y = "Frequency")

ggsave("results/figures/Isotype_usage_expanded.png", p4, width = 7, height = 4)

