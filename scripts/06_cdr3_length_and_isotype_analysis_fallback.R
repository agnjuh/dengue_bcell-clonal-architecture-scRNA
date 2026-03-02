suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(qs)
  library(yaml)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load merged seurat meta (state variables)
merged <- qs::qread("results/seurat/merged_H1H4_S1S4.qs")
meta <- merged@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(
    donor = as.character(donor),
    group = as.character(group),
    cluster = as.character(seurat_clusters),
    is_expanded = ifelse(is_expanded == TRUE, "Expanded", "Non-expanded"),
    barcode = sub("^[^_]+_", "", cell)  # remove "H1_" prefix
  ) %>%
  select(cell, donor, group, cluster, is_expanded, barcode)

# Load VDJ contig annotations from config/samples.yml
samples <- yaml::read_yaml("config/samples.yml")

read_vdj <- function(donor, path) {
  message("Reading VDJ: ", donor)
  dt <- data.table::fread(path)
  dt[, donor := donor]
  dt
}

vdj_list <- lapply(names(samples), function(d) {
  read_vdj(d, samples[[d]]$vdj_csv_gz)
})
vdj <- data.table::rbindlist(vdj_list, fill = TRUE) %>%
  as_tibble()

# Keep productive, high_confidence contigs, and real cells
vdj <- vdj %>%
  filter(
    is_cell == TRUE,
    high_confidence == TRUE,
    productive == TRUE
  ) %>%
  mutate(
    chain = as.character(chain),
    cdr3 = as.character(cdr3),
    cdr3_nt = as.character(cdr3_nt),
    barcode = as.character(barcode),
    c_gene = as.character(c_gene),
    v_gene = as.character(v_gene),
    j_gene = as.character(j_gene),
    d_gene = as.character(d_gene)
  )

# Prefer 1 contig per chain per cell:
# choose the contig with highest UMIs (or reads if UMIs missing)
score_col <- if ("umis" %in% colnames(vdj)) "umis" else if ("reads" %in% colnames(vdj)) "reads" else NULL
if (!is.null(score_col)) {
  vdj <- vdj %>%
    group_by(donor, barcode, chain) %>%
    slice_max(order_by = .data[[score_col]], n = 1, with_ties = FALSE) %>%
    ungroup()
} else {
  vdj <- vdj %>%
    group_by(donor, barcode, chain) %>%
    slice(1) %>%
    ungroup()
}

# Join with seurat metadata
df <- vdj %>%
  inner_join(meta, by = c("donor","barcode"))

# Add CDR3 AA length
df <- df %>%
  mutate(
    cdr3_aa_len = ifelse(!is.na(cdr3) & cdr3 != "", nchar(cdr3), NA_integer_),
    isotype = ifelse(chain == "IGH" & !is.na(c_gene) & c_gene != "",
                     stringr::str_replace(c_gene, "[0-9]+$", ""),
                     NA_character_)
  )

readr::write_csv(df, "results/tables/per_cell_vdj_joined_with_states.csv")

focus_clusters <- c("11","10","8")


# 1) Heavy-chain CDR3 length
# 
heavy <- df %>% filter(chain == "IGH", !is.na(cdr3_aa_len))
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
light <- df %>% filter(chain %in% c("IGK","IGL"), !is.na(cdr3_aa_len))
readr::write_csv(light, "results/tables/cdr3_lengths_light_per_cell.csv")

p3 <- light %>%
  filter(cluster %in% focus_clusters) %>%
  ggplot(aes(x = cdr3_aa_len, fill = is_expanded)) +
  geom_density(alpha = 0.4) +
  facet_grid(chain ~ cluster) +
  labs(title = "Light-chain CDR3 AA length (clusters 11/10/8)",
       x = "CDR3 AA length", y = "Density")
ggsave("results/figures/CDR3_len_light_density_by_cluster.png", p3, width = 11, height = 6)


# 3) Isotype usage (IGH only)
#
isotype_tab <- heavy %>%
  filter(!is.na(isotype)) %>%
  count(is_expanded, isotype) %>%
  group_by(is_expanded) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

readr::write_csv(isotype_tab, "results/tables/isotype_usage_expanded_vs_nonexpanded.csv")

p4 <- ggplot(isotype_tab, aes(x = isotype, y = freq, fill = is_expanded)) +
  geom_col(position = "dodge") +
  labs(title = "Isotype usage (IGH) – Expanded vs Non-expanded",
       x = NULL, y = "Frequency")
ggsave("results/figures/Isotype_usage_expanded.png", p4, width = 7, height = 4)

message("Wrote figures to results/figures and tables to results/tables")
