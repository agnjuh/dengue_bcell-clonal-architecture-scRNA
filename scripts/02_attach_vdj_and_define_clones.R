suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(qs)
  library(yaml)
})

samples <- yaml::read_yaml("config/samples.yml")
params  <- yaml::read_yaml("config/params.yml")
vdj_p   <- params$vdj

dir.create("results/repertoire", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

infer_group <- function(donor) {
  # donor like "H1", "S4"
  substr(donor, 1, 1)
}

normalize_barcode <- function(x) x

pick_top_contig_per_barcode <- function(df) {
  df %>%
    group_by(barcode) %>%
    arrange(desc(umis), desc(reads)) %>%
    slice(1) %>%
    ungroup()
}

per_donor_stats <- list()

for (donor in names(samples)) {
  message("Integrating VDJ for donor: ", donor)

  gex_qs <- file.path("results/seurat", paste0(donor, "_gex.qs"))
  if (!file.exists(gex_qs)) stop("Missing Seurat object: ", gex_qs)

  # Support both:
  # samples[[donor]]$vdj_csv_gz OR samples[[donor]] as a named vector
  vdj_file <- NULL
  if (is.list(samples[[donor]]) && !is.null(samples[[donor]]$vdj_csv_gz)) {
    vdj_file <- samples[[donor]]$vdj_csv_gz
  } else if (!is.null(samples[[donor]][["vdj_csv_gz"]])) {
    vdj_file <- samples[[donor]][["vdj_csv_gz"]]
  } else {
    stop("Could not find 'vdj_csv_gz' for donor ", donor, " in config/samples.yml")
  }
  if (!file.exists(vdj_file)) stop("Missing VDJ file: ", vdj_file)

  group <- infer_group(donor)

  seu <- qs::qread(gex_qs)

  vdj <- data.table::fread(vdj_file) %>% as_tibble()

  if (vdj_p$only_is_cell && "is_cell" %in% names(vdj)) {
    vdj <- vdj %>% filter(is_cell == TRUE)
  }
  if (vdj_p$only_productive && "productive" %in% names(vdj)) {
    vdj <- vdj %>% filter(productive == TRUE)
  }
  if (!is.null(vdj_p$chain) && "chain" %in% names(vdj)) {
    vdj <- vdj %>% filter(chain == vdj_p$chain)
  }

  vdj <- vdj %>%
    mutate(barcode = normalize_barcode(barcode)) %>%
    pick_top_contig_per_barcode()

  required_cols <- c("barcode", "v_gene", "j_gene", "cdr3")
  missing <- setdiff(required_cols, names(vdj))
  if (length(missing) > 0) stop("VDJ file missing columns: ", paste(missing, collapse = ", "))

  vdj <- vdj %>%
    mutate(
      clone_id = paste(v_gene, j_gene, cdr3, sep = "|"),
      cdr3_len = nchar(cdr3)
    ) %>%
    select(
      barcode, clone_id, v_gene, j_gene, c_gene, d_gene, cdr3, cdr3_len,
      reads, umis, raw_clonotype_id, raw_consensus_id, exact_subclonotype_id
    )

  md <- tibble(barcode = normalize_barcode(colnames(seu))) %>%
    left_join(vdj, by = "barcode") %>%
    mutate(has_igh = !is.na(clone_id))

  expanded_min <- vdj_p$expanded_clone_min_cells
  md <- md %>%
    group_by(clone_id) %>%
    mutate(clone_size = ifelse(is.na(clone_id), 0L, n())) %>%
    ungroup() %>%
    mutate(is_expanded = has_igh & (clone_size >= expanded_min))

  # Add metadata in correct cell order
  md_ordered <- md %>%
    mutate(cell = colnames(seu)) %>%
    select(
      cell, has_igh, clone_id, clone_size, is_expanded,
      v_gene, j_gene, c_gene, d_gene, cdr3, cdr3_len,
      reads, umis, raw_clonotype_id, raw_consensus_id, exact_subclonotype_id
    ) %>%
    as.data.frame()
  rownames(md_ordered) <- md_ordered$cell
  md_ordered$cell <- NULL

  seu <- AddMetaData(seu, metadata = md_ordered)

  out_qs <- file.path("results/seurat", paste0(donor, "_gex_vdj.qs"))
  qs::qsave(seu, out_qs)

  donor_rep <- md %>%
    filter(has_igh) %>%
    count(clone_id, sort = TRUE, name = "cells_in_clone") %>%
    mutate(donor = donor, group = group)

  readr::write_csv(donor_rep, file.path("results/repertoire", paste0(donor, "_clone_sizes.csv")))

  per_donor_stats[[donor]] <- tibble(
    donor = donor,
    group = group,
    cells_total = nrow(md),
    cells_with_igh = sum(md$has_igh, na.rm = TRUE),
    frac_with_igh = mean(md$has_igh, na.rm = TRUE),
    expanded_cells = sum(md$is_expanded, na.rm = TRUE),
    frac_expanded = mean(md$is_expanded, na.rm = TRUE)
  )
}

stats_df <- bind_rows(per_donor_stats)
readr::write_csv(stats_df, "results/tables/vdj_mapping_summary.csv")
message("Wrote: results/tables/vdj_mapping_summary.csv")
