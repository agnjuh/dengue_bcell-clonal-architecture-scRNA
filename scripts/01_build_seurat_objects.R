suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(yaml)
})

samples <- yaml::read_yaml("config/samples.yml")
params  <- yaml::read_yaml("config/params.yml")
qc <- params$qc

dir.create("results/seurat", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

qc_summary <- list()

for (donor in names(samples)) {
  message("Processing donor: ", donor)

  gex_h5 <- samples[[donor]]$gex_h5
  group  <- samples[[donor]]$group

  if (!file.exists(gex_h5)) stop("Missing file: ", gex_h5)

  mat <- Read10X_h5(gex_h5)
  # Read10X_h5 can return a list for multi-assay; for GEX it is often "Gene Expression"
  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) {
      mat <- mat[["Gene Expression"]]
    } else {
      # fallback: take the first assay
      mat <- mat[[1]]
    }
  }

  seu <- CreateSeuratObject(counts = mat, project = paste0("GSE288613_", donor))
  seu$donor <- donor
  seu$group <- group

  # Mito percent (human genes typically start with MT-)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  n_before <- ncol(seu)

  seu <- subset(
    seu,
    subset = nFeature_RNA >= qc$min_features &
             nFeature_RNA <= qc$max_features &
             percent.mt <= qc$max_mito_pct
  )

  n_after <- ncol(seu)

  # Standard workflow (kept simple & robust)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- FindNeighbors(seu, dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:30)

  out_qs <- file.path("results/seurat", paste0(donor, "_gex.qs"))
  qs::qsave(seu, out_qs)

  qc_summary[[donor]] <- tibble(
    donor = donor,
    group = group,
    cells_before_qc = n_before,
    cells_after_qc  = n_after,
    median_nFeature = median(seu$nFeature_RNA),
    median_nCount   = median(seu$nCount_RNA),
    median_pct_mt   = median(seu$percent.mt)
  )
}

qc_df <- bind_rows(qc_summary)
readr::write_csv(qc_df, "results/tables/qc_summary.csv")
message("Wrote: results/tables/qc_summary.csv")
