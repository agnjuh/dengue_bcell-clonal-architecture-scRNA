# B cell clonal architecture and somatic diversification across H1, H4, S1 and S4

This repository presents an integrated analysis of B cell receptor repertoires derived from four samples (H1, H4, S1, S4). The analysis combines clonotype reconstruction, somatic hypermutation profiling, V gene usage, isotype distribution and CDR3 amino-acid length analysis within a single-cell transcriptional framework.

---

## Global transcriptional structure

The global UMAP embedding (`UMAP_group_H_vs_S.png`) shows substantial overlap between H and S groups rather than discrete segregation.

![UMAP_group_H_vs_S.png](results/figures/UMAP_group_H_vs_S.png)

When restricting to expanded clones, `UMAP_expanded.png` indicates that clonal expansion is distributed across multiple transcriptional states.

![UMAP_expanded.png](results/figures/UMAP_expanded.png)


## Somatic hypermutation (SHM)

SHM distributions per sample are shown in `SHM_box_by_sample.png`.

![SHM_box_by_sample.png](results/figures/SHM_box_by_sample.png)

Distributional shifts are further resolved in `SHM_ecdf_by_sample.png`.

![SHM_ecdf_by_sample.png](results/figures/SHM_ecdf_by_sample.png)

The relationship between clone size and diversification is visualised in `ALL.clone_size_vs_shm.scatter.png`.

![ALL.clone_size_vs_shm.scatter.png](results/figures/ALL.clone_size_vs_shm.scatter.png)

Sample-specific scatter plots:

- `H1.clone_size_vs_shm.scatter.png`
- `H4.clone_size_vs_shm.scatter.png`
- `S1.clone_size_vs_shm.scatter.png`
- `S4.clone_size_vs_shm.scatter.png`

![H1.clone_size_vs_shm.scatter.png](results/figures/H1.clone_size_vs_shm.scatter.png)
![H4.clone_size_vs_shm.scatter.png](results/figures/H4.clone_size_vs_shm.scatter.png)
![S1.clone_size_vs_shm.scatter.png](results/figures/S1.clone_size_vs_shm.scatter.png)
![S4.clone_size_vs_shm.scatter.png](results/figures/S4.clone_size_vs_shm.scatter.png)


## Isotype and V gene usage

Isotype distribution across expansion states (`Isotype_usage_expanded.png`):

![Isotype_usage_expanded.png](results/figures/Isotype_usage_expanded.png)

Top IGHV gene usage (`V_gene_usage_top15.png`):

![V_gene_usage_top15.png](results/figures/V_gene_usage_top15.png)


## CDR3 architecture

IGH CDR3 amino-acid length distributions (`CDR3_len_IGH_violin.png`):

![CDR3_len_IGH_violin.png](results/figures/CDR3_len_IGH_violin.png)

Cluster-resolved density (`CDR3_len_IGH_density_by_cluster.png`):

![CDR3_len_IGH_density_by_cluster.png](results/figures/CDR3_len_IGH_density_by_cluster.png)

Light chain CDR3 distributions (`CDR3_len_light_density_by_cluster.png`):

![CDR3_len_light_density_by_cluster.png](results/figures/CDR3_len_light_density_by_cluster.png)

---

## Cluster-level expansion

Fraction of expanded cells per cluster (`Cluster_expansion_fraction.png`):

![Cluster_expansion_fraction.png](results/figures/Cluster_expansion_fraction.png)

Group-stratified cluster expansion (`Cluster_group_expansion.png`):

![Cluster_group_expansion.png](results/figures/Cluster_group_expansion.png)