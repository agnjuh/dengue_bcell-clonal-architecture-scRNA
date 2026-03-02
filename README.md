# B cell clonal architecture and somatic diversification across H1, H4, S1 and S4

This repository presents an integrated single-cell analysis of B cell receptor repertoires across four samples (H1, H4, S1, S4). The framework combines clonotype reconstruction, somatic hypermutation (SHM) quantification, V gene usage profiling, isotype distribution and CDR3 amino-acid length analysis within a transcriptionally resolved landscape.

Rather than treating repertoire features independently, the analysis situates clonal expansion and diversification within defined transcriptional states, enabling joint interpretation of differentiation status, mutation burden and clonal architecture.

---

## Global transcriptional structure

The global UMAP embedding (`UMAP_group_H_vs_S.png`) shows substantial overlap between H and S groups, suggesting shared transcriptional structure with quantitative modulation rather than discrete condition-specific segregation.

![UMAP_group_H_vs_S.png](figures/UMAP_group_H_vs_S.png)

When restricting the embedding to expanded clones (`UMAP_expanded.png`), clonal expansion appears distributed across multiple transcriptional states. This indicates that proliferative dominance is not confined to a single phenotypic compartment.

![UMAP_expanded.png](figures/UMAP_expanded.png)

Cluster-level marker structure (`Cluster_marker_panel_dotplot.png`) confirms coherent transcriptional identities across clusters and provides context for interpreting clonal enrichment patterns.

![Cluster_marker_panel_dotplot.png](figures/Cluster_marker_panel_dotplot.png)


## Somatic hypermutation (SHM)

Per-sample SHM distributions (`SHM_box_by_sample.png`) reveal systematic differences between expanded and non-expanded clones. In several samples, expanded clones exhibit elevated median SHM levels, consistent with affinity maturation.

![SHM_box_by_sample.png](figures/SHM_box_by_sample.png)

Distributional differences are further resolved using ECDF curves (`SHM_ecdf_by_sample.png`), which allow comparison of full distributional shifts rather than relying solely on central tendency metrics.

![SHM_ecdf_by_sample.png](figures/SHM_ecdf_by_sample.png)

The global relationship between clone size and diversification is shown in `ALL.clone_size_vs_shm.scatter.png`. This representation evaluates whether increasing clonal dominance is associated with accumulated mutation burden.

![ALL.clone_size_vs_shm.scatter.png](figures/ALL.clone_size_vs_shm.scatter.png)

Sample-specific relationships are shown below:

![H1.clone_size_vs_shm.scatter.png](figures/H1.clone_size_vs_shm.scatter.png)
![H4.clone_size_vs_shm.scatter.png](figures/H4.clone_size_vs_shm.scatter.png)
![S1.clone_size_vs_shm.scatter.png](figures/S1.clone_size_vs_shm.scatter.png)
![S4.clone_size_vs_shm.scatter.png](figures/S4.clone_size_vs_shm.scatter.png)


## Isotype composition and V gene usage

Isotype distribution stratified by expansion state (`Isotype_usage_expanded.png`) reveals compositional shifts associated with clonal dominance, consistent with differential class-switch dynamics.

![Isotype_usage_expanded.png](figures/Isotype_usage_expanded.png)

The top IGHV gene usage profile (`V_gene_usage_top15.png`) highlights non-uniform V gene representation between H and S groups. Differences in IGHV3-23 and related segments suggest biased repertoire selection rather than random sampling effects.

![V_gene_usage_top15.png](figures/V_gene_usage_top15.png)


## CDR3 architecture

IGH CDR3 amino-acid length distributions (`CDR3_len_IGH_violin.png`) demonstrate variability in junctional architecture across clusters.

![CDR3_len_IGH_violin.png](figures/CDR3_len_IGH_violin.png)

Cluster-resolved density profiles (`CDR3_len_IGH_density_by_cluster.png`) refine this observation and identify cluster-specific structural tendencies.

![CDR3_len_IGH_density_by_cluster.png](figures/CDR3_len_IGH_density_by_cluster.png)

Light chain CDR3 distributions (`CDR3_len_light_density_by_cluster.png`) show coordinated yet distinct structural patterns relative to heavy chain architecture.

![CDR3_len_light_density_by_cluster.png](figures/CDR3_len_light_density_by_cluster.png)


## Cluster-level clonal expansion

The fraction of expanded cells per transcriptional cluster (`Cluster_expansion_fraction.png`) quantifies enrichment of proliferative clones within defined cellular states.

![Cluster_expansion_fraction.png](figures/Cluster_expansion_fraction.png)

Group-stratified cluster expansion (`Cluster_group_expansion.png`) reveals differential clonal dominance patterns between H and S samples, suggesting condition-specific modulation of expansion within shared transcriptional states.

![Cluster_group_expansion.png](figures/Cluster_group_expansion.png)


## Summary

Across all four samples, clonal expansion, SHM burden, isotype composition and V gene usage exhibit structured, non-random organisation within the transcriptional landscape. Expanded clones are not restricted to a single phenotype but occupy multiple functional states, with mutation accumulation and V gene bias contributing to repertoire shaping. This integrated approach enables simultaneous evaluation of transcriptional identity and adaptive diversification dynamics.