# Interpreting genomic results

!!! info "Section in progress"
    A compact reference for the metrics used across the dashboard: clade /
    lineage labels, genetic distance, divergence (substs/site), % identity,
    genome coverage, PLS-DA VIP scores, mutation burden, and what counts as a
    "novelty signal".

## Quick glossary (will be expanded)

| Term | Meaning in the dashboard |
|---|---|
| **Clade / lineage** | Named phylogenetic group assigned by Nextclade or the build's own labels (`clade_membership` for seasonal flu, `label_clade` for avian flu, `clade` for Ebola) |
| **Genetic distance** | Distance to the closest reference genome — fraction or raw count depending on data availability; elevated values are a novelty signal |
| **Divergence** | Root-to-tip distance in substitutions/site on the tree |
| **MSA identity** | Mean pairwise identity of the query within the multiple-sequence alignment |
| **VIP score** | Variable-importance-in-projection from PLS-DA — how much each protein separates query from background |
