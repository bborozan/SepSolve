# SepSolve

**SepSolve** is a combinatorial feature selection method for identifying optimal marker genes in single-cell RNA-seq data. Unlike traditional differential expression methods or greedy combinatorial approaches, SepSolve uses a linear programming (LP) framework to select a small set of genes that achieve *c-separation* between all cell types — ensuring robust discrimination that accounts for intra-cell-type variability.

Link to preprint: [biorxiv](https://www.biorxiv.org/content/10.1101/2025.02.12.637849v1)

---

## 🔍 What is c-separation?

Two cell types are said to be *c-separated* in a gene subspace if the distance between their mean expression vectors exceeds a multiple (*c*) of their internal variance. SepSolve formalizes this in an optimisation problem that:

* Accounts for expression variability
* Balances separation between all type pairs
* Yields highly stable and compact marker sets

---

## 🧬 Quick Start

Here's a minimal example using the `paul15` dataset from `scanpy`:

```python
import scanpy as sc

# Load dataset
adata = sc.datasets.paul15()

# Preprocess
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Get 25 marker genes that c-separate cell types
data = adata.X
labels = adata.obs["paul15_clusters"]
markers = sepsolve.get_markers(data, labels, 25)

print("Selected marker indices:", markers)
```

> Use `adata.var_names[markers]` to get actual gene names.

---

## ⚙️ API

```python
sepsolve.get_markers(data, labels, n_markers, c=0.4)
```

* `data`: Preprocessed gene expression matrix (cells × genes)
* `labels`: Cluster or cell type annotations
* `n_markers`: Number of marker genes to select
* `c`: (Optional) Separation parameter, default `0.4`

