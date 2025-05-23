# SepSolve: Optimal Marker Gene Selection via *c*-Separation

**Paper:** *Optimal marker genes for c‑separated cell types*

**Authors:** Bartol Borozan, Luka Borozan, Domagoj Ševerdija, Domagoj Matijević, Stefan Canzar
 

---

**SepSolve** is a combinatorial feature selection method for identifying optimal marker genes in single-cell RNA-seq data. Unlike traditional differential expression methods or greedy combinatorial approaches, SepSolve uses a linear programming (LP) framework to select a small set of genes that achieve *c-separation* between all cell types — ensuring robust discrimination that accounts for intra-cell-type variability.

Link to preprint: [biorxiv](https://www.biorxiv.org/content/10.1101/2025.02.12.637849v1)

Link to RECOMB 2025 submission: [Research in Computational Molecular Biology](https://link.springer.com/chapter/10.1007/978-3-031-90252-9_53)

---

## 🔍 What is c-separation?

Two cell types are said to be *c-separated* in a gene subspace if the distance between their mean expression vectors exceeds a multiple (*c*) of their internal variance. SepSolve formalizes this in an optimisation problem that:

* Accounts for expression variability
* Balances separation between all type pairs
* Yields highly stable and compact marker sets

---

## 🔧 Installation

SepSolve is pure‑Python (≥3.9) and only depends on `numpy`, `scipy`, and the commercial (free academic) solver **Gurobi**.

```bash
# pip install sepsolve # from PyPI – coming soon ✨
# or clone the repo
pip install -e .
```

> **Note:** Gurobi ≥10.0 must be available on your system. Follow the [Gurobi quick‑start](https://www.gurobi.com/documentation/) to obtain an academic licence.

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
sepsolve.get_markers(data, labels, num_markers, c=0.4, ilp=False)
```

* `data`: Preprocessed gene expression matrix (cells × genes)
* `labels`: Cluster or cell type annotations
* `num_markers`: Number of marker genes to select
* `c`: (Optional) Separation parameter, default `0.4`
* `ilp`: (Optional) Use integer linear programming (ILP) instead of a faster LP relaxation

---

## 📝 Citation

If you find SepSolve helpful in your research, please cite

```bibtex
@article{Borozan2025SepSolve,
  title   = {Optimal marker genes for c-separated cell types},
  author  = {Borozan, Bartol and Borozan, Luka and \v{S}everdija, Domagoj and Matijevi\'c, Domagoj and Canzar, Stefan},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.02.12.637849}
}
```
or
```bibtex
@InProceedings{10.1007/978-3-031-90252-9_53,
  author="Borozan, Bartol and Borozan, Luka and {\v{S}}everdija, Domagoj and Matijevi{\'{c}}, Domagoj and Canzar, Stefan", 
  editor="Sankararaman, Sriram",
  title="Optimal Marker Genes for c-Separated Cell Types",
  booktitle="Research in Computational Molecular Biology",
  year="2025",
  publisher="Springer Nature Switzerland",
  address="Cham",
  pages="424--427",
  abstract="The identification of cell types in single-cell RNA-seq studies relies on the distinct expression signature of marker genes. A small set of target genes is also needed to design probes for targeted spatial transcriptomics experiments and to target proteins in single-cell spatial proteomics or for cell sorting. While traditional approaches have relied on testing one gene at a time for differential expression between a given cell type and the rest, more recent methods have highlighted the benefits of a joint selection of markers that together distinguish all pairs of cell types simultaneously. However, existing methods either impose constraints on all pairs of individual cells which becomes intractable even for medium-sized datasets, or ignore intra-cell type expression variation entirely by collapsing all cells of a given type to a single representative. Here we address these limitations and propose to find a small set of genes such that cell types are c-separated in the selected dimensions, a notion introduced previously in learning a mixture of Gaussians. To this end, we formulate a linear program that naturally takes into account expression variation within cell types without including each pair of individual cells in the model, leading to a highly stable set of marker genes that allow to accurately discriminate between cell types and that can be computed to optimality efficiently.",
  isbn="978-3-031-90252-9"
}
```

---

## 📜 License

This project is released under the **MIT License** – see the [LICENSE](LICENSE) file for details.

---
