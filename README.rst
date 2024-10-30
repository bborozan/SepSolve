Basic Example
----------------------
 .. code-block:: python

    import sepsolve
    import scanpy as sc

    # load dataset
    adata = sc.datasets.paul15()

    # preprocess 
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # get 25 marker genes
    data = adata.X
    labels = adata.obs["paul15_clusters"]
    markers = sepsolve.get_markers(adata.X, labels, 25)