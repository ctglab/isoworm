# umap_clustering.py
import anndata
import scanpy as sc
import sys

if len(sys.argv) != 3:
    print("Usage: python qc_filtering.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

data = anndata.read_h5ad(input_file)
sc.pp.neighbors(data)
sc.tl.umap(data)
sc.tl.louvain(data)
sc.pl.umap(data, color='louvain', save=output_file)
