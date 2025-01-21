# qc_filtering.py
import anndata
import scanpy as sc
import sys

if len(sys.argv) != 3:
    print("Usage: python qc_filtering.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

data = anndata.read_h5ad(input_file)
sc.pp.filter_cells(data, min_genes=200)
knee_threshold = 200
data = data[data.obs.n_genes > knee_threshold]
data.write(output_file)