# pca_plot.py
import anndata
import matplotlib.pyplot as plt
import scanpy as sc
import sys

if len(sys.argv) != 3:
    print("Usage: python qc_filtering.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

data = anndata.read_h5ad(input_file)
sc.tl.pca(data, svd_solver='arpack')
sc.pl.pca(data, save=output_file)
