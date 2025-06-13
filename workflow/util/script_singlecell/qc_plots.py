# qc_plots.py
import anndata
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 3:
    print("Usage: python qc_filtering.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

data = anndata.read_h5ad(input_file )
fig, ax = plt.subplots()
ax.scatter(data.X.sum(axis=1), (data.X > 0).sum(axis=1), alpha=0.5)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
plt.xscale('log')
plt.savefig(output_file)