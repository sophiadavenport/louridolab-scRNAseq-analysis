#Intended to take a scanpy object and move it to matrix format for R:
import scanpy as sc
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to .h5ad file")
parser.add_argument("--counts", required=True, help="Output CSV for counts matrix")
parser.add_argument("--metadata", required=True, help="Output CSV for sample metadata")
parser.add_argument("--split_column", required=True, help="Column in metadata to split by")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.X=adata.layers['counts']

for i in adata.obs[args.split_column].unique():
    subset = adata[adata.obs[args.split_column] == i]

    counts_df = pd.DataFrame(
        subset.X.toarray() if hasattr(subset.X, "toarray") else subset.X,
        index=subset.obs_names,
        columns=subset.var_names
    )

    counts_parts=args.counts.split("/")
    counts_path = f"{counts_parts[0]}/{i}_{counts_parts[1]}" #fixing path for the files
    meta_parts=args.metadata.split("/")
    meta_path = f"{meta_parts[0]}/{i}_{meta_parts[1]}" #fixing path for the files

    counts_df.T.to_csv(counts_path)

    subset.obs.to_csv(meta_path)


