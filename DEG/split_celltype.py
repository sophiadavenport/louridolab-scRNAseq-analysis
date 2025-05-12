import scanpy as sc
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True, help="Path to counts file for split by tissue")
parser.add_argument("--cell_type_column", required=True, help="Name of cell type column to split data by again")
args = parser.parse_args()

adata = sc.read_h5ad(args.adata)

#creating safe names 
adata.obs[args.cell_type_column] = adata.obs[args.cell_type_column].str.replace(" ", "_", regex=False) 
adata.obs[args.cell_type_column] = adata.obs[args.cell_type_column].str.replace("+", "", regex=False) 
adata.obs[args.cell_type_column] = adata.obs[args.cell_type_column].str.replace("/", "_", regex=False) 

tissues = ['PE', 'Blood']
tissue_dict = {tissue: adata[adata.obs["Tissue"] == tissue] for tissue in tissues}

for tissue, tissue_adata in tissue_dict.items():
    for cell_type in tissue_adata.obs[args.cell_type_column].unique():
        subset = tissue_adata[tissue_adata.obs[args.cell_type_column] == cell_type]
        subset.X=subset.layers['counts']
        
        counts_df = pd.DataFrame(
            subset.X.toarray() if hasattr(subset.X, "toarray") else subset.X,
            index=subset.obs_names,
            columns=subset.var_names
        )

        # counts_parts=args.counts.split("/")
        # counts_path = f"{counts_parts[0]}/{i}_{counts_parts[1]}" #fixing path for the files
        # meta_parts=args.metadata.split("/")
        # meta_path = f"{meta_parts[0]}/{i}_{meta_parts[1]}" #fixing path for the files
        safe_cell_type = cell_type.replace(" ", "_").replace("+", "").replace("/", "_")

        counts_path = f"intermediate/{tissue}_{safe_cell_type}_counts.csv"
        metadata_path = f"intermediate/{tissue}_{safe_cell_type}_metadata.csv"

        counts_df.T.to_csv(counts_path)

        subset.obs.to_csv(metadata_path)
