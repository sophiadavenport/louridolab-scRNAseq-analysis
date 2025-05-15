import scanpy as sc
from scanpro import scanpro
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("--data_directory", required=True)
args = parser.parse_args()

data_directory = args.data_directory

adata=sc.read_h5ad(data_directory+'lourido_rna_processed.h5ad')

def scanpro_tissue_majorcelltype(adata, tissue_name):
    tissue_adata = adata[adata.obs['Tissue'] == tissue_name]
    out = scanpro(tissue_adata,clusters_col='Major',conds_col='Strains',samples_col='replicate')
    
    #Scanpro Scatterplot results:
    out.plot(save=f"{tissue_name}_proportions_scatterplot.png")

    #Print counts:
    counts = tissue_adata.obs.groupby('Strains').Major.value_counts()
    print(f"{tissue_name} Major Celltype Counts:\n", counts)

    #Plot celltype counts:
    counts_unstack = counts.unstack(fill_value=0)
    counts_unstack.plot(kind='bar',stacked=True, figsize=(10, 6),colormap='tab20')
    plt.ylabel('Cell Count')
    plt.title(f'{tissue_name}: Count of Major Cell Type per Strain')
    plt.legend(title='Major', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{tissue_name}_celltypecounts_barplot.png')
    plt.close()

    #Plot percentages:
    percentages = counts_unstack.div(counts_unstack.sum(axis=1), axis=0) * 100
    percentages.plot(kind='bar', stacked=True, figsize=(10, 6),colormap='tab20')
    plt.ylabel('Percentage of Cells')
    plt.title(f'{tissue_name}: Proportion of Cells per Strain')
    plt.legend(title='Major', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{tissue_name}_proportioncelltypes_barplot.png')
    plt.close()

    return out 

results = {}
for tissue in ['PE', 'Blood']:
    results[tissue] = scanpro_tissue_majorcelltype(adata, tissue)


###Strain Pairwise Comparison:
conditions = ['Type 1', 'Type 2', 'Type 3', 'Uninfected']
tissue_comparison_results={}
PE_adata=adata[adata.obs['Tissue'] == 'PE']
Blood_adata=adata[adata.obs['Tissue'] == 'Blood']
adatas={'PE': PE_adata, 'Blood': Blood_adata}

for tissue in adatas:
    adata=adatas[tissue]
    pairwise_results = {}
    for cond_pair in itertools.combinations(conditions, 2):
        print(f"Running scanpro for {cond_pair[0]} vs {cond_pair[1]}")

        out = scanpro(adata, clusters_col='Major', samples_col='replicate', conds_col='Strains',
            transform='logit', conditions=list(cond_pair), robust=False
        )
    
        pairwise_results[f"{cond_pair[0]}_vs_{cond_pair[1]}"] = out.results
    tissue_comparison_results[tissue]=pairwise_results

with pd.ExcelWriter("tissue_comparison_results.xlsx", engine="openpyxl") as writer:
    for tissue, comparisons in tissue_comparison_results.items():
        row = 0
        for comparison_name, df in comparisons.items():
            sheet = writer.book.create_sheet(title=tissue) if tissue not in writer.book.sheetnames else writer.sheets[tissue]
            worksheet = writer.sheets[tissue]
            worksheet.cell(row=row + 1, column=1).value = comparison_name
            df.to_excel(writer, sheet_name=tissue, startrow=row + 2, index=True)
            row += len(df) + 4