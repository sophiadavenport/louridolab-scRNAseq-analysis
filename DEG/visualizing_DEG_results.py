import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def build_logFC_matrix(path, group_col='comparison'):
    tag_raw = os.path.basename(path).split('_edger_results')[0]
    tag = tag_raw.replace('_', ' ')
    df = pd.read_csv(path)

    grouped = df.groupby(group_col)

    all_significant_genes = set()
    logFC_dict = {}

    for comparison, group_df in grouped:
        comp_clean = comparison.replace('_', ' ').replace('.', ' ')
        significant = group_df[(abs(group_df['logFC']) >= 1) & (group_df['FDR'] < 0.05)]
        genes = significant['gene'].tolist()
        all_significant_genes.update(genes)
        logFC_dict[comp_clean] = significant.set_index('gene')['logFC'].to_dict()

    all_genes = sorted(all_significant_genes)
    comparisons = sorted(logFC_dict.keys())
    matrix = pd.DataFrame(0.0, index=all_genes, columns=comparisons)

    for comp in comparisons:
        for gene in all_genes:
            if gene in logFC_dict[comp]:
                matrix.at[gene, comp] = logFC_dict[comp][gene]

    return matrix, tag

def plot_logFC_heatmap(matrix, tag, output_path, cluster_rows=True, cluster_cols=True):
    num_genes = matrix.shape[0]
    if num_genes <= 5:
        fig_height = 8
    elif num_genes <= 10:
        fig_height = 6
    else:
        fig_height = min(0.25 * num_genes, 30)
    font_size = max(4, min(8, 12 - 0.03 * num_genes))

    g = sns.clustermap(matrix,
                       cmap="bwr",
                       center=0,
                       row_cluster=cluster_rows,
                       col_cluster=cluster_cols,
                       cbar_kws={'shrink': 0.75},
                       linecolor='lightgray'
                       )
    
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)
    g.ax_heatmap.set_xlabel('Comparison', fontsize=10)
    g.ax_heatmap.set_ylabel('Gene', fontsize=10)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=65, fontsize=8)
    g.ax_heatmap.set_yticks(np.arange(matrix.shape[0]))
    g.ax_heatmap.set_yticklabels(matrix.index, fontsize=font_size)
    g.figure.set_size_inches(10, fig_height)
    g.figure.suptitle(f"{tag} Significant Genes", fontsize=14, y=1.05, ha='center')
    cbar_ax = g.cax
    cbar_ax.set_title('logFC', fontsize=10, pad=10, loc='center')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to DEG .csv file")
    parser.add_argument("--output", required=True, help="Path to save the heatmap image")
    args = parser.parse_args()

    matrix, tag = build_logFC_matrix(args.input)
    if matrix.empty:
        print(f"No significant genes found in: {args.input}")
    else:
        plot_logFC_heatmap(matrix, tag, args.output)
        print(f"Heatmap saved to: {args.output}")

