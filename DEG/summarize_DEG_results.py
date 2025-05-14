import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", nargs='+', required=True, help="List of DEG result CSVs.")
parser.add_argument("--output_dot", required=True, help="Path to save the dotplot")
parser.add_argument("--output_report", required=True, help="Path to save the .txt report")

args = parser.parse_args()
summary_lines = []
list_celltype_results=[]
for file in args.input:
    df = pd.read_csv(file)
    tag_raw = os.path.basename(file).split('_edger_results')[0]
    tag = re.sub(r'\b(Blood|PE)\b', '', tag_raw.replace('_', ' ')).strip()
    df_sig = df[(abs(df['logFC']) >= 1) & (df['FDR'] < 0.05)] #filter for significant genes
    counts = df_sig.groupby('comparison')['gene'].count()
    counts.index = counts.index.str.replace('_', ' ').str.replace('.', ' ') #clean up comparison names
    new_df = pd.DataFrame([counts], index=[tag])
    list_celltype_results.append(new_df)

    #generating report info:
    summary_lines.append(f"{tag}")
    grouped = df_sig.groupby('comparison')
    for comparison, group in grouped:
        pos_count = (group['logFC'] > 0).sum()
        neg_count = (group['logFC'] < 0).sum()
        comp_clean = comparison.replace('_', ' ').replace('.', ' ')
        summary_lines.append(f"{comp_clean}")
        summary_lines.append(f"    +: {pos_count}")
        summary_lines.append(f"    -: {neg_count}")
    summary_lines.append("")
with open(args.output_report, 'w') as f:
    f.write('\n'.join(summary_lines))

merged_df = pd.concat(list_celltype_results, sort=False).fillna(0).astype(int)
merged_df = merged_df.sort_index()

tissue = re.search(r'/([^/_]+)_DEG_dotplot', args.output_dot).group(1)

fig, ax = plt.subplots(figsize=(10.5, 8))
max_val = 356
#Flatten the DataFrame for plotting
for i, row_label in enumerate(merged_df.index):
    for j, col_label in enumerate(merged_df.columns):
        value = merged_df.loc[row_label, col_label]
        if value > 0:  #Only plot non-zero values
            dot_size = (value / max_val) * 5000
            ax.scatter(j, i, s=dot_size, color='blue', alpha=0.6)
            if value > 10:
                ax.text(j, i, str(value), ha='center', va='center', fontsize=9, color='white')

# Set ticks and labels
ax.set_xticks(range(len(merged_df.columns)))
ax.set_xticklabels(merged_df.columns, ha='center')
ax.set_yticks(range(len(merged_df.index)))
ax.set_yticklabels(merged_df.index)

ax.set_xlim(-0.5, len(merged_df.columns) - 0.5)
ax.set_ylim(-0.5, len(merged_df.index) - 0.5)
ax.invert_yaxis()
ax.set_title(f"{tissue}: Significant Genes")
ax.set_xlabel("Comparison")
ax.set_ylabel("Major Cell Type")
plt.tight_layout()
plt.savefig(args.output_dot)