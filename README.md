This repository is designed for analyzing single-cell RNA-seq (scRNA-seq) data that has been preprocessed and saved in .h5ad format. It expects two different tissue types (Blood and PE) stored in 'tissue', condition information stored in 'Strains', and celltype annotation in 'Major'. It performs cell type proportion analysis across different tissue types and conducts pairwise condition comparisons.

Differential gene expression analysis is carried out using a pseudobulk approach with edgeR, applied to each major cell type. For each comparison, the pipeline:

Identifies significantly differentially expressed genes, defined as those with an adjusted p-value ≤ 0.05 and absolute log2 fold change (|logFC|) ≥ 1.

Generates heatmaps of significant DEGs for each cell type.

Creates dot plots summarizing significant gene expression differences for each tissue type.

Outputs an Excel file containing all DEG results. 
