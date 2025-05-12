library(tidyverse)
library(limma)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
metadata_file <- args[2]
condition_col <- args[3]
covariates_cols <- args[4]
output_file <- args[5]

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
print(max(counts))
metadata <- read.csv(metadata_file, row.names = 1)
metadata[[condition_col]] <- gsub(" ", ".", metadata[[condition_col]]) #trying to fix the " " issue

group <- as.factor(metadata[[condition_col]])

print(levels(group))
print(table(group))

dge <- DGEList(counts = counts, group = group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

group_levels <- levels(group)
pairwise_comparisons <- combn(group_levels, 2, simplify = FALSE)

all_results <- list()

for (pair in pairwise_comparisons) {
  contrast_string <- paste0(pair[1], " - ", pair[2])
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
  fit2 <- glmQLFTest(fit, contrast = contrast_matrix)
  table <- topTags(fit2, n = Inf)$table
  table$comparison <- paste0(pair[1], "_vs_", pair[2])
  table$gene <- rownames(table)
  all_results[[paste(pair, collapse = "_vs_")]] <- table
}

final_results <- bind_rows(all_results)
write.csv(final_results, output_file, row.names = TRUE)

# run_edgeR <- function(count_dataframe, group_vector, contrast_levels = NULL, use_qlf = TRUE) {
#   library(limma)
#   library(edgeR)

#   group_factor <- factor(group_vector)
#   if (length(group_factor) != ncol(count_dataframe)) {
#     stop("Length of group vector doesn't match the number of barcodes/cells (columns) in counts dataframe.")
#   }

#   dge <- DGEList(counts = count_dataframe, group = group_factor)
#   keep <- filterByExpr(dge, group = group_factor) #Filter lowly expressed genes
#   dge <- dge[keep, , keep.lib.sizes = FALSE]
#   dge <- calcNormFactors(dge)

#   if (nlevels(group_factor) == 2 && is.null(contrast_levels)) {
#     dge <- estimateDisp(dge)
#     result <- exactTest(dge)
#     return(topTags(result, n = Inf)$table)
#   }

#   design <- model.matrix(~0 + group_factor)
#   colnames(design) <- levels(group_factor)
#   dge <- estimateDisp(dge, design)
#   fit <- if (use_qlf) {
#     glmQLFit(dge, design)
#   } else {
#     glmFit(dge, design)
#   }

#   if (!is.null(contrast_levels)) {
#     if (length(contrast_levels) != 2) {
#       stop("contrast_levels must be a vector of two group names, e.g., c('treated', 'control')")
#     }
#     contrast <- makeContrasts(contrasts = paste0(contrast_levels[1], " - ", contrast_levels[2]),
#                               levels = design)
    
#     result <- if (use_qlf) {
#       glmQLFTest(fit, contrast = contrast)
#     } else {
#       glmLRT(fit, contrast = contrast)
#     }
    
#     return(topTags(result, n = Inf)$table)
#   } else {
#     result <- if (use_qlf) {
#       glmQLFTest(fit)
#     } else {
#       glmLRT(fit)
#     }
    
#     return(topTags(result, n = Inf)$table)
#   }
# }

# results <- run_edgeR(counts, group)
# write.csv(results, output_file)