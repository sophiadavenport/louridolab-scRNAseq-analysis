library(tidyverse)
library(limma)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
metadata_file <- args[2]
condition_col <- args[3]
output_file <- args[4]

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
print(max(counts))
metadata <- read.csv(metadata_file, row.names = 1)
metadata[[condition_col]] <- gsub(" ", ".", metadata[[condition_col]])

group <- as.factor(metadata[[condition_col]])

print(levels(group))
print(table(group))

#assuming column 'replicates'...
replicate_covariates <- metadata %>%
  group_by(replicate) %>%
  summarise(
    num_cells = n(),
    avg_mt = mean(percent.mt, na.rm = TRUE),
    avg_features = mean(nFeature_RNA, na.rm = TRUE)
  )

print(replicate_covariates)

metadata <- metadata %>%
  left_join(replicate_covariates, by = "replicate")

dge <- DGEList(counts = counts, group = group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge) 

initial_covariates <- c("Mouse", "avg_features", "avg_mt")
covariates <- initial_covariates
removed_covariates <- c()
design_success <- FALSE

while (!design_success && length(covariates) >= 0) {
  formula_str <- paste("~ 0 + group", if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else "")
  message("Trying design matrix with formula: ", formula_str)

  design <- model.matrix(as.formula(formula_str), data = metadata)
  colnames(design)[seq_along(levels(group))] <- levels(group)

  fit_try <- try({
    dge <- estimateDisp(dge, design)
    glmQLFit(dge, design)
  }, silent = TRUE)

  if (inherits(fit_try, "try-error")) {
    message("Model fitting failed with covariates: ", paste(covariates, collapse = ", "))
    if (length(covariates) == 0) {
      stop("All covariates removed and model still not estimable. Aborting.")
    }
    removed_covariate <- tail(covariates, 1)
    removed_covariates <- c(removed_covariates, removed_covariate)
    covariates <- head(covariates, -1)
    message("Removed covariate: ", removed_covariate)
  } else {
    fit <- fit_try
    design_success <- TRUE
    message("Model fitting succeeded with covariates: ", paste(covariates, collapse = ", "))
  }
}

if (length(removed_covariates) > 0) {
  message("Final model excludes these covariates due to rank deficiency: ", paste(removed_covariates, collapse = ", "))
}

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