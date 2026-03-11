#!/usr/bin/env Rscript

if (!exists("args", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = TRUE)
}
input_file <- if (length(args) >= 1) args[[1]] else "proteome_data.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "Proteome_DE"

ensure_package <- function(pkg, source = c("CRAN", "BIOC")) {
  source <- match.arg(source)
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (source == "CRAN") {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package installation failed: ", pkg)
  }
}

ensure_package("data.table", "CRAN")
ensure_package("limma", "BIOC")

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

parse_sample_name <- function(sample_name) {
  parts <- strsplit(sample_name, "-", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Could not parse sample name: ", sample_name)
  }

  suffix <- parts[[2]]
  data.frame(
    sample = sample_name,
    condition = parts[[1]],
    batch = sub("[0-9]+$", "", suffix),
    replicate = sub("^[A-Za-z]+", "", suffix),
    stringsAsFactors = FALSE
  )
}

extract_gene_symbol <- function(gene_value, fallback_value) {
  raw_value <- trimws(ifelse(is.na(gene_value), "", gene_value))
  if (!nzchar(raw_value)) {
    raw_value <- trimws(ifelse(is.na(fallback_value), "", fallback_value))
  }

  tokens <- unlist(strsplit(raw_value, "[;,| ]+"))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0) {
    return(NA_character_)
  }

  tokens[[1]]
}

median_center_normalize <- function(mat) {
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  global_median <- median(sample_medians, na.rm = TRUE)
  sweep(mat, 2, sample_medians - global_median, "-")
}

impute_for_pca <- function(mat) {
  imputed <- mat
  for (idx in seq_len(nrow(imputed))) {
    row_values <- imputed[idx, ]
    if (all(is.na(row_values))) {
      next
    }
    row_median <- median(row_values, na.rm = TRUE)
    row_values[is.na(row_values)] <- row_median
    imputed[idx, ] <- row_values
  }
  imputed
}

write_pca_plot <- function(mat, sample_info, output_file, title_text) {
  complete_rows <- rowSums(is.na(mat)) < ncol(mat)
  plot_mat <- mat[complete_rows, , drop = FALSE]
  if (nrow(plot_mat) < 2) {
    return(invisible(NULL))
  }

  plot_mat <- impute_for_pca(plot_mat)
  row_sd <- apply(plot_mat, 1, sd)
  variable_rows <- is.finite(row_sd) & row_sd > 0
  plot_mat <- plot_mat[variable_rows, , drop = FALSE]
  if (nrow(plot_mat) < 2) {
    message("Skipping PCA plot because fewer than two variable proteins remained after filtering.")
    return(invisible(NULL))
  }

  pca <- prcomp(t(plot_mat), center = TRUE, scale. = TRUE)
  pc_var <- summary(pca)$importance[2, 1:2] * 100

  pdf(output_file, width = 7, height = 6)
  on.exit(dev.off(), add = TRUE)

  colors <- as.integer(factor(sample_info$condition))
  symbols <- as.integer(factor(sample_info$batch)) + 14

  plot(
    pca$x[, 1],
    pca$x[, 2],
    col = colors,
    pch = symbols,
    cex = 1.4,
    xlab = sprintf("PC1 (%.1f%%)", pc_var[[1]]),
    ylab = sprintf("PC2 (%.1f%%)", pc_var[[2]]),
    main = title_text
  )
  text(pca$x[, 1], pca$x[, 2], labels = sample_info$sample, pos = 3, cex = 0.75)
  legend(
    "topright",
    legend = unique(sample_info$condition),
    col = seq_along(unique(sample_info$condition)),
    pch = 16,
    bty = "n",
    title = "Condition"
  )
  legend(
    "bottomright",
    legend = unique(sample_info$batch),
    pch = seq_along(unique(sample_info$batch)) + 14,
    bty = "n",
    title = "Batch"
  )
}

write_boxplot <- function(mat, output_file, title_text) {
  pdf(output_file, width = 9, height = 5)
  on.exit(dev.off(), add = TRUE)
  boxplot(
    as.data.frame(mat),
    las = 2,
    col = "grey85",
    outline = FALSE,
    main = title_text,
    ylab = "log2 intensity"
  )
}

build_contrasts <- function(condition_levels) {
  if (length(condition_levels) < 2) {
    stop("At least two conditions are required for differential analysis.")
  }

  contrast_strings <- character()
  contrast_names <- character()
  pairs <- combn(condition_levels, 2, simplify = FALSE)
  for (pair in pairs) {
    contrast_names <- c(contrast_names, paste0(pair[[2]], "_vs_", pair[[1]]))
    contrast_strings <- c(contrast_strings, paste0("condition", pair[[2]], "-condition", pair[[1]]))
  }

  names(contrast_strings) <- contrast_names
  contrast_strings
}

cat("============================================\n")
cat(" Step 4: Proteome differential analysis\n")
cat(" Input  : ", normalizePath(input_file, mustWork = FALSE), "\n", sep = "")
cat(" Output : ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
cat("============================================\n")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

proteome_df <- fread(input_file, data.table = FALSE)
sample_cols <- names(proteome_df)[grepl("^[^-]+-[A-Za-z]+[0-9]+$", names(proteome_df))]

if (length(sample_cols) < 4) {
  stop("Could not identify proteome sample columns from the input table.")
}

sample_info <- do.call(rbind, lapply(sample_cols, parse_sample_name))
sample_info$condition <- factor(sample_info$condition)
sample_info$batch <- factor(sample_info$batch)

fwrite(sample_info, file.path(output_dir, "sample_metadata.tsv"), sep = "\t")

annotation_cols <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
present_annotation_cols <- intersect(annotation_cols, names(proteome_df))

intensity_df <- proteome_df[, sample_cols, drop = FALSE]
intensity_df[] <- lapply(intensity_df, function(x) as.numeric(as.character(x)))
intensity_mat <- as.matrix(intensity_df)
rownames(intensity_mat) <- seq_len(nrow(intensity_mat))

positive_values <- intensity_mat[intensity_mat > 0 & !is.na(intensity_mat)]
if (length(positive_values) == 0) {
  stop("No positive quantitative values were found in the proteome matrix.")
}

pseudocount <- min(positive_values) / 2
intensity_mat[intensity_mat <= 0] <- NA_real_
log2_mat <- log2(intensity_mat + pseudocount)
norm_mat <- median_center_normalize(log2_mat)
colnames(norm_mat) <- sample_cols

write_boxplot(log2_mat, file.path(output_dir, "sample_boxplot_before_normalization.pdf"), "Proteome intensities before normalization")
write_boxplot(norm_mat, file.path(output_dir, "sample_boxplot_after_normalization.pdf"), "Proteome intensities after median-centering")
write_pca_plot(norm_mat, sample_info, file.path(output_dir, "sample_pca_after_normalization.pdf"), "PCA after normalization")

gene_values <- if ("Genes" %in% names(proteome_df)) {
  proteome_df$Genes
} else {
  rep(NA_character_, nrow(proteome_df))
}

protein_group_values <- if ("Protein.Group" %in% names(proteome_df)) {
  proteome_df$Protein.Group
} else {
  rep(NA_character_, nrow(proteome_df))
}

gene_symbol <- mapply(
  extract_gene_symbol,
  gene_values,
  protein_group_values,
  USE.NAMES = FALSE
)

annotated_df <- cbind(
  proteome_df[, present_annotation_cols, drop = FALSE],
  data.frame(gene_symbol = gene_symbol, stringsAsFactors = FALSE)
)

min_reps_per_group <- 2
keep_per_condition <- sapply(levels(sample_info$condition), function(cond) {
  cond_cols <- sample_info$sample[sample_info$condition == cond]
  rowSums(!is.na(norm_mat[, cond_cols, drop = FALSE])) >= min_reps_per_group
})

if (is.null(dim(keep_per_condition))) {
  keep_per_condition <- matrix(keep_per_condition, ncol = 1)
}

keep_rows <- rowSums(keep_per_condition) == length(levels(sample_info$condition))
keep_rows <- keep_rows & !is.na(gene_symbol) & nzchar(gene_symbol)

filtered_mat <- norm_mat[keep_rows, , drop = FALSE]
filtered_annot <- annotated_df[keep_rows, , drop = FALSE]

if (nrow(filtered_mat) == 0) {
  stop("No proteins remained after filtering. Check gene symbols and missing-value thresholds.")
}

row_missing <- rowSums(is.na(filtered_mat))
row_mean <- rowMeans(filtered_mat, na.rm = TRUE)
priority_df <- data.frame(
  row_index = seq_len(nrow(filtered_mat)),
  gene_symbol = filtered_annot$gene_symbol,
  missing_n = row_missing,
  mean_intensity = row_mean,
  stringsAsFactors = FALSE
)
priority_df <- priority_df[order(priority_df$gene_symbol, priority_df$missing_n, -priority_df$mean_intensity), ]
selected_idx <- priority_df$row_index[!duplicated(priority_df$gene_symbol)]

collapsed_mat <- filtered_mat[selected_idx, , drop = FALSE]
collapsed_annot <- filtered_annot[selected_idx, , drop = FALSE]
rownames(collapsed_mat) <- collapsed_annot$gene_symbol

if (nrow(collapsed_mat) == 0) {
  stop("No unique gene symbols remained after collapsing duplicated proteins.")
}

fwrite(
  cbind(collapsed_annot, as.data.frame(collapsed_mat, check.names = FALSE)),
  file.path(output_dir, "normalized_proteome_matrix.tsv"),
  sep = "\t"
)

design <- model.matrix(~ 0 + condition + batch, data = sample_info)
if (qr(design)$rank < ncol(design)) {
  design <- model.matrix(~ 0 + condition, data = sample_info)
}

contrast_strings <- build_contrasts(levels(sample_info$condition))
contrast_matrix <- makeContrasts(contrasts = unname(contrast_strings), levels = design)
colnames(contrast_matrix) <- names(contrast_strings)

fit <- lmFit(collapsed_mat, design, na.action = na.exclude)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)

summary_lines <- c(
  "Proteome pseudo-GSEA preprocessing summary",
  paste0("Input file: ", normalizePath(input_file, mustWork = FALSE)),
  paste0("Detected samples: ", paste(sample_info$sample, collapse = ", ")),
  paste0("Conditions: ", paste(levels(sample_info$condition), collapse = ", ")),
  paste0("Batches: ", paste(levels(sample_info$batch), collapse = ", ")),
  paste0("Pseudocount for log2 transform: ", signif(pseudocount, 4)),
  paste0("Rows before filtering: ", nrow(proteome_df)),
  paste0("Rows after filtering: ", nrow(filtered_mat)),
  paste0("Unique genes after collapsing duplicates: ", nrow(collapsed_mat)),
  "Bias mitigation:",
  "- Median-centering was applied across samples to reduce loading and global intensity bias.",
  "- EXP/REQ was modeled as a batch factor when the design matrix was full rank.",
  "- No global missing-value imputation was used for differential testing.",
  "- Ranking for downstream GSEA is provided as the moderated t statistic."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

for (contrast_name in colnames(contrast_matrix)) {
  result_df <- topTable(fit, coef = contrast_name, number = Inf, sort.by = "none")
  result_df$geneID <- rownames(result_df)
  result_df$rank_metric <- result_df$t
  result_df <- result_df[, c("geneID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "rank_metric")]
  names(result_df)[names(result_df) == "logFC"] <- "log2FoldChange"
  names(result_df)[names(result_df) == "AveExpr"] <- "baseMean"
  names(result_df)[names(result_df) == "P.Value"] <- "pvalue"
  names(result_df)[names(result_df) == "adj.P.Val"] <- "padj"

  fwrite(
    result_df,
    file.path(output_dir, paste0(contrast_name, "_DE_Results.tsv")),
    sep = "\t"
  )
}

cat("Differential analysis complete.\n")
