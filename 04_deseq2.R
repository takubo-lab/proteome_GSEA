#!/usr/bin/env Rscript
# =============================================================================
# Step 4: DESeq2 Differential Expression Analysis
# =============================================================================
# Performs all-group pairwise comparisons using DESeq2.
# Reads configuration from config.sh and sample info from samples.tsv.
#
# Cross-platform:
#   Linux : bash run_pipeline.sh (full pipeline)
#   Win/Mac: Rscript scripts/04_deseq2.R [config.sh]
#            (count.txt と samples.tsv が事前に用意されている前提)
#
# 可視化設定: plot_config.default.yml / plot_config.yml
# =============================================================================

# --- Parse config ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.sh")

parse_config <- function(config_path) {
  lines <- readLines(config_path)
  lines <- lines[!grepl("^\\s*#", lines) & !grepl("^\\s*$", lines)]
  lines <- lines[!grepl("^(if|elif|else|fi|source|set|SCRIPT_DIR|PROJECT_DIR)", lines)]
  env <- list()
  for (line in lines) {
    line <- sub("\\s+#.*$", "", line)
    m <- regmatches(line, regexpr("^([A-Za-z_][A-Za-z0-9_]*)=[\"']?([^\"']*)[\"']?", line))
    if (length(m) == 1 && nchar(m) > 0) {
      parts <- strsplit(m, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      val <- paste(parts[-1], collapse = "=")
      val <- gsub('^["\047]|["\047]$', '', val)
      env[[key]] <- val
    }
  }
  env
}

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_sample_name <- function(x) {
  gsub("-", "_", x)
}

resolve_input_path <- function(path_value, project_dir, cfg = list()) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return("")
  }

  resolved <- trimws(path_value)
  replacements <- c(
    PROJECT_DIR = project_dir,
    PIPELINE_DIR = cfg$PIPELINE_DIR %||% project_dir,
    HOME = Sys.getenv("HOME", unset = "")
  )

  for (nm in names(replacements)) {
    resolved <- gsub(paste0("\\$\\{", nm, "\\}"), replacements[[nm]], resolved)
    resolved <- gsub(paste0("\\$", nm), replacements[[nm]], resolved)
  }

  path.expand(resolved)
}

parse_path_list <- function(x, project_dir, cfg = list()) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }

  unique(vapply(strsplit(x, ",", fixed = TRUE)[[1]], function(p) {
    resolve_input_path(p, project_dir, cfg)
  }, character(1)))
}

read_sample_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  required_cols <- c("sample_name", "group", "fq_prefix", "lane_suffix")

  if (length(lines) == 0) {
    stop("ERROR: samples.tsv is empty: ", path)
  }

  normalized <- sub("^\\s*#\\s*", "", lines)
  header_idx <- which(vapply(normalized, function(line) {
    fields <- strsplit(trimws(line), "\\t")[[1]]
    length(fields) >= length(required_cols) && identical(fields[seq_along(required_cols)], required_cols)
  }, logical(1)))[1]

  if (!is.na(header_idx) && grepl("^\\s*#", lines[header_idx])) {
    lines[header_idx] <- normalized[header_idx]
  }

  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^\\s*#", lines)]

  dat <- fread(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    header = TRUE,
    blank.lines.skip = TRUE
  )

  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("ERROR: Missing required columns in samples.tsv: ",
         paste(missing_cols, collapse = ", "),
         " (file: ", path, ")")
  }

  dat[, c("sample_name", "group", "fq_prefix", "lane_suffix"), with = FALSE]
}

read_sample_tables <- function(sample_files) {
  sample_list <- lapply(sample_files, function(path) {
    if (!file.exists(path)) {
      stop("ERROR: samples.tsv not found: ", path)
    }
    read_sample_table(path)
  })

  sample_info <- dplyr::bind_rows(sample_list)
  dup <- sample_info$sample_name[duplicated(sample_info$sample_name)]
  if (length(dup) > 0) {
    stop("ERROR: Duplicate sample names found across merged sample sheets: ",
         paste(unique(dup), collapse = ", "))
  }
  sample_info
}

read_count_tables <- function(count_files) {
  count_list <- lapply(count_files, function(path) {
    if (!file.exists(path)) {
      stop("ERROR: Count file not found: ", path)
    }
    fread(path, header = TRUE)
  })

  geneid_col <- colnames(count_list[[1]])[1]
  merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = geneid_col), count_list)
  merged[is.na(merged)] <- 0

  dup_cols <- colnames(merged)[duplicated(colnames(merged))]
  if (length(dup_cols) > 0) {
    stop("ERROR: Duplicate count columns found across merged count files: ",
         paste(unique(dup_cols), collapse = ", "))
  }

  merged
}

cfg <- parse_config(config_file)

# Resolve project directory
project_dir <- getwd()
output_dir  <- file.path(project_dir, "DESeq2")

merge_samples_env <- Sys.getenv("PIPELINE_MERGE_SAMPLES", unset = "")
merge_counts_env  <- Sys.getenv("PIPELINE_MERGE_COUNTS", unset = "")

sample_files <- parse_path_list(merge_samples_env, project_dir, cfg)
count_files  <- parse_path_list(merge_counts_env, project_dir, cfg)

if (length(sample_files) == 0) {
  sample_files <- file.path(project_dir, "samples.tsv")
}
if (length(count_files) == 0) {
  count_files <- file.path(project_dir, "counts", "count.txt")
}

samples_tsv <- paste(sample_files, collapse = ", ")
count_file  <- paste(count_files, collapse = ", ")

min_count   <- as.integer(cfg$MIN_COUNT_FILTER %||% 50)
pca_top     <- as.integer(cfg$PCA_TOP_GENES %||% 500)
padj_thr    <- as.numeric(cfg$VOLCANO_PADJ %||% 0.001)
lfc_thr     <- as.numeric(cfg$VOLCANO_LOG2FC %||% 1)
gene_list   <- unlist(strsplit(cfg$HIGHLIGHT_GENES %||% "", ","))

# --- Load packages ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(umap)
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(ggrepel)
})

# --- Load plot config ---
# Determine script directory for sourcing plot_utils.R
script_dir <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE)
}, error = function(e) {
  candidates <- c(
    file.path(project_dir, "scripts"),
    file.path(cfg$PIPELINE_DIR %||% ".", "scripts")
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "plot_utils.R"))) return(normalizePath(d))
  }
  "scripts"
})
source(file.path(script_dir, "plot_utils.R"))
pcfg <- load_plot_config()

set.seed(123)

cat("============================================\n")
cat(" Step 4: DESeq2 Analysis\n")
cat(" Count file : ", count_file, "\n")
cat(" Samples    : ", samples_tsv, "\n")
cat(" Output     : ", output_dir, "\n")
cat("============================================\n")

# --- Read sample info ---
sample_info <- read_sample_tables(sample_files)

# --- Read count matrix ---
raw <- read_count_tables(count_files)
geneid_col <- colnames(raw)[1]
rawCounts  <- raw %>%
  tibble::column_to_rownames(geneid_col)

# Clean column names (remove BAM path prefixes if present)
clean_names <- colnames(rawCounts)
clean_names <- gsub(".*[/\\\\]", "", clean_names)
clean_names <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", clean_names)
colnames(rawCounts) <- clean_names

# Match samples to columns
count_keys <- normalize_sample_name(clean_names)
sample_keys <- normalize_sample_name(sample_info$sample_name)
matched_idx <- match(sample_keys, count_keys)

if (all(is.na(matched_idx))) {
  count_preview <- clean_names[seq_len(min(5, length(clean_names)))]
  sample_preview <- sample_info$sample_name[seq_len(min(5, nrow(sample_info)))]
  stop("ERROR: No sample names in samples.tsv match count matrix columns.\n",
       "  Count columns: ", paste(count_preview, collapse = ", "), "\n",
       "  Sample names : ", paste(sample_preview, collapse = ", "))
}

sample_info <- sample_info[!is.na(matched_idx), ]
matched_idx <- matched_idx[!is.na(matched_idx)]

rawCounts <- rawCounts[, matched_idx, drop = FALSE]
colnames(rawCounts) <- sample_info$sample_name

samples <- factor(sample_info$sample_name, levels = sample_info$sample_name)
groups  <- factor(sample_info$group)

sampleData <- data.frame(
  samples = samples,
  class   = groups,
  row.names = as.character(samples)
)

# --- Filter low-count genes ---
rawCounts_filt <- rawCounts[rowSums(rawCounts) > min_count, ]
rawCounts_filt <- rawCounts_filt[!duplicated(rawCounts_filt), ]

cat("Genes after filtering: ", nrow(rawCounts_filt), "\n")

# --- Create DESeq2 object ---
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dds <- DESeqDataSetFromMatrix(
  countData = rawCounts_filt,
  colData   = sampleData,
  design    = ~class
)

# Normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
norm_df <- data.frame(normalized_counts) %>%
  rownames_to_column("GeneName") %>%
  tibble()
fwrite(norm_df, file.path(output_dir, "Counts.normalized.txt"), sep = "\t")

# --- Run DESeq2 ---
dds <- DESeq(dds)

# --- PCA ---
vst_out  <- vst(dds)
stdev    <- apply(assay(vst_out), 1, sd)
top_n    <- min(pca_top, nrow(assay(vst_out)))
vst_short <- assay(vst_out)[rev(order(stdev))[1:top_n], ]

pcaData    <- data.frame(prcomp(t(vst_short))$x)
percentVar <- round(100 * summary(prcomp(t(vst_short)))$importance[2, ])

group_cols <- get_group_colors(as.character(sampleData$class), pcfg)

p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, fill = sampleData$class)) +
  geom_point(size = pcfg$pca$point_size, alpha = pcfg$pca$point_alpha,
             pch = 21, color = "Black") +
  scale_fill_manual(values = group_cols) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(fill = "Group") +
  ggtitle("PCA") +
  theme_pipeline(pcfg) +
  theme(axis.text = element_text(size = pcfg$theme$base_size + 3),
        axis.title = element_text(size = pcfg$theme$base_size + 3))

if (pcfg$pca$show_labels) {
  p_pca <- p_pca + geom_text_repel(aes(label = rownames(pcaData)),
                                    size = pcfg$pca$label_size, max.overlaps = 20)
}

save_plot(p_pca, "PCA_class", type = "pca", outdir = output_dir, cfg = pcfg)

# --- UMAP ---
config_umap <- umap.defaults
config_umap$n_neighbors <- min(pcfg$umap$n_neighbors, ncol(vst_short) - 1)
u <- umap(t(vst_short), config = config_umap)
u_df <- data.frame(u$layout)

p_umap <- ggplot(u_df, aes(x = X1, y = X2, fill = sampleData$class)) +
  geom_point(size = pcfg$umap$point_size, alpha = pcfg$umap$point_alpha,
             pch = 21, color = "Black") +
  scale_fill_manual(values = group_cols) +
  xlab("UMAP_1") + ylab("UMAP_2") +
  labs(fill = "Group") +
  ggtitle("UMAP") +
  theme_pipeline(pcfg) +
  theme(axis.text = element_text(size = pcfg$theme$base_size + 3),
        axis.title = element_text(size = pcfg$theme$base_size + 3))

if (pcfg$umap$show_labels) {
  p_umap <- p_umap + geom_text_repel(aes(label = rownames(u_df)),
                                      size = pcfg$umap$label_size, max.overlaps = 20)
}

save_plot(p_umap, "UMAP_all", type = "umap", outdir = output_dir, cfg = pcfg)

# --- PC Loadings ---
loadings <- prcomp(t(vst_short))$rotation
loadings_df <- data.frame(PC1 = loadings[, 1], PC2 = loadings[, 2]) %>%
  rownames_to_column("Gene")
fwrite(loadings_df, file.path(output_dir, "PC_loadings.tsv"), sep = "\t")

# --- Pairwise comparisons (all groups) ---
all_levels <- levels(groups)
all_pairs  <- combn(all_levels, 2, simplify = FALSE)

cat("Performing ", length(all_pairs), " pairwise comparisons...\n")

# Volcano plot function (uses plot_config)
plot_volcano <- function(deseq2Results, gene_list = c(""), file_out,
                         padj_thr = 0.001, lfc_thr = 1, pcfg = NULL) {
  if (is.null(pcfg)) pcfg <- load_plot_config()
  vcfg <- pcfg$volcano
  vcols <- get_volcano_colors(pcfg)

  set.seed(42)
  df <- data.frame(deseq2Results) %>%
    mutate(sig = case_when(
      geneID %in% gene_list ~ "marked",
      padj < padj_thr & log2FoldChange >  lfc_thr & !(geneID %in% gene_list) ~ "up",
      padj < padj_thr & log2FoldChange < -lfc_thr & !(geneID %in% gene_list) ~ "down",
      TRUE ~ "non"
    )) %>%
    mutate(mark = geneID %in% gene_list) %>%
    arrange(mark)

  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), fill = sig)) +
    geom_point(pch = 21, aes(size = mark), alpha = vcfg$alpha_ns) +
    scale_fill_manual(values = c(
      down = unname(vcols["down"]),
      marked = unname(vcols["marked"]),
      non = unname(vcols["ns"]),
      up = unname(vcols["up"])
    )) +
    scale_size_manual(values = c(vcfg$point_size, vcfg$marked_point_size)) +
    geom_text_repel(data = filter(df, mark & log2FoldChange > 0),
                    aes(label = geneID), size = vcfg$label_size, max.overlaps = 10,
                    nudge_y = vcfg$nudge_y, nudge_x = vcfg$nudge_x_pos,
                    box.padding = 0.5, color = vcfg$label_color) +
    geom_text_repel(data = filter(df, mark & log2FoldChange < 0),
                    aes(label = geneID), size = vcfg$label_size, max.overlaps = 10,
                    nudge_y = vcfg$nudge_y, nudge_x = vcfg$nudge_x_neg,
                    box.padding = 0.5, color = vcfg$label_color) +
    geom_vline(xintercept = lfc_thr, linetype = "dashed") +
    geom_vline(xintercept = -lfc_thr, linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    theme_pipeline(pcfg)

  save_plot(p, paste0(basename(file_out), "_volcano"), type = "volcano",
            outdir = dirname(file_out), cfg = pcfg)
}

# Paired test function (for paired experimental designs)
my_paired_test <- function(norm_counts, sampleData, group1, group2) {
  sampleData <- as.data.frame(sampleData)
  colnames(sampleData) <- c("samples", "class")
  sampleData$patient <- gsub(paste0("_", group1, "$|_", group2, "$"), "",
                             as.character(sampleData$samples))

  df_g1 <- sampleData %>% filter(class == group1) %>% arrange(patient)
  df_g2 <- sampleData %>% filter(class == group2) %>% arrange(patient)
  common_patients <- intersect(df_g1$patient, df_g2$patient)

  if (length(common_patients) < 2) {
    message("  Not enough paired samples, skipping paired test.")
    return(NULL)
  }

  df_g1 <- df_g1 %>% filter(patient %in% common_patients) %>% arrange(patient)
  df_g2 <- df_g2 %>% filter(patient %in% common_patients) %>% arrange(patient)
  message(sprintf("  Paired test: %d patients (%s vs %s)",
                  length(common_patients), group1, group2))

  mat_g1 <- log2(norm_counts[, as.character(df_g1$samples), drop = FALSE] + 1)
  mat_g2 <- log2(norm_counts[, as.character(df_g2$samples), drop = FALSE] + 1)

  n_genes     <- nrow(norm_counts)
  ttest_p     <- numeric(n_genes)
  ttest_stat  <- numeric(n_genes)
  wilcox_p    <- numeric(n_genes)
  wilcox_stat <- numeric(n_genes)
  mean_log2fc <- numeric(n_genes)

  for (i in seq_len(n_genes)) {
    x_g1 <- as.numeric(mat_g1[i, ])
    x_g2 <- as.numeric(mat_g2[i, ])
    mean_log2fc[i] <- mean(x_g1 - x_g2, na.rm = TRUE)
    t_res <- tryCatch(t.test(x_g1, x_g2, paired = TRUE),
                      error = function(e) list(statistic = NA_real_, p.value = NA_real_))
    ttest_p[i]    <- t_res$p.value
    ttest_stat[i] <- as.numeric(t_res$statistic)
    w_res <- tryCatch(wilcox.test(x_g1, x_g2, paired = TRUE, exact = FALSE),
                      error = function(e) list(statistic = NA_real_, p.value = NA_real_))
    wilcox_p[i]    <- w_res$p.value
    wilcox_stat[i] <- as.numeric(w_res$statistic)
  }

  data.frame(
    geneID            = rownames(norm_counts),
    paired_mean_log2FC = mean_log2fc,
    ttest_stat        = ttest_stat,
    ttest_p           = ttest_p,
    ttest_padj_BH     = p.adjust(ttest_p, method = "BH"),
    wilcox_stat       = wilcox_stat,
    wilcox_p          = wilcox_p,
    wilcox_padj_BH    = p.adjust(wilcox_p, method = "BH"),
    stringsAsFactors  = FALSE
  )
}

# Run all pairwise comparisons
for (pair in all_pairs) {
  g1 <- pair[1]
  g2 <- pair[2]
  comp_name <- paste0(g1, "_vs_", g2)
  message(sprintf("[DESeq2] %s", comp_name))

  res <- results(dds, contrast = c("class", g1, g2))
  res_df <- as_tibble(res, rownames = "geneID") %>% arrange(padj)
  fwrite(res_df, file.path(output_dir, paste0(comp_name, "_DE_Results.tsv")), sep = "\t")

  # Paired test (if applicable)
  norm_sub <- counts(dds, normalized = TRUE)
  sd_sub   <- sampleData[sampleData$class %in% c(g1, g2), , drop = FALSE]
  norm_sub <- norm_sub[, as.character(sd_sub$samples), drop = FALSE]

  paired_res <- my_paired_test(norm_sub, sd_sub, group1 = g1, group2 = g2)
  if (!is.null(paired_res)) {
    merged <- left_join(as.data.frame(res_df), paired_res, by = "geneID") %>%
      arrange(padj)
    fwrite(merged, file.path(output_dir, paste0(comp_name, "_DE_Results_with_paired.tsv")),
           sep = "\t")
  }

  # Volcano
  plot_volcano(res_df, gene_list,
               file.path(output_dir, comp_name),
               padj_thr = padj_thr, lfc_thr = lfc_thr, pcfg = pcfg)
}

cat("============================================\n")
cat(" Step 4: DESeq2 Analysis Complete\n")
cat("============================================\n")
