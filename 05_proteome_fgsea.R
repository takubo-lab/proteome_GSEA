#!/usr/bin/env Rscript

if (!exists("args", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = TRUE)
}
de_dir <- if (length(args) >= 1) args[[1]] else "Proteome_DE"
gmt_file <- if (length(args) >= 2) args[[2]] else file.path("Gene_set", "Human_old_HSC_set2.gmt")
output_dir <- if (length(args) >= 3) args[[3]] else "Proteome_fGSEA"

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
ensure_package("fgsea", "BIOC")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  gmt_list <- vector("list", length(lines))
  gmt_names <- character(length(lines))

  for (idx in seq_along(lines)) {
    fields <- strsplit(lines[[idx]], "\t", fixed = FALSE)[[1]]
    gmt_names[[idx]] <- fields[[1]]
    genes <- unique(toupper(fields[-c(1, 2)]))
    genes <- genes[nzchar(genes)]
    gmt_list[[idx]] <- genes
  }

  names(gmt_list) <- gmt_names
  gmt_list
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

write_barplot <- function(fgsea_df, output_file, title_text) {
  sig_df <- fgsea_df[!is.na(fgsea_df$padj) & fgsea_df$padj < 0.05, , drop = FALSE]
  if (nrow(sig_df) == 0) {
    return(invisible(NULL))
  }

  up_df <- sig_df[sig_df$NES > 0, , drop = FALSE]
  down_df <- sig_df[sig_df$NES < 0, , drop = FALSE]
  up_df <- up_df[order(up_df$padj), , drop = FALSE]
  down_df <- down_df[order(down_df$padj), , drop = FALSE]
  plot_df <- rbind(head(down_df, 10), head(up_df, 10))
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- plot_df[order(plot_df$NES), , drop = FALSE]
  colors <- ifelse(plot_df$NES > 0, "firebrick3", "steelblue4")

  pdf(output_file, width = 9, height = 7)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(5, 12, 4, 2))
  barplot(
    plot_df$NES,
    names.arg = plot_df$pathway,
    horiz = TRUE,
    las = 1,
    col = colors,
    border = NA,
    main = title_text,
    xlab = "Normalized Enrichment Score"
  )
}

cat("============================================\n")
cat(" Step 5: Proteome fGSEA\n")
cat(" DE dir : ", normalizePath(de_dir, mustWork = FALSE), "\n", sep = "")
cat(" GMT    : ", normalizePath(gmt_file, mustWork = FALSE), "\n", sep = "")
cat(" Output : ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
cat("============================================\n")

if (!dir.exists(de_dir)) {
  stop("Differential analysis directory not found: ", de_dir)
}

if (!file.exists(gmt_file)) {
  stop("GMT file not found: ", gmt_file)
}

pathways <- read_gmt(gmt_file)
de_files <- list.files(de_dir, pattern = "_DE_Results\\.tsv$", full.names = TRUE)

if (length(de_files) == 0) {
  stop("No differential result files were found in: ", de_dir)
}

for (de_file in de_files) {
  de_df <- fread(de_file, data.table = FALSE)
  if (!("geneID" %in% names(de_df))) {
    stop("geneID column is missing in: ", de_file)
  }

  rank_col <- if ("rank_metric" %in% names(de_df)) {
    "rank_metric"
  } else if ("t" %in% names(de_df)) {
    "t"
  } else {
    "log2FoldChange"
  }

  stats <- de_df[[rank_col]]
  names(stats) <- toupper(de_df$geneID)
  keep_idx <- !is.na(stats) & nzchar(names(stats))
  stats <- stats[keep_idx]
  stats <- tapply(stats, names(stats), max)
  stats <- sort(unlist(stats), decreasing = TRUE)

  fgsea_res <- fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )

  overlap_n <- vapply(fgsea_res$pathway, function(pathway_name) {
    length(intersect(pathways[[pathway_name]], names(stats)))
  }, integer(1))

  fgsea_df <- as.data.frame(fgsea_res)
  fgsea_df$overlap_n <- overlap_n
  fgsea_df$leadingEdge <- vapply(fgsea_df$leadingEdge, function(x) paste(x, collapse = ","), character(1))
  fgsea_df <- fgsea_df[order(fgsea_df$padj, -abs(fgsea_df$NES)), ]

  comparison_name <- sub("_DE_Results\\.tsv$", "", basename(de_file))
  comparison_dir <- file.path(output_dir, paste0("fGSEA_", sanitize_filename(comparison_name)))
  dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

  fwrite(
    fgsea_df,
    file.path(comparison_dir, paste0("fgsea_Results_", comparison_name, ".tsv")),
    sep = "\t"
  )

  write_barplot(
    fgsea_df,
    file.path(comparison_dir, paste0("fgsea_barplot_", comparison_name, ".pdf")),
    paste0("Proteome pseudo-GSEA: ", comparison_name)
  )
}

writeLines(
  c(
    "Proteome pseudo-GSEA notes",
    paste0("Gene set file: ", normalizePath(gmt_file, mustWork = FALSE)),
    "Gene symbols were harmonized by uppercasing both the proteome-derived ranks and GMT members.",
    "This is a practical approximation for mouse-to-human symbol matching, not a full ortholog conversion.",
    "Interpret pathways with low overlap_n cautiously."
  ),
  file.path(output_dir, "fgsea_summary.txt")
)

cat("Proteome fGSEA complete.\n")
