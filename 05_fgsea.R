#!/usr/bin/env Rscript
# =============================================================================
# Step 5: fGSEA (fast Gene Set Enrichment Analysis)
# =============================================================================
# Runs fGSEA on DESeq2 results for MSigDB gene sets + custom gene sets.
#
# Cross-platform:
#   Linux : bash run_pipeline.sh (full pipeline)
#   Win/Mac: Rscript scripts/05_fgsea.R [config.sh]
#
# 可視化設定: plot_config.default.yml / plot_config.yml
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config.sh")

# --- Parse config ---
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

resolve_config_path <- function(path_value, project_dir, cfg = list()) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return("")
  }

  resolved <- path_value
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

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

cfg <- parse_config(config_file)

project_dir  <- getwd()
deseq2_dir   <- file.path(project_dir, "DESeq2")
fgsea_dir    <- file.path(project_dir, "fGSEA")
genome       <- cfg$GENOME %||% "hg38"

fgsea_min    <- as.integer(cfg$FGSEA_MIN_SIZE %||% 15)
fgsea_max    <- as.integer(cfg$FGSEA_MAX_SIZE %||% 500)
fgsea_nperm  <- as.integer(cfg$FGSEA_NPERM %||% 10000)
fgsea_nproc  <- as.integer(cfg$FGSEA_NPROC %||% 12)
fgsea_padj   <- as.numeric(cfg$FGSEA_PADJ %||% 0.05)
custom_gmt   <- resolve_config_path(cfg$CUSTOM_GMT %||% "", project_dir, cfg)

species <- ifelse(genome == "hg38", "Homo sapiens", "Mus musculus")

# --- Load packages ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(magrittr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# --- Load plot config ---
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
dir.create(fgsea_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Step 5: fGSEA\n")
cat(" Species : ", species, "\n")
cat(" DESeq2  : ", deseq2_dir, "\n")
cat(" Output  : ", fgsea_dir, "\n")
cat("============================================\n")

# --- Build gene set collection ---
msig_categories <- c("H", "C2", "C3", "C4", "C5", "C6", "C7")
msigdb_list <- lapply(msig_categories, function(cat) {
  msigdb <- msigdbr::msigdbr(species = species, category = cat)
  split(toupper(msigdb$gene_symbol), msigdb$gs_name) %>%
    lapply(function(x) unique(x[nzchar(x)]))
})
names(msigdb_list) <- msig_categories

# Add custom gene set if available
if (nchar(custom_gmt) > 0 && file.exists(custom_gmt)) {
  cat("Loading custom gene set: ", custom_gmt, "\n")
  gmt_lines <- readLines(custom_gmt, warn = FALSE)
  gmt_lines <- gmt_lines[nzchar(trimws(gmt_lines))]
  gmt_split <- strsplit(gmt_lines, "\t", fixed = TRUE)
  gmt_list <- setNames(vector("list", length(gmt_split)), vapply(gmt_split, `[[`, character(1), 1))
  for (i in seq_along(gmt_split)) {
    entries <- gmt_split[[i]]
    genes <- unique(toupper(entries[-c(1, 2)]))
    genes <- genes[nzchar(genes)]
    gmt_list[[i]] <- genes
  }
  custom_name <- tools::file_path_sans_ext(basename(custom_gmt))
  msigdb_list[[custom_name]] <- gmt_list
  cat("  Added ", length(gmt_list), " custom sets from ", custom_name, "\n")
} else if (nchar(custom_gmt) > 0) {
  warning("Custom gene set file not found: ", custom_gmt)
}

save_enrichment_plots <- function(fgsea_res, pathways, ranks, comp_dir, comparison_name, category_name) {
  if (!(category_name == "H" || !(category_name %in% msig_categories))) {
    return(invisible(NULL))
  }

  plot_dir <- file.path(comp_dir, paste0("enrichment_", sanitize_filename(category_name)))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  targets <- fgsea_res[0, , drop = FALSE]
  if (category_name == "H") {
    targets <- fgsea_res %>% filter(!is.na(padj) & padj < fgsea_padj) %>% arrange(padj)
  } else {
    targets <- fgsea_res %>% filter(!is.na(NES)) %>% arrange(padj)
  }

  if (nrow(targets) == 0) {
    return(invisible(NULL))
  }

  for (pathway_name in targets$pathway) {
    genes <- pathways[[pathway_name]]
    if (is.null(genes) || length(genes) == 0) next

    p <- fgsea::plotEnrichment(genes, ranks) +
      labs(
        title = paste0("Enrichment: ", pathway_name),
        subtitle = paste0(comparison_name, " — ", category_name),
        x = "Rank in ordered gene list",
        y = "Enrichment score"
      ) +
      theme_pipeline(pcfg)

    save_plot(
      p,
      paste0(
        "fgsea_enrichment_",
        sanitize_filename(comparison_name),
        "_",
        sanitize_filename(category_name),
        "_",
        sanitize_filename(pathway_name)
      ),
      type = "fgsea_bar",
      outdir = plot_dir,
      cfg = pcfg
    )
  }
}

# --- fGSEA function ---
run_fgsea <- function(de_file, gene_set_list) {
  de <- fread(de_file, sep = "\t")
  ranks <- de$log2FoldChange
  names(ranks) <- toupper(de$geneID)
  ranks <- na.omit(ranks)
  ranks <- sort(ranks, decreasing = TRUE)

  # Get comparison name
  file_name <- basename(de_file)
  file_name <- gsub("_DE_Results\\.tsv$", "", file_name)
  comp_dir  <- file.path(fgsea_dir, paste0("fGSEA_", file_name))
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  cat_names <- names(gene_set_list)

  for (i in seq_along(gene_set_list)) {
    cat("  ", file_name, " - ", cat_names[i], "\n")
    fgsea_res <- fgseaMultilevel(
      pathways      = gene_set_list[[i]],
      stats         = ranks,
      eps           = 0,
      minSize       = fgsea_min,
      maxSize       = fgsea_max,
      nPermSimple   = fgsea_nperm,
      nproc         = fgsea_nproc
    )
    fgsea_res %<>% mutate(direction = case_when(
      padj < fgsea_padj & NES > 0 ~ "upReg",
      padj < fgsea_padj & NES < 0 ~ "downReg",
      TRUE ~ "none"
    ))

    out_file <- file.path(comp_dir,
                          paste0("fgsea_Results_", file_name, "_", cat_names[i], ".tsv"))
    fwrite(fgsea_res, out_file, sep = "\t")

    # --- Bar plot of top pathways (using plot_config) ---
    sig_res <- fgsea_res %>% filter(padj < fgsea_padj)
    if (nrow(sig_res) > 0) {
      top_n <- pcfg$fgsea$top_n_plot %||% 20
      top_up <- sig_res %>% filter(NES > 0) %>% arrange(padj) %>% head(top_n / 2)
      top_down <- sig_res %>% filter(NES < 0) %>% arrange(padj) %>% head(top_n / 2)
      top_paths <- bind_rows(top_up, top_down) %>% arrange(NES)

      if (nrow(top_paths) > 0) {
        vcols <- get_volcano_colors(pcfg)
        p_bar <- ggplot(top_paths, aes(x = reorder(pathway, NES), y = NES,
                                        fill = ifelse(NES > 0, "up", "down"))) +
          geom_col(width = pcfg$fgsea$bar_width %||% 0.7) +
          scale_fill_manual(values = c(up = unname(vcols["up"]), down = unname(vcols["down"])),
                            guide = "none") +
          coord_flip() +
          labs(x = "", y = "Normalized Enrichment Score (NES)",
               title = paste0("fGSEA: ", file_name, " — ", cat_names[i])) +
          theme_pipeline(pcfg) +
          theme(axis.text.y = element_text(size = 8))

        save_plot(p_bar,
                  paste0("fgsea_barplot_", file_name, "_", cat_names[i]),
                  type = "fgsea_bar", outdir = comp_dir, cfg = pcfg)
      }
    }

    save_enrichment_plots(fgsea_res, gene_set_list[[i]], ranks, comp_dir, file_name, cat_names[i])
  }
}

# --- Run on all DE result files ---
de_files <- list.files(deseq2_dir, pattern = "_DE_Results\\.tsv$", full.names = TRUE)
# Exclude "_with_paired" files
de_files <- de_files[!grepl("_with_paired", de_files)]

if (length(de_files) == 0) {
  stop("ERROR: No DE result files found in ", deseq2_dir)
}

cat("Found ", length(de_files), " DE result files.\n")

for (de_file in de_files) {
  cat("Processing: ", basename(de_file), "\n")
  run_fgsea(de_file, msigdb_list)
}

cat("============================================\n")
cat(" Step 5: fGSEA Complete\n")
cat("============================================\n")
