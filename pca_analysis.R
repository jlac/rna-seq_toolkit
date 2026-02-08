#!/usr/bin/env Rscript

# =============================================================================
# PCA Analysis of RNA-seq Raw Counts
# Filters low-count genes, normalizes with DESeq2 VST, runs PCA,
# and outputs top genes driving each principal component.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
})

# =============================================================================
# Command-line arguments
# =============================================================================
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Path to raw counts matrix (genes x samples, tab or comma delimited) [required]"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata file (tab or comma delimited). First column = sample IDs matching counts column names. Optional but recommended for colored PCA."),
  make_option(c("--color_by"), type = "character", default = NULL,
              help = "Column name in metadata to color PCA points by [default: first column after sample ID]"),
  make_option(c("--shape_by"), type = "character", default = NULL,
              help = "Column name in metadata to shape PCA points by [default: none]"),
  make_option(c("--dup_method"), type = "character", default = "sum",
              help = "How to handle duplicate gene IDs: 'sum' (combine counts) or 'make_unique' (append suffix) [default: %default]"),
  make_option(c("--min_count"), type = "integer", default = 10,
              help = "Minimum count threshold for filtering [default: %default]"),
  make_option(c("--min_samples"), type = "integer", default = 3,
              help = "Minimum number of samples that must meet --min_count [default: %default]"),
  make_option(c("--min_samples_frac"), type = "double", default = NULL,
              help = "Fraction of samples that must meet --min_count (overrides --min_samples if set) [default: NULL]"),
  make_option(c("--norm_method"), type = "character", default = "vst",
              help = "Normalization method: 'vst' or 'rlog' [default: %default]"),
  make_option(c("--n_top_genes"), type = "integer", default = 100,
              help = "Number of top genes to report per PC [default: %default]"),
  make_option(c("--n_pcs"), type = "integer", default = 5,
              help = "Number of PCs to extract top genes for [default: %default]"),
  make_option(c("--n_top_var"), type = "integer", default = NULL,
              help = "If set, use only the top N most variable genes for PCA [default: use all filtered genes]"),
  make_option(c("--label_samples"), action = "store_true", default = FALSE,
              help = "Label sample points on PCA plot [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "pca_results",
              help = "Output directory [default: %default]"),
  make_option(c("--prefix"), type = "character", default = "pca",
              help = "Output file prefix [default: %default]"),
  make_option(c("--width"), type = "double", default = 8,
              help = "Plot width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "Plot height in inches [default: %default]"),
  make_option(c("--point_size"), type = "double", default = 3,
              help = "Point size on PCA plot [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "PCA analysis of RNA-seq raw counts with top gene loadings per PC")
opt <- parse_args(opt_parser)

if (is.null(opt$counts)) {
  print_help(opt_parser)
  stop("--counts is required.", call. = FALSE)
}

# =============================================================================
# Helper: detect delimiter and read file
# =============================================================================
read_delim_auto <- function(filepath, resolve_dup_genes = FALSE, dup_method = "sum") {
  first_line <- readLines(filepath, n = 1)
  if (grepl("\t", first_line)) {
    sep <- "\t"
  } else if (grepl(",", first_line)) {
    sep <- ","
  } else {
    sep <- "\t"  # default
  }
  
  # Read without row.names to handle potential duplicates
  df <- read.table(filepath, header = TRUE, sep = sep,
                   check.names = FALSE, stringsAsFactors = FALSE)
  
  if (resolve_dup_genes) {
    gene_col <- df[, 1]
    data_cols <- df[, -1, drop = FALSE]
    
    n_dup <- sum(duplicated(gene_col))
    if (n_dup > 0) {
      cat(sprintf("  Found %d duplicate gene IDs (%d unique of %d total).\n",
                  n_dup, length(unique(gene_col)), length(gene_col)))
      
      if (dup_method == "sum") {
        cat("  Resolving by summing counts across duplicates...\n")
        data_cols$`.gene_id` <- gene_col
        agg <- aggregate(. ~ `.gene_id`, data = data_cols, FUN = sum)
        rownames(agg) <- agg$`.gene_id`
        agg$`.gene_id` <- NULL
        cat(sprintf("  Collapsed %d rows -> %d unique genes\n", nrow(df), nrow(agg)))
        return(as.data.frame(agg))
      } else {
        cat("  Resolving by appending numeric suffix to duplicates...\n")
        gene_col <- make.unique(gene_col, sep = "_dup")
        rownames(data_cols) <- gene_col
        return(data_cols)
      }
    } else {
      rownames(data_cols) <- gene_col
      return(data_cols)
    }
  } else {
    # For metadata: first column as rownames (no duplicates expected)
    rownames(df) <- df[, 1]
    df <- df[, -1, drop = FALSE]
    return(df)
  }
}

# =============================================================================
# 1. Load data
# =============================================================================
cat("=== Loading counts matrix ===\n")
counts_raw <- read_delim_auto(opt$counts, resolve_dup_genes = TRUE, dup_method = opt$dup_method)
cat(sprintf("  Raw matrix: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# Ensure integer counts
counts_raw <- round(as.matrix(counts_raw))
mode(counts_raw) <- "integer"

# Load metadata if provided
meta <- NULL
if (!is.null(opt$metadata)) {
  cat("=== Loading metadata ===\n")
  meta <- read_delim_auto(opt$metadata, resolve_dup_genes = FALSE)
  
  # Match sample order
  common_samples <- intersect(colnames(counts_raw), rownames(meta))
  if (length(common_samples) == 0) {
    stop("No matching sample IDs between counts columns and metadata rows.", call. = FALSE)
  }
  if (length(common_samples) < ncol(counts_raw)) {
    cat(sprintf("  Warning: %d/%d samples matched metadata. Using matched subset.\n",
                length(common_samples), ncol(counts_raw)))
  }
  counts_raw <- counts_raw[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  cat(sprintf("  Metadata: %d samples, %d variables\n", nrow(meta), ncol(meta)))
}

# =============================================================================
# 2. Filter low-count genes
# =============================================================================
cat("=== Filtering genes ===\n")
min_samp <- if (!is.null(opt$min_samples_frac)) {
  ceiling(opt$min_samples_frac * ncol(counts_raw))
} else {
  opt$min_samples
}

keep <- rowSums(counts_raw >= opt$min_count) >= min_samp
counts_filt <- counts_raw[keep, ]
cat(sprintf("  Filter: count >= %d in >= %d samples\n", opt$min_count, min_samp))
cat(sprintf("  Retained: %d / %d genes (%.1f%%)\n",
            nrow(counts_filt), nrow(counts_raw),
            100 * nrow(counts_filt) / nrow(counts_raw)))

if (nrow(counts_filt) < 50) {
  stop("Too few genes passed filtering. Consider relaxing --min_count or --min_samples.", call. = FALSE)
}

# =============================================================================
# 3. Normalize with DESeq2
# =============================================================================
cat(sprintf("=== Normalizing with %s ===\n", toupper(opt$norm_method)))

# Create a minimal DESeqDataSet
if (!is.null(meta) && ncol(meta) > 0) {
  col_data <- meta
} else {
  col_data <- data.frame(sample = colnames(counts_filt), row.names = colnames(counts_filt))
}

dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = col_data,
                              design = ~ 1)

if (opt$norm_method == "rlog") {
  norm_data <- assay(rlog(dds, blind = TRUE))
} else {
  norm_data <- assay(vst(dds, blind = TRUE))
}
cat(sprintf("  Normalized matrix: %d genes x %d samples\n", nrow(norm_data), ncol(norm_data)))

# =============================================================================
# 4. Optionally subset to top variable genes
# =============================================================================
if (!is.null(opt$n_top_var)) {
  cat(sprintf("=== Selecting top %d most variable genes ===\n", opt$n_top_var))
  rv <- rowVars(norm_data)
  names(rv) <- rownames(norm_data)
  top_var_genes <- names(sort(rv, decreasing = TRUE))[1:min(opt$n_top_var, length(rv))]
  pca_input <- norm_data[top_var_genes, ]
} else {
  pca_input <- norm_data
}

# =============================================================================
# 5. Run PCA
# =============================================================================
cat("=== Running PCA ===\n")
pca_res <- prcomp(t(pca_input), center = TRUE, scale. = TRUE)

# Variance explained
var_pct <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
n_pcs <- min(opt$n_pcs, ncol(pca_res$x))
cat("  Variance explained:\n")
for (i in 1:n_pcs) {
  cat(sprintf("    PC%d: %.2f%%\n", i, var_pct[i]))
}

# =============================================================================
# 6. Create output directory
# =============================================================================
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 7. PCA plots
# =============================================================================
cat("=== Generating PCA plots ===\n")

pca_df <- as.data.frame(pca_res$x[, 1:n_pcs])
pca_df$Sample <- rownames(pca_df)

if (!is.null(meta)) {
  pca_df <- cbind(pca_df, meta[rownames(pca_df), , drop = FALSE])
}

# Determine color and shape aesthetics
color_col <- opt$color_by
shape_col <- opt$shape_by

if (!is.null(meta) && is.null(color_col)) {
  color_col <- colnames(meta)[1]
}

# Generate all pairwise PC plots for first few PCs
pc_pairs <- list(c(1, 2))
if (n_pcs >= 3) pc_pairs <- c(pc_pairs, list(c(1, 3), c(2, 3)))

for (pair in pc_pairs) {
  pc_x <- paste0("PC", pair[1])
  pc_y <- paste0("PC", pair[2])
  
  p <- ggplot(pca_df, aes_string(x = pc_x, y = pc_y))
  
  if (!is.null(color_col) && color_col %in% colnames(pca_df)) {
    if (!is.null(shape_col) && shape_col %in% colnames(pca_df)) {
      p <- p + geom_point(aes_string(color = color_col, shape = shape_col),
                          size = opt$point_size)
    } else {
      p <- p + geom_point(aes_string(color = color_col), size = opt$point_size)
    }
  } else {
    p <- p + geom_point(size = opt$point_size, color = "steelblue")
  }
  
  if (opt$label_samples) {
    p <- p + geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 20)
  }
  
  p <- p +
    labs(x = sprintf("%s (%.1f%%)", pc_x, var_pct[pair[1]]),
         y = sprintf("%s (%.1f%%)", pc_y, var_pct[pair[2]]),
         title = sprintf("PCA: %s vs %s", pc_x, pc_y)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor = element_blank())
  
  fname <- file.path(opt$output_dir,
                     sprintf("%s_%s_vs_%s.pdf", opt$prefix, pc_x, pc_y))
  ggsave(fname, p, width = opt$width, height = opt$height)
  cat(sprintf("  Saved: %s\n", fname))
  
  # Also save PNG
  fname_png <- sub("\\.pdf$", ".png", fname)
  ggsave(fname_png, p, width = opt$width, height = opt$height, dpi = 150)
}

# =============================================================================
# 8. Scree plot
# =============================================================================
scree_n <- min(20, length(var_pct))
scree_df <- data.frame(PC = factor(1:scree_n, levels = 1:scree_n),
                       Variance = var_pct[1:scree_n],
                       Cumulative = cumsum(var_pct[1:scree_n]))

p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(aes(x = as.numeric(PC), y = Cumulative / max(Cumulative) * max(Variance)),
            color = "firebrick", linewidth = 0.8) +
  geom_point(aes(x = as.numeric(PC), y = Cumulative / max(Cumulative) * max(Variance)),
             color = "firebrick", size = 2) +
  scale_y_continuous(
    name = "Variance Explained (%)",
    sec.axis = sec_axis(~ . / max(scree_df$Variance) * max(scree_df$Cumulative),
                        name = "Cumulative (%)")
  ) +
  labs(title = "Scree Plot", x = "Principal Component") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

scree_fname <- file.path(opt$output_dir, sprintf("%s_scree.pdf", opt$prefix))
ggsave(scree_fname, p_scree, width = opt$width, height = opt$height)
ggsave(sub("\\.pdf$", ".png", scree_fname), p_scree,
       width = opt$width, height = opt$height, dpi = 150)
cat(sprintf("  Saved: %s\n", scree_fname))

# =============================================================================
# 9. Extract top genes per PC (by rotation/loading magnitude)
# =============================================================================
cat(sprintf("=== Extracting top %d genes per PC ===\n", opt$n_top_genes))

loadings <- pca_res$rotation  # genes x PCs

all_top_genes <- list()

for (i in 1:n_pcs) {
  pc_name <- paste0("PC", i)
  pc_loads <- loadings[, i]
  
  # Sort by absolute loading
  ord <- order(abs(pc_loads), decreasing = TRUE)
  top_idx <- ord[1:min(opt$n_top_genes, length(ord))]
  
  top_df <- data.frame(
    Gene = names(pc_loads)[top_idx],
    Loading = pc_loads[top_idx],
    Abs_Loading = abs(pc_loads[top_idx]),
    Rank = 1:length(top_idx),
    stringsAsFactors = FALSE
  )
  
  all_top_genes[[pc_name]] <- top_df
  
  # Write per-PC file
  out_file <- file.path(opt$output_dir,
                        sprintf("%s_top%d_genes_%s.tsv", opt$prefix, opt$n_top_genes, pc_name))
  write.table(top_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  %s: top gene = %s (loading = %.4f)\n",
              pc_name, top_df$Gene[1], top_df$Loading[1]))
}

# Also write a combined file (wide format - gene lists side by side)
max_rows <- opt$n_top_genes
combined <- data.frame(Rank = 1:max_rows)
for (pc_name in names(all_top_genes)) {
  df <- all_top_genes[[pc_name]]
  combined[[paste0(pc_name, "_Gene")]] <- df$Gene[1:max_rows]
  combined[[paste0(pc_name, "_Loading")]] <- round(df$Loading[1:max_rows], 6)
}
combined_file <- file.path(opt$output_dir,
                           sprintf("%s_top%d_genes_all_PCs.tsv", opt$prefix, opt$n_top_genes))
write.table(combined, combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Combined file: %s\n", combined_file))

# =============================================================================
# 10. Top loading genes heatmap (union of top 20 from each PC)
# =============================================================================
cat("=== Generating top loading genes heatmap ===\n")
top_union <- unique(unlist(lapply(all_top_genes, function(df) df$Gene[1:min(20, nrow(df))])))
cat(sprintf("  Union of top 20 genes per PC: %d unique genes\n", length(top_union)))

heatmap_mat <- norm_data[top_union, ]
# Scale by row
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# Annotation
ann_col <- NULL
if (!is.null(meta) && !is.null(color_col) && color_col %in% colnames(meta)) {
  ann_col <- meta[colnames(heatmap_mat_scaled), color_col, drop = FALSE]
}

heatmap_fname <- file.path(opt$output_dir, sprintf("%s_top_loadings_heatmap.pdf", opt$prefix))
pdf(heatmap_fname, width = max(opt$width, 10), height = max(opt$height, 8))
pheatmap(heatmap_mat_scaled,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         clustering_method = "ward.D2",
         annotation_col = ann_col,
         show_rownames = (length(top_union) <= 60),
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         main = "Top PCA Loading Genes (row-scaled)")
dev.off()
cat(sprintf("  Saved: %s\n", heatmap_fname))

# =============================================================================
# 11. Save PCA coordinates & variance explained
# =============================================================================
pca_coords_file <- file.path(opt$output_dir, sprintf("%s_coordinates.tsv", opt$prefix))
write.table(pca_df, pca_coords_file, sep = "\t", row.names = FALSE, quote = FALSE)

var_file <- file.path(opt$output_dir, sprintf("%s_variance_explained.tsv", opt$prefix))
var_df <- data.frame(PC = paste0("PC", 1:length(var_pct)),
                     Variance_Pct = var_pct,
                     Cumulative_Pct = cumsum(var_pct))
write.table(var_df, var_file, sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 12. Save normalized matrix
# =============================================================================
norm_file <- file.path(opt$output_dir, sprintf("%s_normalized_matrix.tsv", opt$prefix))
write.table(data.frame(Gene = rownames(norm_data), norm_data, check.names = FALSE),
            norm_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved normalized matrix: %s\n", norm_file))

# =============================================================================
# Summary
# =============================================================================
cat("\n=== Done! ===\n")
cat(sprintf("  Input: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  After filtering: %d genes\n", nrow(counts_filt)))
cat(sprintf("  Normalization: %s\n", toupper(opt$norm_method)))
cat(sprintf("  Output directory: %s/\n", opt$output_dir))
cat("  Files generated:\n")
cat("    - PCA scatter plots (PDF + PNG)\n")
cat("    - Scree plot (PDF + PNG)\n")
cat("    - Top loading genes per PC (TSV)\n")
cat("    - Combined top genes all PCs (TSV)\n")
cat("    - Top loadings heatmap (PDF)\n")
cat("    - PCA coordinates (TSV)\n")
cat("    - Variance explained (TSV)\n")
cat("    - Normalized expression matrix (TSV)\n")
