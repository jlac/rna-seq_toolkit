#!/usr/bin/env Rscript

################################################################################
# Expression Heatmap Generator - Command Line Version
# 
# Creates customizable heatmaps from normalized expression data
# All options controllable via command-line arguments
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(viridis)
  library(ggplot2)
  library(reshape2)
})

################################################################################
# COMMAND LINE ARGUMENTS
################################################################################

option_list <- list(
  # Input/Output
  make_option(c("-e", "--expression"), type="character", default="Reformat_TPMCountFile_rsemgenes.txt",
              help="Expression matrix file [default: %default]"),
  make_option(c("-d", "--design"), type="character", default="evitta_design.txt",
              help="Design/metadata file [default: %default]"),
  # Output Options
  make_option(c("-o", "--output"), type="character", default="expression_heatmap.pdf",
              help="Output file (PDF or PNG) [default: %default]"),
  make_option(c("--save_matrix"), action="store_true", default=FALSE,
              help="Save the expression matrix used in the figure to a file"),
  make_option(c("--matrix_file"), type="character", default=NULL,
              help="Output file for expression matrix (default: <output>_matrix.txt)"),
  make_option(c("--duplicate_genes"), type="character", default="make_unique",
              help="How to handle duplicate gene IDs: make_unique, sum, mean, first [default: %default]"),
  
  # Plot Type
  make_option(c("--plot_type"), type="character", default="heatmap",
              help="Type of plot: heatmap, violin, boxplot, volcano [default: %default]"),
  
  # Volcano Plot Input (only used when --plot_type volcano)
  make_option(c("--de_file"), type="character", default=NULL,
              help="Differential expression results table (required for volcano plots)"),
  make_option(c("--de_gene_col"), type="character", default=NULL,
              help="Column with gene IDs in DE table (default: auto-detect, else first column)"),
  make_option(c("--de_lfc_col"), type="character", default=NULL,
              help="Column with log2 fold changes (default: auto-detect log2FoldChange/logFC)"),
  make_option(c("--de_pval_col"), type="character", default=NULL,
              help="Column with p-values to plot (default: auto-detect padj/FDR/adj.P.Val)"),
  
  # Volcano Plot Thresholds
  make_option(c("--lfc_threshold"), type="double", default=1,
              help="Absolute log2FC cutoff for significance [default: %default]"),
  make_option(c("--pval_threshold"), type="double", default=0.05,
              help="P-value cutoff for significance [default: %default]"),
  
  # Volcano Plot Highlighting
  make_option(c("--highlight_genes"), type="character", default=NULL,
              help="Comma-separated genes to highlight and label on the volcano plot"),
  make_option(c("--highlight_file"), type="character", default=NULL,
              help="File with genes to highlight, one per line"),
  make_option(c("--highlight_color"), type="character", default="#D62728",
              help="Color for highlighted genes [default: %default]"),
  make_option(c("--highlight_only"), action="store_true", default=FALSE,
              help="Color only the highlighted genes; draw all others in grey"),
  make_option(c("--label_top_n"), type="integer", default=0,
              help="Also label the N most significant genes [default: %default]"),
  make_option(c("--no_labels"), action="store_true", default=FALSE,
              help="Highlight genes with color but do not draw text labels"),
  
  # Volcano Plot Appearance
  make_option(c("--volcano_colors"), type="character", default="#3C6FBF,#BFBFBF,#C0392B",
              help="Down,NS,Up colors as comma-separated hex [default: %default]"),
  make_option(c("--point_size"), type="double", default=1.2,
              help="Point size for volcano plot [default: %default]"),
  make_option(c("--label_size"), type="double", default=3.5,
              help="Text size for volcano gene labels [default: %default]"),
  make_option(c("--plot_title"), type="character", default=NULL,
              help="Title for the plot (default: auto-generated)"),
  make_option(c("--xlim"), type="character", default=NULL,
              help="X-axis limits for volcano as 'min,max' (default: auto)"),
  make_option(c("--ylim"), type="character", default=NULL,
              help="Y-axis limits for volcano as 'min,max' (default: auto)"),
  
  # Sample Filtering
  make_option(c("-f", "--filter_column"), type="character", default="Pop",
              help="Column name to filter samples on (use 'none' for no filtering) [default: %default]"),
  make_option(c("-v", "--filter_values"), type="character", default="PopNeg,PopPos",
              help="Comma-separated values to keep (e.g., 'PopNeg,PopPos') [default: %default]"),
  
  # Violin/Box Plot Options (only used when --plot_type is violin or boxplot)
  make_option(c("--group_by"), type="character", default="Group",
              help="Column to group samples by for violin/boxplot [default: %default]"),
  make_option(c("--violin_groups"), type="character", default=NULL,
              help="Comma-separated list of groups to include in violin/boxplot (default: all groups)"),
  make_option(c("--violin_colors"), type="character", default=NULL,
              help="Comma-separated hex colors for violin/boxplot groups (optional)"),
  make_option(c("--add_points"), action="store_true", default=FALSE,
              help="Add individual data points to violin/boxplot"),
  make_option(c("--add_stats"), action="store_true", default=FALSE,
              help="Add statistical comparisons to violin/boxplot (requires ggsignif package)"),
  
  # Gene Selection
  make_option(c("-n", "--n_genes"), type="integer", default=1000,
              help="Number of top variable genes to plot [default: %default]"),
  make_option(c("-g", "--gene_list"), type="character", default=NULL,
              help="Comma-separated list of specific genes to plot (overrides -n)"),
  make_option(c("--gene_file"), type="character", default=NULL,
              help="File with gene list (one per line, overrides -n and -g)"),
  
  # DE-based Gene Selection (heatmap/violin/boxplot, requires --de_file)
  make_option(c("--top_de_genes"), type="integer", default=NULL,
              help="Select top N up AND top N down genes from --de_file (gives up to 2N genes)"),
  make_option(c("--top_up"), type="integer", default=NULL,
              help="Number of upregulated genes to select (overrides --top_de_genes)"),
  make_option(c("--top_down"), type="integer", default=NULL,
              help="Number of downregulated genes to select (overrides --top_de_genes)"),
  make_option(c("--de_rank_by"), type="character", default="pvalue",
              help="Rank DE genes by: pvalue or lfc [default: %default]"),
  make_option(c("--de_apply_thresholds"), action="store_true", default=FALSE,
              help="Only consider genes passing --lfc_threshold and --pval_threshold before ranking"),
  make_option(c("--de_order_rows"), action="store_true", default=FALSE,
              help="Order heatmap rows up-then-down by rank instead of clustering genes"),
  
  # Expression Transformation
  make_option(c("-t", "--transform"), type="character", default="log2",
              help="Expression transformation: log2, log10, zscore, none [default: %default]"),
  make_option(c("-p", "--pseudocount"), type="numeric", default=1,
              help="Pseudocount for log transformation [default: %default]"),
  make_option(c("-s", "--scale"), type="character", default="zscore",
              help="Row scaling: zscore, none [default: %default]"),
  
  # Annotations
  make_option(c("-a", "--annotations"), type="character", default="BigGroup,Pop",
              help="Comma-separated column names for annotations (use 'none' for no annotations) [default: %default]"),
  
  # Group Averaging (for heatmaps only)
  make_option(c("--average_groups"), action="store_true", default=FALSE,
              help="Average expression within groups instead of showing individual samples"),
  make_option(c("--average_by"), type="character", default="Group",
              help="Column to group samples by for averaging [default: %default]"),
  make_option(c("--average_function"), type="character", default="mean",
              help="Function for averaging: mean, median [default: %default]"),
  
  # Clustering
  make_option(c("--cluster_rows"), type="logical", default=TRUE,
              help="Cluster genes [default: %default]"),
  make_option(c("--cluster_columns"), type="logical", default=TRUE,
              help="Cluster samples [default: %default]"),
  make_option(c("--no_cluster_rows"), action="store_false", dest="cluster_rows",
              help="Don't cluster genes"),
  make_option(c("--no_cluster_columns"), action="store_false", dest="cluster_columns",
              help="Don't cluster samples"),
  
  # Display Options
  make_option(c("--show_row_names"), action="store_true", default=FALSE,
              help="Show gene names on heatmap"),
  make_option(c("--row_fontsize"), type="numeric", default=6,
              help="Gene name font size [default: %default]"),
  make_option(c("--col_fontsize"), type="numeric", default=8,
              help="Sample name font size [default: %default]"),
  make_option(c("--annotation_fontsize"), type="numeric", default=10,
              help="Annotation label font size [default: %default]"),
  
  # Dimensions
  make_option(c("--width"), type="numeric", default=10,
              help="Heatmap width in inches [default: %default]"),
  make_option(c("--height"), type="numeric", default=12,
              help="Heatmap height in inches [default: %default]"),
  make_option(c("--dpi"), type="numeric", default=300,
              help="DPI for PNG output [default: %default]"),
  
  # Colors
  make_option(c("--color_scheme"), type="character", default="blue_white_red",
              help="Color scheme: blue_white_red, viridis, green_black_red, purple_white_orange [default: %default]"),
  make_option(c("--color_min"), type="numeric", default=-2,
              help="Minimum value for color scale [default: %default]"),
  make_option(c("--color_max"), type="numeric", default=2,
              help="Maximum value for color scale [default: %default]"),
  
  # Other
  make_option(c("--split_by"), type="character", default=NULL,
              help="Column name to split heatmap by (creates separate panels)"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Print detailed progress messages")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="\nGenerate expression heatmaps with flexible filtering and annotation options",
  epilogue=paste(
    "\nExamples:",
    "  # Basic heatmap with defaults",
    "  Rscript create_expression_heatmap_cli.R",
    "",
    "  # Heatmap with groups averaged (cleaner visualization)",
    "  Rscript create_expression_heatmap_cli.R --average_groups --average_by BigGroup -n 1000",
    "",
    "  # Save expression matrix alongside the figure",
    "  Rscript create_expression_heatmap_cli.R -n 1000 --save_matrix -o my_heatmap.pdf",
    "  # Creates: my_heatmap.pdf and my_heatmap_matrix.txt",
    "",
    "  # Heatmap of the top 50 up and top 50 down genes from a DE table",
    "  Rscript create_expression_heatmap_cli.R --de_file deseq2_results.txt --top_de_genes 50 --de_order_rows -o top_de_heatmap.pdf",
    "",
    "  # Rank by raw p-value instead of padj, 25 each direction",
    "  Rscript create_expression_heatmap_cli.R --de_file results.txt --top_de_genes 25 --de_pval_col pvalue -o top25.pdf",
    "",
    "  # Only upregulated genes, restricted to significant ones",
    "  Rscript create_expression_heatmap_cli.R --de_file results.txt --top_up 40 --top_down 0 --de_apply_thresholds -o top_up.pdf",
    "",
    "  # Volcano plot highlighting specific genes",
    "  Rscript create_expression_heatmap_cli.R --plot_type volcano --de_file deseq2_results.txt --highlight_genes \"IFNG,IL6,TNF\" -o volcano.pdf",
    "",
    "  # Volcano with a highlight gene file, everything else greyed out",
    "  Rscript create_expression_heatmap_cli.R --plot_type volcano --de_file results.csv --highlight_file panel.txt --highlight_only -o volcano.pdf",
    "",
    "  # Volcano with custom thresholds and top-10 labeling",
    "  Rscript create_expression_heatmap_cli.R --plot_type volcano --de_file results.txt --lfc_threshold 0.58 --pval_threshold 0.01 --label_top_n 10 -o volcano.pdf",
    "",
    "  # Save matrix to custom filename",
    "  Rscript create_expression_heatmap_cli.R --save_matrix --matrix_file data/expression.txt",
    "",
    "  # Media samples only, top 500 genes",
    "  Rscript create_expression_heatmap_cli.R -f BigGroup -v \"Media_pop-neg,Media_Pop-POS\" -n 500",
    "",
    "  # Average replicates by Pop, then plot",
    "  Rscript create_expression_heatmap_cli.R --average_groups --average_by Pop -f Pop -v \"PopNeg,PopPos\"",
    "",
    "  # Violin plot of specific genes across groups",
    "  Rscript create_expression_heatmap_cli.R --plot_type violin -g \"IFNG,IL6,TNF\" --group_by Group",
    "",
    "  # Boxplot comparing two specific groups",
    "  Rscript create_expression_heatmap_cli.R --plot_type boxplot --gene_file my_genes.txt --group_by Pop --violin_groups \"PopNeg,PopPos\" --add_points",
    "",
    "  # Violin plot with custom colors and statistics",
    "  Rscript create_expression_heatmap_cli.R --plot_type violin -g \"IFNG,IL6,TNF,IL1B\" --group_by BigGroup --violin_colors \"#E41A1C,#377EB8,#4DAF4A\" --add_points --add_stats",
    "",
    "  # All samples, no filtering",
    "  Rscript create_expression_heatmap_cli.R -f none -a \"Group,BigGroup,Pop\"",
    "",
    "  # High-res PNG output",
    "  Rscript create_expression_heatmap_cli.R -o heatmap.png --dpi 600 --width 12 --height 16",
    "",
    "  # Handle duplicate genes by summing",
    "  Rscript create_expression_heatmap_cli.R --duplicate_genes sum",
    sep="\n"
  )
)

opt <- parse_args(opt_parser)

# Verbose logging function
vcat <- function(...) {
  if (opt$verbose) {
    cat(...)
  }
}

################################################################################
# PROCESS ARGUMENTS
################################################################################

cat("Expression Heatmap Generator\n")
cat("============================\n\n")

################################################################################
# DE TABLE LOADER
#
# Shared by volcano plots and by DE-based gene selection (--top_de_genes), so
# both use identical column detection and NA handling.
#
# Returns a list with:
#   $data     data.frame of gene, log2FC, pvalue (NA rows removed)
#   $adjusted TRUE if the p-value column is an adjusted/FDR column
#   $pval_col name of the p-value column that was used
################################################################################

load_de_table <- function(de_path) {
  
  if (!file.exists(de_path)) {
    stop(sprintf("DE results file not found: %s", de_path))
  }
  
  cat(sprintf("Loading DE results from %s...\n", de_path))
  
  # Sniff the delimiter so .csv and .tsv/.txt both work
  first_line <- readLines(de_path, n = 1)
  de_sep <- if (grepl(",", first_line) && !grepl("\t", first_line)) "," else "\t"
  
  de <- read.table(de_path,
                   header = TRUE,
                   sep = de_sep,
                   check.names = FALSE,
                   stringsAsFactors = FALSE,
                   quote = "\"",
                   comment.char = "")
  
  vcat(sprintf("  Columns found: %s\n", paste(colnames(de), collapse = ", ")))
  
  # ---- Identify the gene, log2FC, and p-value columns ----
  
  pick_column <- function(user_choice, candidates, what, fallback = NULL) {
    if (!is.null(user_choice)) {
      if (!user_choice %in% colnames(de)) {
        stop(sprintf("%s column '%s' not found in DE table. Available: %s",
                     what, user_choice, paste(colnames(de), collapse = ", ")))
      }
      return(user_choice)
    }
    hit <- candidates[tolower(candidates) %in% tolower(colnames(de))]
    if (length(hit) > 0) {
      # Return the column as it is actually spelled in the file
      return(colnames(de)[tolower(colnames(de)) == tolower(hit[1])][1])
    }
    if (!is.null(fallback)) return(fallback)
    stop(sprintf("Could not auto-detect the %s column. Specify it explicitly.\nAvailable columns: %s",
                 what, paste(colnames(de), collapse = ", ")))
  }
  
  # DESeq2, edgeR, and limma naming conventions, in preference order
  gene_col_name <- pick_column(
    opt$de_gene_col,
    c("gene", "gene_id", "gene_name", "geneid", "genes", "symbol", "gene_symbol", "id"),
    "Gene ID",
    fallback = colnames(de)[1]
  )
  lfc_col_name <- pick_column(
    opt$de_lfc_col,
    c("log2FoldChange", "logFC", "log2FC", "log2_fold_change", "lfc", "fold_change_log2"),
    "log2 fold change"
  )
  pval_col_name <- pick_column(
    opt$de_pval_col,
    c("padj", "FDR", "adj.P.Val", "qvalue", "q_value", "adj_pvalue", "p_adj",
      "pvalue", "PValue", "P.Value", "p_value", "pval"),
    "p-value"
  )
  
  cat(sprintf("  Gene column:     %s\n", gene_col_name))
  cat(sprintf("  log2FC column:   %s\n", lfc_col_name))
  cat(sprintf("  P-value column:  %s\n", pval_col_name))
  
  # Note whether this is an adjusted or raw p-value, for axis labels and reporting
  is_adjusted <- tolower(pval_col_name) %in%
    tolower(c("padj", "FDR", "adj.P.Val", "qvalue", "q_value", "adj_pvalue", "p_adj"))
  
  # ---- Build the standardized data frame ----
  
  out <- data.frame(
    gene = as.character(de[[gene_col_name]]),
    log2FC = suppressWarnings(as.numeric(de[[lfc_col_name]])),
    pvalue = suppressWarnings(as.numeric(de[[pval_col_name]])),
    stringsAsFactors = FALSE
  )
  
  n_start <- nrow(out)
  
  # Drop rows with missing values. In DESeq2 output, NA padj means the gene was
  # filtered out by independent filtering or flagged as an outlier, so these
  # rows carry no significance call and cannot be ranked or plotted.
  out <- out[!is.na(out$log2FC) & !is.na(out$pvalue), ]
  n_dropped <- n_start - nrow(out)
  if (n_dropped > 0) {
    cat(sprintf("  Dropped %d genes with NA log2FC or p-value\n", n_dropped))
  }
  
  if (nrow(out) == 0) {
    stop("No genes remain after removing NA values. Check the column selections.")
  }
  
  list(data = out, adjusted = is_adjusted, pval_col = pval_col_name)
}

################################################################################
# VOLCANO PLOT
#
# Volcano plots are built from a differential expression results table, not
# from the expression matrix, so this branch runs on its own and exits when
# finished. The expression matrix and design file are not required.
################################################################################

if (opt$plot_type == "volcano") {
  
  if (is.null(opt$de_file)) {
    stop("Volcano plots require a DE results table. Use --de_file.")
  }
  if (!file.exists(opt$de_file)) {
    stop(sprintf("DE results file not found: %s", opt$de_file))
  }
  
  ##############################################################################
  # Load DE results
  ##############################################################################
  
  de_loaded <- load_de_table(opt$de_file)
  volcano_data <- de_loaded$data
  adjusted_pval <- de_loaded$adjusted
  
  # P-values of exactly zero become Inf after -log10. Floor them at the smallest
  # nonzero p-value in the table so those genes stay on the plot; note it, since
  # their height is then a floor rather than a measured value.
  n_zero <- sum(volcano_data$pvalue == 0)
  if (n_zero > 0) {
    min_nonzero <- min(volcano_data$pvalue[volcano_data$pvalue > 0])
    cat(sprintf("  NOTE: %d genes have p-value = 0; capping at %.3g for plotting\n",
                n_zero, min_nonzero))
    volcano_data$pvalue[volcano_data$pvalue == 0] <- min_nonzero
  }
  
  volcano_data$neglog10p <- -log10(volcano_data$pvalue)
  
  ##############################################################################
  # Assign significance categories
  ##############################################################################
  
  volcano_data$category <- "Not significant"
  volcano_data$category[volcano_data$pvalue < opt$pval_threshold &
                          volcano_data$log2FC >= opt$lfc_threshold] <- "Up"
  volcano_data$category[volcano_data$pvalue < opt$pval_threshold &
                          volcano_data$log2FC <= -opt$lfc_threshold] <- "Down"
  volcano_data$category <- factor(volcano_data$category,
                                  levels = c("Down", "Not significant", "Up"))
  
  cat(sprintf("  %d genes plotted: %d up, %d down, %d not significant\n",
              nrow(volcano_data),
              sum(volcano_data$category == "Up"),
              sum(volcano_data$category == "Down"),
              sum(volcano_data$category == "Not significant")))
  
  ##############################################################################
  # Resolve the highlight gene list
  ##############################################################################
  
  highlight_genes <- NULL
  if (!is.null(opt$highlight_file)) {
    if (!file.exists(opt$highlight_file)) {
      stop(sprintf("Highlight gene file not found: %s", opt$highlight_file))
    }
    highlight_genes <- readLines(opt$highlight_file)
    highlight_genes <- trimws(highlight_genes)
    highlight_genes <- highlight_genes[nchar(highlight_genes) > 0]
  } else if (!is.null(opt$highlight_genes)) {
    highlight_genes <- trimws(strsplit(opt$highlight_genes, ",")[[1]])
    highlight_genes <- highlight_genes[nchar(highlight_genes) > 0]
  }
  
  volcano_data$highlight <- FALSE
  
  if (!is.null(highlight_genes)) {
    volcano_data$highlight <- volcano_data$gene %in% highlight_genes
    
    n_found <- sum(volcano_data$highlight)
    cat(sprintf("  Highlighting %d of %d requested genes\n",
                n_found, length(highlight_genes)))
    
    # Report requested genes that are not in the plotted table, separating genes
    # absent from the file entirely from genes dropped above for NA values
    missing <- setdiff(highlight_genes, volcano_data$gene)
    if (length(missing) > 0) {
      in_table <- intersect(missing, as.character(de[[gene_col_name]]))
      not_in_table <- setdiff(missing, in_table)
      if (length(not_in_table) > 0) {
        cat(sprintf("  WARNING: %d highlight genes not in DE table: %s\n",
                    length(not_in_table),
                    paste(head(not_in_table, 10), collapse = ", ")))
      }
      if (length(in_table) > 0) {
        cat(sprintf("  WARNING: %d highlight genes dropped for NA values: %s\n",
                    length(in_table),
                    paste(head(in_table, 10), collapse = ", ")))
      }
    }
    
    if (n_found == 0) {
      warning("None of the highlight genes were found. Check that gene IDs match the DE table.")
    }
  }
  
  # Optionally add the most significant genes to the label set
  volcano_data$label_me <- volcano_data$highlight
  if (opt$label_top_n > 0) {
    ranked <- order(volcano_data$pvalue)
    top_idx <- head(ranked[volcano_data$category[ranked] != "Not significant"],
                    opt$label_top_n)
    volcano_data$label_me[top_idx] <- TRUE
    vcat(sprintf("  Labeling top %d significant genes\n", length(top_idx)))
  }
  
  ##############################################################################
  # Build the plot
  ##############################################################################
  
  cat("Generating volcano plot...\n")
  
  vcolors <- trimws(strsplit(opt$volcano_colors, ",")[[1]])
  if (length(vcolors) != 3) {
    warning("--volcano_colors needs exactly 3 colors (down,ns,up); using defaults")
    vcolors <- c("#3C6FBF", "#BFBFBF", "#C0392B")
  }
  names(vcolors) <- c("Down", "Not significant", "Up")
  
  background <- volcano_data[!volcano_data$highlight, ]
  foreground <- volcano_data[volcano_data$highlight, ]
  
  p <- ggplot(volcano_data, aes(x = log2FC, y = neglog10p))
  
  # Threshold guides
  if (opt$lfc_threshold > 0) {
    p <- p + geom_vline(xintercept = c(-opt$lfc_threshold, opt$lfc_threshold),
                        linetype = "dashed", color = "grey50", linewidth = 0.4)
  }
  p <- p + geom_hline(yintercept = -log10(opt$pval_threshold),
                      linetype = "dashed", color = "grey50", linewidth = 0.4)
  
  if (opt$highlight_only && nrow(foreground) > 0) {
    # Everything grey except the highlighted genes
    p <- p +
      geom_point(data = background, color = "grey80",
                 size = opt$point_size, alpha = 0.6) +
      geom_point(data = foreground, color = opt$highlight_color,
                 size = opt$point_size * 2, alpha = 0.9)
  } else {
    # Standard three-color volcano, with highlighted genes drawn on top
    p <- p +
      geom_point(data = background, aes(color = category),
                 size = opt$point_size, alpha = 0.6)
    if (nrow(foreground) > 0) {
      p <- p +
        geom_point(data = foreground, color = opt$highlight_color,
                   size = opt$point_size * 2.2, alpha = 1) +
        geom_point(data = foreground, shape = 21, fill = NA, color = "black",
                   size = opt$point_size * 2.2, stroke = 0.4)
    }
    p <- p + scale_color_manual(values = vcolors, name = NULL, drop = FALSE)
  }
  
  # Gene labels
  label_data <- volcano_data[volcano_data$label_me, ]
  if (nrow(label_data) > 0 && !opt$no_labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        aes(label = gene),
        size = opt$label_size,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.size = 0.3,
        box.padding = 0.4,
        point.padding = 0.3
      )
    } else {
      # ggrepel gives non-overlapping labels; fall back to plain text without it
      cat("  NOTE: ggrepel not installed, labels may overlap.\n")
      cat("        Install with: install.packages(\"ggrepel\")\n")
      p <- p + geom_text(data = label_data, aes(label = gene),
                         size = opt$label_size, vjust = -0.8, check_overlap = TRUE)
    }
  }
  
  # Axis limits
  if (!is.null(opt$xlim)) {
    xl <- as.numeric(trimws(strsplit(opt$xlim, ",")[[1]]))
    p <- p + coord_cartesian(xlim = xl)
  }
  if (!is.null(opt$ylim)) {
    yl <- as.numeric(trimws(strsplit(opt$ylim, ",")[[1]]))
    p <- p + coord_cartesian(ylim = yl)
  }
  
  plot_title <- if (!is.null(opt$plot_title)) {
    opt$plot_title
  } else {
    sprintf("Volcano plot (|log2FC| > %.2g, %s < %.2g)",
            opt$lfc_threshold,
            if (adjusted_pval) "adj. p" else "p",
            opt$pval_threshold)
  }
  
  p <- p +
    labs(
      x = expression(log[2] ~ "fold change"),
      y = if (adjusted_pval) {
        expression(-log[10] ~ "adjusted p-value")
      } else {
        expression(-log[10] ~ "p-value")
      },
      title = plot_title
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 10)
    )
  
  ##############################################################################
  # Save
  ##############################################################################
  
  cat(sprintf("Saving volcano plot to %s...\n", opt$output))
  
  if (tolower(tools::file_ext(opt$output)) == "png") {
    ggsave(opt$output, p, width = opt$width, height = opt$height, dpi = opt$dpi)
  } else {
    ggsave(opt$output, p, width = opt$width, height = opt$height)
  }
  
  # Save the plotted table, matching the --save_matrix behavior of the other modes
  if (opt$save_matrix) {
    matrix_file <- if (is.null(opt$matrix_file)) {
      paste0(tools::file_path_sans_ext(opt$output), "_data.txt")
    } else {
      opt$matrix_file
    }
    cat(sprintf("Saving volcano data to %s...\n", matrix_file))
    out <- volcano_data[, c("gene", "log2FC", "pvalue", "neglog10p",
                            "category", "highlight")]
    write.table(out, file = matrix_file, sep = "\t",
                quote = FALSE, row.names = FALSE)
  }
  
  ##############################################################################
  # Summary
  ##############################################################################
  
  cat("\n=== SUMMARY ===\n")
  cat(sprintf("Plot type: volcano\n"))
  cat(sprintf("Output file: %s\n", opt$output))
  cat(sprintf("DE table: %s\n", opt$de_file))
  cat(sprintf("Genes plotted: %d\n", nrow(volcano_data)))
  cat(sprintf("Thresholds: |log2FC| >= %.2g, %s < %.2g\n",
              opt$lfc_threshold,
              if (adjusted_pval) "adj. p" else "p",
              opt$pval_threshold))
  cat(sprintf("  Up:   %d\n", sum(volcano_data$category == "Up")))
  cat(sprintf("  Down: %d\n", sum(volcano_data$category == "Down")))
  cat(sprintf("  NS:   %d\n", sum(volcano_data$category == "Not significant")))
  if (sum(volcano_data$highlight) > 0) {
    cat(sprintf("Highlighted genes: %d\n", sum(volcano_data$highlight)))
  }
  cat("\nDone!\n")
  
  quit(save = "no", status = 0)
}

# Check input files exist (heatmap, violin, and boxplot modes)
if (!file.exists(opt$expression)) {
  stop(sprintf("Expression file not found: %s", opt$expression))
}
if (!file.exists(opt$design)) {
  stop(sprintf("Design file not found: %s", opt$design))
}

# Process filter values
if (tolower(opt$filter_column) == "none") {
  filter_column <- NULL
  filter_values <- NULL
  vcat("Sample filtering: DISABLED\n")
} else {
  filter_column <- opt$filter_column
  filter_values <- strsplit(opt$filter_values, ",")[[1]]
  filter_values <- trimws(filter_values)  # Remove whitespace
  vcat(sprintf("Sample filtering: %s in [%s]\n", 
              filter_column, paste(filter_values, collapse=", ")))
}

# Process gene selection
# Precedence: --gene_file > -g > --top_de_genes/--top_up/--top_down > -n
use_top_variable <- FALSE
use_de_selection <- FALSE
gene_list <- NULL

# Resolve the requested up/down counts
n_up   <- if (!is.null(opt$top_up))   opt$top_up   else opt$top_de_genes
n_down <- if (!is.null(opt$top_down)) opt$top_down else opt$top_de_genes
de_selection_requested <- !is.null(n_up) || !is.null(n_down)
if (de_selection_requested) {
  if (is.null(n_up))   n_up   <- 0
  if (is.null(n_down)) n_down <- 0
  if (n_up < 0 || n_down < 0) {
    stop("--top_de_genes / --top_up / --top_down must not be negative")
  }
  if (n_up == 0 && n_down == 0) {
    stop("--top_up and --top_down are both 0; nothing to select")
  }
}

if (!is.null(opt$gene_file)) {
  # Gene list from file
  gene_list <- readLines(opt$gene_file)
  gene_list <- trimws(gene_list)
  gene_list <- gene_list[nchar(gene_list) > 0]  # Remove empty lines
  vcat(sprintf("Gene selection: %d genes from file %s\n", 
              length(gene_list), opt$gene_file))
  if (de_selection_requested) {
    cat("  NOTE: --gene_file takes precedence; ignoring DE-based gene selection\n")
  }
} else if (!is.null(opt$gene_list)) {
  # Gene list from command line
  gene_list <- strsplit(opt$gene_list, ",")[[1]]
  gene_list <- trimws(gene_list)
  vcat(sprintf("Gene selection: %d genes from command line\n", length(gene_list)))
  if (de_selection_requested) {
    cat("  NOTE: -g takes precedence; ignoring DE-based gene selection\n")
  }
} else if (de_selection_requested) {
  # Top up/down genes from a DE results table
  if (is.null(opt$de_file)) {
    stop("DE-based gene selection requires a DE results table. Use --de_file.")
  }
  if (!tolower(opt$de_rank_by) %in% c("pvalue", "lfc")) {
    stop(sprintf("--de_rank_by must be 'pvalue' or 'lfc', got '%s'", opt$de_rank_by))
  }
  use_de_selection <- TRUE
  vcat(sprintf("Gene selection: top %d up / %d down from %s, ranked by %s\n",
              n_up, n_down, opt$de_file, tolower(opt$de_rank_by)))
} else {
  # Top variable genes
  use_top_variable <- TRUE
  vcat(sprintf("Gene selection: Top %d variable genes\n", opt$n_genes))
}

# Process annotations
if (tolower(opt$annotations) == "none") {
  annotation_columns <- NULL
  vcat("Annotations: DISABLED\n")
} else {
  annotation_columns <- strsplit(opt$annotations, ",")[[1]]
  annotation_columns <- trimws(annotation_columns)
  vcat(sprintf("Annotations: %s\n", paste(annotation_columns, collapse=", ")))
}

# Process color scheme
color_scheme <- switch(
  tolower(opt$color_scheme),
  "blue_white_red" = colorRamp2(
    c(opt$color_min, 0, opt$color_max),
    c("blue", "white", "red")
  ),
  "viridis" = viridis(100),
  "green_black_red" = colorRamp2(
    c(opt$color_min, 0, opt$color_max),
    c("green", "black", "red")
  ),
  "purple_white_orange" = colorRamp2(
    c(opt$color_min, 0, opt$color_max),
    c("#8E44AD", "white", "#E67E22")
  ),
  {
    warning(sprintf("Unknown color scheme '%s', using blue_white_red", opt$color_scheme))
    colorRamp2(c(opt$color_min, 0, opt$color_max), c("blue", "white", "red"))
  }
)

vcat(sprintf("Color scheme: %s [%.1f to %.1f]\n", 
            opt$color_scheme, opt$color_min, opt$color_max))

################################################################################
# LOAD DATA
################################################################################

cat("Loading data...\n")

# Load expression data without setting rownames initially
expr_raw <- read.table(opt$expression, 
                       header = TRUE, 
                       sep = "\t",
                       check.names = FALSE,
                       stringsAsFactors = FALSE)

# Get gene column (first column)
gene_col <- expr_raw[, 1]
expr_values <- expr_raw[, -1, drop = FALSE]

# Check for duplicates
duplicated_genes <- gene_col[duplicated(gene_col)]

if (length(duplicated_genes) > 0) {
  cat(sprintf("  WARNING: Found %d duplicate gene IDs\n", length(unique(duplicated_genes))))
  vcat(sprintf("  Examples: %s\n", paste(head(unique(duplicated_genes), 5), collapse=", ")))
  
  if (opt$duplicate_genes == "make_unique") {
    cat("  Making gene IDs unique by adding suffixes...\n")
    gene_col <- make.unique(gene_col, sep = "_")
  } else if (opt$duplicate_genes == "sum") {
    cat("  Aggregating duplicate genes by summing...\n")
    expr_values <- as.data.frame(expr_values)
    expr_values$gene <- gene_col
    expr_values <- aggregate(. ~ gene, data = expr_values, FUN = sum)
    gene_col <- expr_values$gene
    expr_values$gene <- NULL
  } else if (opt$duplicate_genes == "mean") {
    cat("  Aggregating duplicate genes by averaging...\n")
    expr_values <- as.data.frame(expr_values)
    expr_values$gene <- gene_col
    expr_values <- aggregate(. ~ gene, data = expr_values, FUN = mean)
    gene_col <- expr_values$gene
    expr_values$gene <- NULL
  } else if (opt$duplicate_genes == "first") {
    cat("  Keeping first occurrence of duplicate genes...\n")
    keep_idx <- !duplicated(gene_col)
    gene_col <- gene_col[keep_idx]
    expr_values <- expr_values[keep_idx, , drop = FALSE]
  }
}

# Set rownames
rownames(expr_values) <- gene_col
expr_data <- expr_values

design <- read.table(opt$design, 
                    header = TRUE, 
                    sep = "\t",
                    stringsAsFactors = FALSE)

rownames(design) <- design$Sample

cat(sprintf("  Loaded %d genes x %d samples\n", nrow(expr_data), ncol(expr_data)))

################################################################################
# FILTER SAMPLES
################################################################################

if (!is.null(filter_column) && !is.null(filter_values)) {
  cat("Filtering samples...\n")
  
  if (!filter_column %in% colnames(design)) {
    stop(sprintf("Filter column '%s' not found in design file. Available: %s",
                filter_column, paste(colnames(design), collapse=", ")))
  }
  
  samples_to_keep <- design$Sample[design[[filter_column]] %in% filter_values]
  
  if (length(samples_to_keep) == 0) {
    stop(sprintf("No samples match filter criteria: %s in [%s]",
                filter_column, paste(filter_values, collapse=", ")))
  }
  
  expr_data <- expr_data[, samples_to_keep, drop = FALSE]
  design <- design[design$Sample %in% samples_to_keep, ]
  
  cat(sprintf("  Retained %d samples\n", ncol(expr_data)))
}

# Ensure design and expression have matching samples in same order
common_samples <- intersect(colnames(expr_data), design$Sample)
expr_data <- expr_data[, common_samples, drop = FALSE]
design <- design[match(common_samples, design$Sample), ]

################################################################################
# SELECT GENES
################################################################################

cat("Selecting genes...\n")

de_row_order <- NULL   # set by DE selection, used by --de_order_rows

if (use_top_variable) {
  gene_vars <- apply(expr_data, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(opt$n_genes, length(gene_vars))])
  expr_matrix <- expr_data[top_genes, , drop = FALSE]
  cat(sprintf("  Selected %d most variable genes\n", nrow(expr_matrix)))
  
} else if (use_de_selection) {
  
  de_loaded <- load_de_table(opt$de_file)
  de_tab <- de_loaded$data
  rank_label <- if (de_loaded$adjusted) "adj. p" else "p"
  
  # Keep only DE genes that are actually present in the expression matrix, so
  # the counts reported below reflect what can really be drawn
  n_before <- nrow(de_tab)
  de_tab <- de_tab[de_tab$gene %in% rownames(expr_data), ]
  if (nrow(de_tab) == 0) {
    stop(paste("No genes in the DE table match the expression matrix.",
               "Check that both use the same gene ID type (symbols vs Ensembl IDs)."))
  }
  if (nrow(de_tab) < n_before) {
    cat(sprintf("  %d of %d DE genes found in the expression matrix\n",
               nrow(de_tab), n_before))
  }
  
  # Optionally restrict to genes passing the significance thresholds
  if (opt$de_apply_thresholds) {
    n_pre <- nrow(de_tab)
    de_tab <- de_tab[de_tab$pvalue < opt$pval_threshold &
                       abs(de_tab$log2FC) >= opt$lfc_threshold, ]
    cat(sprintf("  %d of %d genes pass |log2FC| >= %.2g and %s < %.2g\n",
               nrow(de_tab), n_pre, opt$lfc_threshold, rank_label, opt$pval_threshold))
    if (nrow(de_tab) == 0) {
      stop("No genes pass the thresholds. Loosen --lfc_threshold / --pval_threshold, or drop --de_apply_thresholds.")
    }
  }
  
  up_pool   <- de_tab[de_tab$log2FC > 0, ]
  down_pool <- de_tab[de_tab$log2FC < 0, ]
  
  # Rank within each direction. Ties in p-value are common in DE output
  # (especially with adjusted p-values, where many genes share a value), so
  # break them by fold-change magnitude to make the selection deterministic.
  rank_pool <- function(pool) {
    if (nrow(pool) == 0) return(pool)
    if (tolower(opt$de_rank_by) == "lfc") {
      pool[order(-abs(pool$log2FC), pool$pvalue), ]
    } else {
      pool[order(pool$pvalue, -abs(pool$log2FC)), ]
    }
  }
  
  up_sel   <- head(rank_pool(up_pool),   n_up)
  down_sel <- head(rank_pool(down_pool), n_down)
  
  # Warn when a pool could not supply as many genes as requested
  if (n_up > 0 && nrow(up_sel) < n_up) {
    cat(sprintf("  WARNING: only %d upregulated genes available (%d requested)\n",
               nrow(up_sel), n_up))
  }
  if (n_down > 0 && nrow(down_sel) < n_down) {
    cat(sprintf("  WARNING: only %d downregulated genes available (%d requested)\n",
               nrow(down_sel), n_down))
  }
  
  selected <- rbind(up_sel, down_sel)
  if (nrow(selected) == 0) {
    stop("DE-based selection produced no genes.")
  }
  
  expr_matrix <- expr_data[selected$gene, , drop = FALSE]
  de_row_order <- selected$gene   # up first, then down, each in rank order
  
  cat(sprintf("  Selected %d genes: %d up, %d down (ranked by %s)\n",
             nrow(expr_matrix), nrow(up_sel), nrow(down_sel),
             if (tolower(opt$de_rank_by) == "lfc") "|log2FC|" else rank_label))
  
  # If thresholds were not applied, say how many of the chosen genes would
  # actually pass them, since ranking alone does not guarantee significance
  if (!opt$de_apply_thresholds) {
    n_sig <- sum(selected$pvalue < opt$pval_threshold &
                   abs(selected$log2FC) >= opt$lfc_threshold)
    if (n_sig < nrow(selected)) {
      cat(sprintf("  NOTE: %d of %d selected genes pass |log2FC| >= %.2g and %s < %.2g\n",
                 n_sig, nrow(selected), opt$lfc_threshold,
                 rank_label, opt$pval_threshold))
    }
  }
  
  vcat(sprintf("  Top up:   %s\n", paste(head(up_sel$gene, 5), collapse = ", ")))
  vcat(sprintf("  Top down: %s\n", paste(head(down_sel$gene, 5), collapse = ", ")))
  
} else {
  available_genes <- intersect(gene_list, rownames(expr_data))
  
  if (length(available_genes) == 0) {
    stop("None of the specified genes found in expression data!")
  }
  
  if (length(available_genes) < length(gene_list)) {
    missing <- setdiff(gene_list, available_genes)
    warning(sprintf("Missing %d genes: %s", 
                   length(missing), 
                   paste(head(missing, 10), collapse=", ")))
  }
  
  expr_matrix <- expr_data[available_genes, , drop = FALSE]
  cat(sprintf("  Using %d genes (%d requested)\n", 
             nrow(expr_matrix), length(gene_list)))
}

################################################################################
# TRANSFORM EXPRESSION VALUES
################################################################################

if (opt$transform != "none") {
  cat(sprintf("Applying %s transformation...\n", opt$transform))
  
  if (opt$transform == "log2") {
    expr_matrix <- log2(expr_matrix + opt$pseudocount)
  } else if (opt$transform == "log10") {
    expr_matrix <- log10(expr_matrix + opt$pseudocount)
  } else if (opt$transform == "zscore") {
    expr_matrix <- t(scale(t(expr_matrix)))
  }
}

################################################################################
# SCALE ROWS
################################################################################

if (opt$scale == "zscore") {
  cat("Applying z-score scaling per gene...\n")
  
  # Check for genes with zero variance (will cause Inf in z-score)
  gene_sds <- apply(expr_matrix, 1, sd, na.rm = TRUE)
  zero_var_genes <- gene_sds == 0 | is.na(gene_sds)
  
  if (any(zero_var_genes)) {
    n_zero <- sum(zero_var_genes)
    cat(sprintf("  WARNING: Removing %d genes with zero variance\n", n_zero))
    vcat(sprintf("  Examples: %s\n", paste(head(rownames(expr_matrix)[zero_var_genes], 5), collapse=", ")))
    expr_matrix <- expr_matrix[!zero_var_genes, , drop = FALSE]
    
    if (nrow(expr_matrix) == 0) {
      stop("All genes have zero variance! Try using -t log2 or -s none")
    }
  }
  
  expr_matrix <- t(scale(t(expr_matrix)))
}

# Handle any remaining infinite or NA values
n_inf <- sum(is.infinite(expr_matrix))
n_na <- sum(is.na(expr_matrix))

if (n_inf > 0) {
  cat(sprintf("  WARNING: Replacing %d infinite values with NA\n", n_inf))
  expr_matrix[is.infinite(expr_matrix)] <- NA
}

if (n_na > 0) {
  cat(sprintf("  WARNING: Found %d NA values (%.2f%% of matrix)\n", 
             n_na, 100 * n_na / length(expr_matrix)))
  
  # If too many NAs, warn the user
  if (n_na / length(expr_matrix) > 0.1) {
    warning("More than 10% of values are NA. Consider using -t log2 transformation.")
  }
  
  # Replace NAs with 0 for visualization
  expr_matrix[is.na(expr_matrix)] <- 0
}

################################################################################
# AVERAGE GROUPS (if requested, heatmap only)
################################################################################

if (opt$average_groups && opt$plot_type == "heatmap") {
  cat("Averaging expression within groups...\n")
  
  # Check that average_by column exists
  if (!opt$average_by %in% colnames(design)) {
    stop(sprintf("Averaging column '%s' not found in design. Available: %s",
                opt$average_by, paste(colnames(design), collapse=", ")))
  }
  
  # Get grouping information for current samples
  sample_groups <- design[colnames(expr_matrix), opt$average_by]
  
  # Determine averaging function
  avg_func <- switch(
    tolower(opt$average_function),
    "mean" = function(x) mean(x, na.rm = TRUE),
    "median" = function(x) median(x, na.rm = TRUE),
    {
      warning(sprintf("Unknown averaging function '%s', using mean", opt$average_function))
      function(x) mean(x, na.rm = TRUE)
    }
  )
  
  # Get unique groups
  unique_groups <- unique(sample_groups)
  
  cat(sprintf("  Collapsing %d samples into %d groups\n", 
             ncol(expr_matrix), length(unique_groups)))
  vcat(sprintf("  Groups: %s\n", paste(unique_groups, collapse=", ")))
  
  # Create new matrix with averaged values
  averaged_matrix <- matrix(NA, 
                            nrow = nrow(expr_matrix), 
                            ncol = length(unique_groups))
  rownames(averaged_matrix) <- rownames(expr_matrix)
  colnames(averaged_matrix) <- unique_groups
  
  # Average expression for each group
  for (group in unique_groups) {
    group_samples <- colnames(expr_matrix)[sample_groups == group]
    
    if (length(group_samples) == 1) {
      # Only one sample, just copy
      averaged_matrix[, group] <- expr_matrix[, group_samples]
    } else {
      # Multiple samples, average
      averaged_matrix[, group] <- apply(expr_matrix[, group_samples, drop = FALSE], 
                                       1, avg_func)
    }
    
    vcat(sprintf("  %s: %d samples averaged\n", group, length(group_samples)))
  }
  
  # Replace expr_matrix with averaged version
  expr_matrix <- averaged_matrix
  
  # Update design to have one row per group
  design_averaged <- design[match(unique_groups, design[[opt$average_by]]), ]
  rownames(design_averaged) <- unique_groups
  design <- design_averaged
  
  cat(sprintf("  Final matrix: %d genes x %d groups\n", 
             nrow(expr_matrix), ncol(expr_matrix)))
}

################################################################################
# SAVE EXPRESSION MATRIX (if requested)
################################################################################

if (opt$save_matrix) {
  # Determine output filename
  if (is.null(opt$matrix_file)) {
    # Generate default filename from output file
    base_name <- tools::file_path_sans_ext(opt$output)
    matrix_file <- paste0(base_name, "_matrix.txt")
  } else {
    matrix_file <- opt$matrix_file
  }
  
  cat(sprintf("Saving expression matrix to %s...\n", matrix_file))
  
  # Create output dataframe with gene names as first column
  output_matrix <- data.frame(
    gene_id = rownames(expr_matrix),
    expr_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Write to file
  write.table(output_matrix, 
              file = matrix_file, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
  
  cat(sprintf("  Saved %d genes x %d %s\n", 
             nrow(expr_matrix), 
             ncol(expr_matrix),
             if (opt$average_groups && opt$plot_type == "heatmap") "groups" else "samples"))
}

################################################################################
# CREATE ANNOTATIONS
################################################################################

ha <- NULL
if (!is.null(annotation_columns) && opt$plot_type == "heatmap") {
  cat("Creating column annotations...\n")
  
  # Check that annotation columns exist
  missing_cols <- setdiff(annotation_columns, colnames(design))
  if (length(missing_cols) > 0) {
    stop(sprintf("Annotation columns not found in design: %s\nAvailable: %s",
                paste(missing_cols, collapse=", "),
                paste(colnames(design), collapse=", ")))
  }
  
  annotation_df <- design[, annotation_columns, drop = FALSE]
  rownames(annotation_df) <- design$Sample
  
  # Auto-generate colors for annotations
  ha_colors <- list()
  for (col in annotation_columns) {
    unique_vals <- unique(annotation_df[[col]])
    n_levels <- length(unique_vals)
    
    # Use different color palettes for different numbers of levels
    if (n_levels <= 8) {
      colors <- brewer.pal(max(3, n_levels), "Set2")[1:n_levels]
    } else if (n_levels <= 12) {
      colors <- brewer.pal(n_levels, "Set3")
    } else {
      colors <- rainbow(n_levels)
    }
    
    ha_colors[[col]] <- setNames(colors, unique_vals)
  }
  
  ha <- HeatmapAnnotation(
    df = annotation_df,
    col = ha_colors,
    annotation_name_gp = gpar(fontsize = opt$annotation_fontsize),
    simple_anno_size = unit(0.5, "cm")
  )
}

################################################################################
# GENERATE PLOT BASED ON TYPE
################################################################################

if (opt$plot_type %in% c("violin", "boxplot")) {
  
  ############################################################################
  # VIOLIN / BOXPLOT
  ############################################################################
  
  cat(sprintf("Generating %s plot...\n", opt$plot_type))
  
  # Violin/boxplot requires gene list
  if (use_top_variable) {
    stop("Violin/boxplot requires a specific gene list. Use -g or --gene_file option.")
  }
  
  # Check group_by column exists
  if (!opt$group_by %in% colnames(design)) {
    stop(sprintf("Group column '%s' not found in design. Available: %s",
                opt$group_by, paste(colnames(design), collapse=", ")))
  }
  
  # Filter groups if specified
  if (!is.null(opt$violin_groups)) {
    violin_groups <- strsplit(opt$violin_groups, ",")[[1]]
    violin_groups <- trimws(violin_groups)
    
    # Filter design and expression
    keep_samples <- design$Sample[design[[opt$group_by]] %in% violin_groups]
    design <- design[design$Sample %in% keep_samples, ]
    expr_matrix <- expr_matrix[, keep_samples, drop = FALSE]
    
    cat(sprintf("  Filtered to %d groups: %s\n", 
               length(violin_groups), paste(violin_groups, collapse=", ")))
  }
  
  # Prepare data for plotting
  plot_data <- as.data.frame(t(expr_matrix))
  plot_data$Sample <- rownames(plot_data)
  plot_data$Group <- design[plot_data$Sample, opt$group_by]
  
  # Melt to long format
  plot_data_long <- melt(plot_data, 
                         id.vars = c("Sample", "Group"),
                         variable.name = "Gene",
                         value.name = "Expression")
  
  # Ensure Group is a factor with proper ordering
  plot_data_long$Group <- factor(plot_data_long$Group)
  
  # Parse colors if provided
  if (!is.null(opt$violin_colors)) {
    custom_colors <- strsplit(opt$violin_colors, ",")[[1]]
    custom_colors <- trimws(custom_colors)
  } else {
    # Auto-generate colors
    n_groups <- length(unique(plot_data_long$Group))
    if (n_groups <= 8) {
      custom_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
    } else {
      custom_colors <- rainbow(n_groups)
    }
  }
  
  # Create base plot
  p <- ggplot(plot_data_long, aes(x = Gene, y = Expression, fill = Group))
  
  if (opt$plot_type == "violin") {
    p <- p + geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(0.9))
  } else {
    p <- p + geom_boxplot(alpha = 0.7, position = position_dodge(0.9), outlier.shape = NA)
  }
  
  # Add individual points if requested
  if (opt$add_points) {
    p <- p + geom_point(aes(group = Group), 
                       position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
                       alpha = 0.5, size = 1)
  }
  
  # Add statistical comparisons if requested
  if (opt$add_stats) {
    if (requireNamespace("ggsignif", quietly = TRUE)) {
      # Simple pairwise comparisons for first two groups
      p <- p + ggsignif::geom_signif(
        comparisons = list(levels(plot_data_long$Group)[1:2]),
        map_signif_level = TRUE,
        test = "t.test"
      )
    } else {
      warning("ggsignif package not installed, skipping statistical annotations")
    }
  }
  
  # Styling
  p <- p + 
    scale_fill_manual(values = custom_colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    ) +
    labs(
      x = "Gene",
      y = if (opt$scale == "zscore") "Z-score" else if (opt$transform == "log2") "log2(TPM + 1)" else "Expression",
      fill = opt$group_by,
      title = sprintf("%s Plot - %s by %s", 
                     tools::toTitleCase(opt$plot_type), 
                     nrow(expr_matrix), 
                     opt$group_by)
    )
  
  # Save plot
  cat(sprintf("Saving %s plot to %s...\n", opt$plot_type, opt$output))
  
  output_ext <- tolower(tools::file_ext(opt$output))
  
  # Calculate dimensions based on number of genes
  n_genes <- length(unique(plot_data_long$Gene))
  plot_width <- max(opt$width, 2 + n_genes * 0.8)  # Scale width with genes
  
  if (output_ext == "png") {
    ggsave(opt$output, p, width = plot_width, height = opt$height, dpi = opt$dpi)
  } else {
    ggsave(opt$output, p, width = plot_width, height = opt$height)
  }
  
} else {
  
  ############################################################################
  # HEATMAP
  ############################################################################

cat("Generating heatmap...\n")

# Determine column split if specified
column_split <- NULL
if (!is.null(opt$split_by)) {
  if (!opt$split_by %in% colnames(design)) {
    warning(sprintf("Split column '%s' not found in design, ignoring", opt$split_by))
  } else {
    column_split <- design[[opt$split_by]]
    cat(sprintf("  Splitting columns by: %s\n", opt$split_by))
  }
}

# Check if clustering can be performed
can_cluster_rows <- opt$cluster_rows
can_cluster_cols <- opt$cluster_columns

# Order rows by DE rank (up first, then down) instead of clustering them.
# Restrict to rows still present: zero-variance genes may have been dropped
# during scaling after selection.
row_split <- NULL
if (opt$de_order_rows) {
  if (is.null(de_row_order)) {
    warning("--de_order_rows has no effect without DE-based gene selection, ignoring")
  } else {
    kept <- de_row_order[de_row_order %in% rownames(expr_matrix)]
    expr_matrix <- expr_matrix[kept, , drop = FALSE]
    # Label each row's direction so the two blocks are visually separated
    row_split <- factor(
      ifelse(kept %in% up_sel$gene, "Up", "Down"),
      levels = c("Up", "Down")
    )
    can_cluster_rows <- FALSE
    cat("  Rows ordered by DE rank (up, then down); row clustering disabled\n")
  }
}

if (can_cluster_rows && nrow(expr_matrix) < 2) {
  warning("Cannot cluster rows with less than 2 genes, disabling row clustering")
  can_cluster_rows <- FALSE
}

if (can_cluster_cols && ncol(expr_matrix) < 2) {
  warning("Cannot cluster columns with less than 2 samples, disabling column clustering")
  can_cluster_cols <- FALSE
}

# Final check for remaining issues
if (can_cluster_rows || can_cluster_cols) {
  if (any(is.na(expr_matrix)) || any(is.infinite(expr_matrix))) {
    warning("Matrix contains NA or Inf values, clustering may fail. Consider using -t log2 or --no_cluster_rows --no_cluster_columns")
  }
}

ht <- Heatmap(
  expr_matrix,
  name = "Expression",
  
  # Colors
  col = color_scheme,
  
  # Clustering
  cluster_rows = can_cluster_rows,
  cluster_columns = can_cluster_cols,
  
  # Splitting
  column_split = column_split,
  row_split = row_split,
  
  # Annotations
  top_annotation = ha,
  
  # Row names (genes)
  show_row_names = opt$show_row_names,
  row_names_gp = gpar(fontsize = opt$row_fontsize),
  
  # Column names (samples)
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = opt$col_fontsize),
  column_names_rot = 45,
  
  # Legend
  heatmap_legend_param = list(
    title = if (opt$scale == "zscore") "Z-score" else "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

################################################################################
# SAVE OUTPUT
################################################################################

cat(sprintf("Saving heatmap to %s...\n", opt$output))

# Determine output format
output_ext <- tolower(tools::file_ext(opt$output))

if (output_ext == "png") {
  png(opt$output, 
      width = opt$width, 
      height = opt$height, 
      units = "in", 
      res = opt$dpi)
  draw(ht)
  dev.off()
} else {
  # Default to PDF
  pdf(opt$output, width = opt$width, height = opt$height)
  draw(ht)
  dev.off()
}

}  # End of if/else for plot_type

################################################################################
# SUMMARY
################################################################################

cat("\n=== SUMMARY ===\n")
cat(sprintf("Plot type: %s\n", opt$plot_type))
cat(sprintf("Output file: %s\n", opt$output))
if (opt$save_matrix) {
  if (is.null(opt$matrix_file)) {
    matrix_file <- paste0(tools::file_path_sans_ext(opt$output), "_matrix.txt")
  } else {
    matrix_file <- opt$matrix_file
  }
  cat(sprintf("Matrix file: %s\n", matrix_file))
}
cat(sprintf("Genes: %d\n", nrow(expr_matrix)))
if (opt$average_groups && opt$plot_type == "heatmap") {
  cat(sprintf("Groups: %d (averaged from individual samples)\n", ncol(expr_matrix)))
  cat(sprintf("Averaging method: %s by %s\n", opt$average_function, opt$average_by))
} else {
  cat(sprintf("Samples: %d\n", ncol(expr_matrix)))
}
cat(sprintf("Expression range: [%.2f, %.2f]\n", 
           min(expr_matrix, na.rm=TRUE), 
           max(expr_matrix, na.rm=TRUE)))

if (opt$plot_type == "heatmap") {
  if (!is.null(annotation_columns)) {
    cat("\nSample distribution:\n")
    for (col in annotation_columns) {
      cat(sprintf("  %s:\n", col))
      tab <- table(design[[col]])
      for (i in 1:length(tab)) {
        cat(sprintf("    %s: %d\n", names(tab)[i], tab[i]))
      }
    }
  }
} else {
  cat(sprintf("\nGrouping by: %s\n", opt$group_by))
  cat("Group distribution:\n")
  tab <- table(design[[opt$group_by]])
  for (i in 1:length(tab)) {
    cat(sprintf("  %s: %d samples\n", names(tab)[i], tab[i]))
  }
}

cat("\nDone!\n")
