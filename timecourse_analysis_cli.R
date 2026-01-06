#!/usr/bin/env Rscript

################################################################################
# Timecourse Bulk RNA-seq Analysis Script
# Using limma-voom and/or DESeq2
# With full command-line argument support
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

################################################################################
# COMMAND LINE ARGUMENT PARSING
################################################################################

option_list <- list(
  # Required arguments
  make_option(c("-c", "--counts"), type="character", default=NULL, dest="counts",
              help="Path to counts matrix file (required) [genes x samples]", 
              metavar="FILE"),
  
  make_option(c("-d", "--design"), type="character", default=NULL, dest="design",
              help="Path to design/metadata file (required) [samples x variables]", 
              metavar="FILE"),
  
  # Time variable settings
  make_option(c("-t", "--time-column"), type="character", default="Time", dest="time_column",
              help="Column name for time variable in design file [default: %default]",
              metavar="STRING"),
  
  make_option(c("--time-numeric"), action="store_true", default=FALSE, dest="time_numeric",
              help="Treat time as continuous numeric variable [default: %default]"),
  
  make_option(c("--time-categorical"), action="store_false", dest="time_numeric",
              help="Treat time as categorical factor"),
  
  # Modeling approach for continuous time
  make_option(c("--spline-df"), type="integer", default=3, dest="spline_df",
              help="Degrees of freedom for natural splines (use 0 for polynomial) [default: %default]",
              metavar="INT"),
  
  # Method selection
  make_option(c("--method"), type="character", default="both", dest="method",
              help="Analysis method: 'limma', 'deseq2', or 'both' [default: %default]",
              metavar="STRING"),
  
  # Covariates
  make_option(c("--covariates"), type="character", default=NULL, dest="covariates",
              help="Comma-separated list of covariate column names (e.g., 'Treatment,Batch')",
              metavar="STRING"),
  
  # Differential timecourse testing
  make_option(c("--group-column"), type="character", default=NULL, dest="group_column",
              help="Column name for grouping variable (for differential timecourse analysis)",
              metavar="STRING"),
  
  make_option(c("--groups"), type="character", default=NULL, dest="groups",
              help="Two groups to compare, comma-separated (e.g., 'WT,KO' or 'Control,Treatment')",
              metavar="STRING"),
  
  make_option(c("--test-interaction"), action="store_true", default=FALSE, dest="test_interaction",
              help="Test for time × group interaction (differential timecourse) [default: %default]"),
  
  # Filtering parameters
  make_option(c("--min-count"), type="integer", default=10, dest="min_count",
              help="Minimum count threshold for filtering [default: %default]",
              metavar="INT"),
  
  make_option(c("--min-samples"), type="integer", default=3, dest="min_samples",
              help="Minimum number of samples with min-count [default: %default]",
              metavar="INT"),
  
  # Output settings
  make_option(c("-o", "--output-dir"), type="character", default="timecourse_results", dest="output_dir",
              help="Output directory path [default: %default]",
              metavar="DIR"),
  
  make_option(c("--top-genes"), type="integer", default=50, dest="top_genes",
              help="Number of top genes for visualizations [default: %default]",
              metavar="INT"),
  
  make_option(c("--fdr-threshold"), type="double", default=0.05, dest="fdr_threshold",
              help="FDR threshold for significance [default: %default]",
              metavar="FLOAT"),
  
  # Additional options
  make_option(c("--no-plots"), action="store_true", default=FALSE, dest="no_plots",
              help="Skip generating plots (faster) [default: generate plots]"),
  
  make_option(c("--duplicate-genes"), type="character", default="sum", dest="duplicate_genes",
              help="How to handle duplicate gene IDs: 'sum' (sum counts), 'first' (keep first), or 'unique' (make unique) [default: %default]",
              metavar="STRING"),
  
  make_option(c("--round-counts"), action="store_true", default=FALSE, dest="round_counts",
              help="Round counts to integers (required for DESeq2 with RSEM/Salmon data) [default: auto-detect]"),
  
  make_option(c("--cluster-genes"), action="store_true", default=FALSE, dest="cluster_genes",
              help="Cluster significant genes by temporal pattern [default: FALSE]"),
  
  make_option(c("--n-clusters"), type="integer", default=6, dest="n_clusters",
              help="Number of clusters for gene clustering [default: %default]",
              metavar="INT"),
  
  make_option(c("--cluster-method"), type="character", default="kmeans", dest="cluster_method",
              help="Clustering method: 'kmeans', 'hierarchical', or 'fuzzy' [default: %default]",
              metavar="STRING"),
  
  # Pathway enrichment
  make_option(c("--run-enrichment"), action="store_true", default=FALSE, dest="run_enrichment",
              help="Run pathway enrichment analysis on gene clusters [default: FALSE]"),
  
  make_option(c("--enrichment-db"), type="character", default="reactome", dest="enrichment_db",
              help="Enrichment database: 'reactome', 'gobp', or 'gmt' [default: %default]",
              metavar="STRING"),
  
  make_option(c("--organism"), type="character", default="human", dest="organism",
              help="Organism for enrichment: 'human' or 'mouse' [default: %default]",
              metavar="STRING"),
  
  make_option(c("--gmt-file"), type="character", default=NULL, dest="gmt_file",
              help="Path to GMT file for custom enrichment (required if enrichment-db=gmt)",
              metavar="FILE"),
  
  make_option(c("--enrichment-pval"), type="double", default=0.05, dest="enrichment_pval",
              help="P-value cutoff for enrichment [default: %default]",
              metavar="FLOAT"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, dest="verbose",
              help="Print verbose output [default: %default]")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="\nTimecourse RNA-seq differential expression analysis using limma-voom and/or DESeq2.",
  epilogue=paste(
    "\nExamples:",
    "  # Basic usage with natural splines:",
    "  Rscript timecourse_analysis.R -c counts.txt -d design.txt --time-numeric",
    "",
    "  # With covariates:",
    "  Rscript timecourse_analysis.R -c counts.txt -d design.txt --time-numeric \\",
    "    --covariates Treatment,Batch",
    "",
    "  # Use DESeq2 only with categorical time:",
    "  Rscript timecourse_analysis.R -c counts.txt -d design.txt \\",
    "    --time-categorical --method deseq2",
    "",
    "  # Custom filtering and output:",
    "  Rscript timecourse_analysis.R -c counts.txt -d design.txt --time-numeric \\",
    "    --min-count 20 --min-samples 5 -o my_results",
    "",
    sep="\n"
  )
)

opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$counts)) {
  print_help(opt_parser)
  stop("Error: --counts argument is required", call.=FALSE)
}

if (is.null(opt$design)) {
  print_help(opt_parser)
  stop("Error: --design argument is required", call.=FALSE)
}

# Validate method
if (!opt$method %in% c("limma", "deseq2", "both")) {
  stop("Error: --method must be 'limma', 'deseq2', or 'both'", call.=FALSE)
}

# Parse covariates
covariates <- if (!is.null(opt$covariates)) {
  strsplit(opt$covariates, ",")[[1]]
} else {
  c()
}

# Parse groups for comparison
compare_groups <- if (!is.null(opt$groups)) {
  groups_vec <- strsplit(opt$groups, ",")[[1]]
  if (length(groups_vec) != 2) {
    stop("Error: --groups must specify exactly 2 groups (e.g., 'WT,KO')", call.=FALSE)
  }
  groups_vec
} else {
  NULL
}

# Validate group comparison settings
if (opt$test_interaction && is.null(opt$group_column)) {
  stop("Error: --test-interaction requires --group-column to be specified", call.=FALSE)
}

if (opt$test_interaction && is.null(compare_groups)) {
  stop("Error: --test-interaction requires --groups to be specified", call.=FALSE)
}

# Validate enrichment settings
if (opt$run_enrichment && !opt$cluster_genes) {
  stop("Error: --run-enrichment requires --cluster-genes to be enabled", call.=FALSE)
}

if (opt$run_enrichment && opt$enrichment_db == "gmt" && is.null(opt$gmt_file)) {
  stop("Error: --enrichment-db gmt requires --gmt-file to be specified", call.=FALSE)
}

if (opt$run_enrichment && !opt$enrichment_db %in% c("reactome", "gobp", "gmt")) {
  stop("Error: --enrichment-db must be 'reactome', 'gobp', or 'gmt'", call.=FALSE)
}

if (opt$run_enrichment && !opt$organism %in% c("human", "mouse")) {
  stop("Error: --organism must be 'human' or 'mouse'", call.=FALSE)
}

# Handle spline_df (NULL if 0)
spline_df <- if (opt$spline_df == 0) NULL else opt$spline_df

# Set flags for which methods to run
run_limma <- opt$method %in% c("limma", "both")
run_deseq2 <- opt$method %in% c("deseq2", "both")

# Store parameters in named list for easy access
params <- list(
  counts_file = opt$counts,
  design_file = opt$design,
  time_column = opt$time_column,
  time_as_numeric = opt$time_numeric,
  spline_df = spline_df,
  covariates = covariates,
  group_column = opt$group_column,
  compare_groups = compare_groups,
  test_interaction = opt$test_interaction,
  min_count = opt$min_count,
  min_samples = opt$min_samples,
  output_dir = opt$output_dir,
  top_n_genes = opt$top_genes,
  fdr_threshold = opt$fdr_threshold,
  generate_plots = !opt$no_plots,
  duplicate_genes = opt$duplicate_genes,
  round_counts = opt$round_counts,
  cluster_genes = opt$cluster_genes,
  n_clusters = opt$n_clusters,
  cluster_method = opt$cluster_method,
  run_enrichment = opt$run_enrichment,
  enrichment_db = opt$enrichment_db,
  organism = opt$organism,
  gmt_file = opt$gmt_file,
  enrichment_pval = opt$enrichment_pval,
  verbose = opt$verbose
)

################################################################################
# HELPER FUNCTIONS
################################################################################

vcat <- function(...) {
  if (params$verbose) {
    cat(...)
  }
}

# Create output directories
create_output_dirs <- function(base_dir) {
  dirs <- c(
    base_dir,
    file.path(base_dir, "limma"),
    file.path(base_dir, "deseq2"),
    file.path(base_dir, "figures"),
    file.path(base_dir, "comparison")
  )
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

# Load data
load_data <- function(counts_file, design_file, duplicate_strategy = "sum") {
  cat("Loading data...\n")
  
  # Check files exist
  if (!file.exists(counts_file)) {
    stop(sprintf("Counts file not found: %s", counts_file))
  }
  if (!file.exists(design_file)) {
    stop(sprintf("Design file not found: %s", design_file))
  }
  
  # Validate duplicate strategy
  if (!duplicate_strategy %in% c("sum", "first", "unique")) {
    stop("duplicate_strategy must be 'sum', 'first', or 'unique'")
  }
  
  # Load counts - first without setting row names to check for duplicates
  vcat("  Reading counts file: ", counts_file, "\n")
  counts_raw <- read.table(counts_file, header = TRUE, check.names = FALSE, 
                           sep = "\t", stringsAsFactors = FALSE)
  
  # Check for duplicate gene IDs
  gene_ids <- counts_raw[, 1]
  if (any(duplicated(gene_ids))) {
    n_dup <- sum(duplicated(gene_ids))
    cat(sprintf("  Warning: Found %d duplicate gene IDs\n", n_dup))
    
    if (duplicate_strategy == "sum") {
      cat("  Summing counts for duplicate genes...\n")
      # Sum counts for duplicates
      counts_matrix <- as.matrix(counts_raw[, -1])
      rownames(counts_matrix) <- gene_ids
      
      # Aggregate by summing
      counts_summed <- aggregate(counts_matrix, by = list(gene_ids), FUN = sum)
      rownames(counts_summed) <- counts_summed[, 1]
      counts <- as.matrix(counts_summed[, -1])
      
      cat(sprintf("  Collapsed to %d unique genes\n", nrow(counts)))
      
    } else if (duplicate_strategy == "first") {
      cat("  Keeping first occurrence of duplicate genes...\n")
      keep_idx <- !duplicated(gene_ids)
      counts_raw <- counts_raw[keep_idx, ]
      gene_ids <- counts_raw[, 1]
      rownames(counts_raw) <- gene_ids
      counts <- as.matrix(counts_raw[, -1])
      
      cat(sprintf("  Kept %d unique genes\n", nrow(counts)))
      
    } else if (duplicate_strategy == "unique") {
      cat("  Making gene IDs unique by adding suffixes...\n")
      gene_ids <- make.unique(gene_ids)
      rownames(counts_raw) <- gene_ids
      counts <- as.matrix(counts_raw[, -1])
      
      cat(sprintf("  Created %d unique gene IDs\n", nrow(counts)))
    }
  } else {
    # No duplicates, proceed normally
    rownames(counts_raw) <- gene_ids
    counts <- as.matrix(counts_raw[, -1])
  }
  
  # Load design
  vcat("  Reading design file: ", design_file, "\n")
  design <- read.table(design_file, header = TRUE, row.names = 1, 
                      check.names = FALSE, sep = "\t")
  
  # Check that samples match
  common_samples <- intersect(colnames(counts), rownames(design))
  if (length(common_samples) == 0) {
    stop("No matching samples between counts and design files!")
  }
  
  if (length(common_samples) < ncol(counts)) {
    warning(sprintf("Only %d of %d samples in counts file found in design file",
                   length(common_samples), ncol(counts)))
  }
  
  counts <- counts[, common_samples]
  design <- design[common_samples, , drop = FALSE]
  
  cat(sprintf("  Loaded %d genes and %d samples\n", nrow(counts), ncol(counts)))
  
  return(list(counts = counts, design = design))
}

# Filter low-expressed genes
filter_genes <- function(counts, min_count, min_samples) {
  cat(sprintf("Filtering genes (min %d counts in at least %d samples)...\n", 
              min_count, min_samples))
  
  keep <- rowSums(counts >= min_count) >= min_samples
  counts_filtered <- counts[keep, ]
  
  cat(sprintf("  Kept %d of %d genes (%.1f%%)\n", 
              sum(keep), length(keep), 100 * sum(keep) / length(keep)))
  
  return(counts_filtered)
}

# Create design matrix for limma
create_limma_design <- function(design, time_column, time_as_numeric, 
                               spline_df, covariates) {
  require(splines)
  
  # Get time variable
  time_var <- design[[time_column]]
  
  if (time_as_numeric) {
    # Convert to numeric if not already
    time_var <- as.numeric(as.character(time_var))
    
    if (!is.null(spline_df)) {
      # Natural spline
      vcat(sprintf("  Using natural splines with df=%d\n", spline_df))
      time_matrix <- ns(time_var, df = spline_df)
      colnames(time_matrix) <- paste0("time_ns", 1:spline_df)
    } else {
      # Linear + quadratic terms
      vcat("  Using polynomial terms (linear + quadratic)\n")
      time_matrix <- cbind(time_var, time_var^2)
      colnames(time_matrix) <- c("time_linear", "time_quadratic")
    }
  } else {
    # Categorical time
    time_var <- factor(time_var)
    vcat(sprintf("  Using categorical time with %d levels\n", nlevels(time_var)))
    time_matrix <- model.matrix(~ 0 + time_var)
    colnames(time_matrix) <- levels(time_var)
  }
  
  # Add covariates if specified
  if (length(covariates) > 0) {
    vcat(sprintf("  Adding %d covariate(s): %s\n", 
                length(covariates), paste(covariates, collapse=", ")))
    
    covariate_matrices <- lapply(covariates, function(cov) {
      if (!cov %in% colnames(design)) {
        stop(sprintf("Covariate '%s' not found in design file!", cov))
      }
      var <- design[[cov]]
      if (is.numeric(var)) {
        mat <- matrix(var, ncol = 1)
        colnames(mat) <- cov
      } else {
        mat <- model.matrix(~ 0 + factor(var))
        colnames(mat) <- paste0(cov, "_", levels(factor(var)))
      }
      return(mat)
    })
    design_matrix <- cbind(time_matrix, do.call(cbind, covariate_matrices))
  } else {
    design_matrix <- time_matrix
  }
  
  rownames(design_matrix) <- rownames(design)
  
  return(design_matrix)
}

# Create interaction design matrix for limma (time × group)
create_interaction_design_limma <- function(design, params) {
  require(splines)
  
  # Get group variable as factor with specified levels
  group_var <- factor(design[[params$group_column]], 
                     levels = params$compare_groups)
  
  # Get time variable
  time_var <- design[[params$time_column]]
  
  if (params$time_as_numeric) {
    time_var <- as.numeric(as.character(time_var))
    
    if (!is.null(params$spline_df)) {
      vcat(sprintf("  Using natural splines with df=%d for interaction model\n", 
                  params$spline_df))
      # Create time basis
      time_basis <- ns(time_var, df = params$spline_df)
      
      # Create design with interaction: group + time + group:time
      design_matrix <- model.matrix(~ group_var * time_basis)
      
      # Rename columns for clarity
      colnames(design_matrix) <- gsub("group_var", "group", colnames(design_matrix))
      colnames(design_matrix) <- gsub("time_basis", "time_ns", colnames(design_matrix))
    } else {
      vcat("  Using polynomial terms (linear + quadratic) for interaction model\n")
      time_linear <- time_var
      time_quad <- time_var^2
      
      design_matrix <- model.matrix(~ group_var * (time_linear + time_quad))
      colnames(design_matrix) <- gsub("group_var", "group", colnames(design_matrix))
    }
  } else {
    # Categorical time
    time_var <- factor(time_var)
    vcat(sprintf("  Using categorical time (%d levels) for interaction model\n", 
                nlevels(time_var)))
    
    design_matrix <- model.matrix(~ group_var * time_var)
    colnames(design_matrix) <- gsub("group_var", "group", colnames(design_matrix))
    colnames(design_matrix) <- gsub("time_var", "", colnames(design_matrix))
  }
  
  rownames(design_matrix) <- rownames(design)
  
  return(design_matrix)
}

################################################################################
# LIMMA-VOOM ANALYSIS
################################################################################

run_limma_analysis <- function(counts, design, params) {
  cat("\n", rep("=", 70), "\n", sep = "")
  if (params$test_interaction) {
    cat("RUNNING LIMMA-VOOM ANALYSIS - DIFFERENTIAL TIMECOURSE\n")
    cat("Testing for Time × Group Interaction\n")
  } else {
    cat("RUNNING LIMMA-VOOM ANALYSIS\n")
  }
  cat(rep("=", 70), "\n", sep = "")
  
  # Filter design to include only the two groups being compared
  if (params$test_interaction) {
    group_var <- design[[params$group_column]]
    keep_samples <- group_var %in% params$compare_groups
    
    if (sum(keep_samples) == 0) {
      stop(sprintf("No samples found for groups: %s", 
                  paste(params$compare_groups, collapse=", ")))
    }
    
    design <- design[keep_samples, , drop = FALSE]
    counts <- counts[, keep_samples, drop = FALSE]
    
    cat(sprintf("Comparing groups: %s vs %s\n", 
                params$compare_groups[1], params$compare_groups[2]))
    cat(sprintf("  %s: %d samples\n", params$compare_groups[1], 
                sum(group_var[keep_samples] == params$compare_groups[1])))
    cat(sprintf("  %s: %d samples\n", params$compare_groups[2],
                sum(group_var[keep_samples] == params$compare_groups[2])))
  }
  
  # Create DGEList object
  dge <- DGEList(counts = counts)
  
  # Calculate normalization factors
  vcat("  Calculating normalization factors...\n")
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  if (params$test_interaction) {
    design_matrix <- create_interaction_design_limma(design, params)
  } else {
    design_matrix <- create_limma_design(design, params$time_column, 
                                         params$time_as_numeric, params$spline_df, 
                                         params$covariates)
  }
  
  cat("\nDesign matrix dimensions:", nrow(design_matrix), "x", ncol(design_matrix), "\n")
  vcat("Design matrix:\n")
  if (params$verbose) print(head(design_matrix))
  
  # Voom transformation
  vcat("  Running voom transformation...\n")
  v <- voom(dge, design_matrix, plot = FALSE)
  
  # Save voom plot
  if (params$generate_plots) {
    plot_file <- if (params$test_interaction) {
      file.path(params$output_dir, "figures", "limma_voom_plot_interaction.pdf")
    } else {
      file.path(params$output_dir, "figures", "limma_voom_plot.pdf")
    }
    pdf(plot_file, width = 8, height = 6)
    v <- voom(dge, design_matrix, plot = TRUE)
    dev.off()
  }
  
  # Fit linear model
  vcat("  Fitting linear model...\n")
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)
  
  # Test appropriate coefficients
  if (params$test_interaction) {
    # Test interaction terms
    if (params$time_as_numeric) {
      if (!is.null(params$spline_df)) {
        interaction_coefs <- paste0("group", params$compare_groups[2], ":time_ns", 1:params$spline_df)
      } else {
        interaction_coefs <- c(paste0("group", params$compare_groups[2], ":time_linear"),
                              paste0("group", params$compare_groups[2], ":time_quad"))
      }
    } else {
      time_levels <- levels(factor(design[[params$time_column]]))
      interaction_coefs <- paste0("group", params$compare_groups[2], ":", time_levels[-1])
    }
    
    vcat("  Testing interaction coefficients:", paste(interaction_coefs, collapse=", "), "\n")
    results <- topTable(fit, coef = interaction_coefs, number = Inf, sort.by = "F")
    
    result_file <- file.path(params$output_dir, "limma", 
                            "limma_interaction_results.txt")
  } else {
    # Test for time effect (all time coefficients)
    if (params$time_as_numeric) {
      if (!is.null(params$spline_df)) {
        time_coefs <- paste0("time_ns", 1:params$spline_df)
      } else {
        time_coefs <- c("time_linear", "time_quadratic")
      }
    } else {
      time_coefs <- levels(factor(design[[params$time_column]]))
    }
    
    vcat("  Testing time effect with coefficients:", paste(time_coefs, collapse=", "), "\n")
    results <- topTable(fit, coef = time_coefs, number = Inf, sort.by = "F")
    
    result_file <- file.path(params$output_dir, "limma", 
                            "limma_time_effect_results.txt")
  }
  
  # Save results
  write.table(results, file = result_file,
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  
  cat(sprintf("\nFound %d genes with FDR < %.2f\n", 
              sum(results$adj.P.Val < params$fdr_threshold), params$fdr_threshold))
  
  # Return results
  return(list(
    fit = fit,
    results = results,
    voom = v,
    design_matrix = design_matrix,
    design = design,
    counts = counts
  ))
}

################################################################################
# DESEQ2 ANALYSIS
################################################################################

run_deseq2_analysis <- function(counts, design, params) {
  cat("\n", rep("=", 70), "\n", sep = "")
  if (params$test_interaction) {
    cat("RUNNING DESEQ2 ANALYSIS - DIFFERENTIAL TIMECOURSE\n")
    cat("Testing for Time × Group Interaction\n")
  } else {
    cat("RUNNING DESEQ2 ANALYSIS\n")
  }
  cat(rep("=", 70), "\n", sep = "")
  
  # Filter design to include only the two groups being compared
  if (params$test_interaction) {
    group_var <- design[[params$group_column]]
    keep_samples <- group_var %in% params$compare_groups
    
    if (sum(keep_samples) == 0) {
      stop(sprintf("No samples found for groups: %s", 
                  paste(params$compare_groups, collapse=", ")))
    }
    
    design <- design[keep_samples, , drop = FALSE]
    counts <- counts[, keep_samples, drop = FALSE]
    
    cat(sprintf("Comparing groups: %s vs %s\n", 
                params$compare_groups[1], params$compare_groups[2]))
    cat(sprintf("  %s: %d samples\n", params$compare_groups[1], 
                sum(group_var[keep_samples] == params$compare_groups[1])))
    cat(sprintf("  %s: %d samples\n", params$compare_groups[2],
                sum(group_var[keep_samples] == params$compare_groups[2])))
  }
  
  # Check if counts are integers (DESeq2 requirement)
  has_decimals <- any(counts != floor(counts))
  
  if (has_decimals) {
    cat("  Note: Counts contain non-integer values (common with RSEM/Salmon)\n")
    cat("  Rounding counts to integers for DESeq2...\n")
    counts <- round(counts)
    vcat("  Min rounded count: ", min(counts), "\n")
    vcat("  Max rounded count: ", max(counts), "\n")
  }
  
  # Prepare design formula
  time_var <- design[[params$time_column]]
  
  if (params$test_interaction) {
    # Interaction model
    group_var <- factor(design[[params$group_column]], 
                       levels = params$compare_groups)
    design$Group <- group_var
    
    if (params$time_as_numeric) {
      time_var <- as.numeric(as.character(time_var))
      design[[params$time_column]] <- time_var
      
      if (!is.null(params$spline_df)) {
        require(splines)
        vcat(sprintf("  Using natural splines with df=%d for interaction model\n", 
                    params$spline_df))
        formula_str <- paste0("~ Group * ns(", params$time_column, ", df = ", params$spline_df, ")")
        reduced_formula_str <- paste0("~ Group + ns(", params$time_column, ", df = ", params$spline_df, ")")
      } else {
        vcat("  Using polynomial terms (linear + quadratic) for interaction model\n")
        formula_str <- paste0("~ Group * (", params$time_column, " + I(", params$time_column, "^2))")
        reduced_formula_str <- paste0("~ Group + ", params$time_column, " + I(", params$time_column, "^2)")
      }
    } else {
      design[[params$time_column]] <- factor(time_var)
      vcat(sprintf("  Using categorical time (%d levels) for interaction model\n", 
                  nlevels(factor(time_var))))
      formula_str <- paste0("~ Group * ", params$time_column)
      reduced_formula_str <- paste0("~ Group + ", params$time_column)
    }
    
    result_file <- file.path(params$output_dir, "deseq2", "deseq2_interaction_results.txt")
    
  } else {
    # Time effect model
    if (params$time_as_numeric) {
      time_var <- as.numeric(as.character(time_var))
      design[[params$time_column]] <- time_var
      
      if (!is.null(params$spline_df)) {
        require(splines)
        vcat(sprintf("  Using natural splines with df=%d\n", params$spline_df))
        design$time_ns <- ns(time_var, df = params$spline_df)
        formula_str <- paste0("~ ns(", params$time_column, ", df = ", params$spline_df, ")")
      } else {
        vcat("  Using polynomial terms (linear + quadratic)\n")
        formula_str <- paste0("~ ", params$time_column, " + I(", params$time_column, "^2)")
      }
    } else {
      design[[params$time_column]] <- factor(time_var)
      vcat(sprintf("  Using categorical time with %d levels\n", 
                  nlevels(factor(time_var))))
      formula_str <- paste0("~ ", params$time_column)
    }
    
    # Add covariates
    if (length(params$covariates) > 0) {
      vcat(sprintf("  Adding %d covariate(s): %s\n", 
                  length(params$covariates), paste(params$covariates, collapse=", ")))
      formula_str <- paste0(formula_str, " + ", paste(params$covariates, collapse = " + "))
    }
    
    reduced_formula_str <- "~ 1"
    result_file <- file.path(params$output_dir, "deseq2", "deseq2_time_effect_results.txt")
  }
  
  cat("\nDESeq2 formula:", formula_str, "\n")
  if (params$test_interaction) {
    cat("Reduced formula:", reduced_formula_str, "\n")
  }
  
  # Create DESeq2 dataset
  vcat("  Creating DESeqDataSet...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = design,
    design = as.formula(formula_str)
  )
  
  # Run DESeq2
  vcat("  Running DESeq2 (this may take a while)...\n")
  dds <- DESeq(dds, test = "LRT", reduced = as.formula(reduced_formula_str), 
               quiet = !params$verbose)
  
  # Get results
  vcat("  Extracting results...\n")
  res <- results(dds)
  res <- res[order(res$padj), ]
  
  # Convert to data frame
  results_df <- as.data.frame(res)
  
  # Save results
  write.table(results_df, file = result_file,
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  
  cat(sprintf("\nFound %d genes with FDR < %.2f\n",
              sum(results_df$padj < params$fdr_threshold, na.rm = TRUE), 
              params$fdr_threshold))
  
  # Normalized counts for visualization
  norm_counts <- counts(dds, normalized = TRUE)
  
  return(list(
    dds = dds,
    results = results_df,
    norm_counts = norm_counts,
    design = design,
    counts = counts
  ))
}

################################################################################
# CLUSTERING AND PATTERN ANALYSIS FUNCTIONS
################################################################################

# Perform temporal pattern clustering
cluster_temporal_patterns <- function(expr_data, design, time_column, 
                                     n_clusters, method = "kmeans") {
  vcat("  Performing temporal pattern clustering...\n")
  vcat(sprintf("    Method: %s\n", method))
  vcat(sprintf("    Number of clusters: %d\n", n_clusters))
  
  # Get time points
  time_points <- sort(unique(design[[time_column]]))
  
  # Calculate mean expression at each timepoint for each gene
  mean_profiles <- t(sapply(rownames(expr_data), function(gene) {
    sapply(time_points, function(tp) {
      samples <- rownames(design)[design[[time_column]] == tp]
      mean(expr_data[gene, samples])
    })
  }))
  colnames(mean_profiles) <- paste0("T", time_points)
  
  # Scale profiles (each gene centered and scaled)
  scaled_profiles <- t(scale(t(mean_profiles)))
  
  # Remove genes with NA (e.g., no variation)
  valid_genes <- complete.cases(scaled_profiles)
  scaled_profiles <- scaled_profiles[valid_genes, ]
  
  # Perform clustering
  if (method == "kmeans") {
    set.seed(123)  # For reproducibility
    clustering <- kmeans(scaled_profiles, centers = n_clusters, nstart = 25)
    clusters <- clustering$cluster
    centers <- clustering$centers
    
  } else if (method == "hierarchical") {
    hc <- hclust(dist(scaled_profiles), method = "ward.D2")
    clusters <- cutree(hc, k = n_clusters)
    
    # Calculate centers
    centers <- t(sapply(1:n_clusters, function(k) {
      if (sum(clusters == k) > 1) {
        colMeans(scaled_profiles[clusters == k, , drop = FALSE])
      } else {
        scaled_profiles[clusters == k, ]
      }
    }))
    
  } else if (method == "fuzzy") {
    # Fuzzy c-means clustering
    if (requireNamespace("e1071", quietly = TRUE)) {
      library(e1071)
      fcm_result <- cmeans(scaled_profiles, centers = n_clusters)
      clusters <- fcm_result$cluster
      centers <- fcm_result$centers
    } else {
      cat("  Warning: e1071 package not available, falling back to kmeans\n")
      set.seed(123)
      clustering <- kmeans(scaled_profiles, centers = n_clusters, nstart = 25)
      clusters <- clustering$cluster
      centers <- clustering$centers
    }
  }
  
  # Classify pattern types
  pattern_types <- classify_patterns(centers, time_points)
  
  return(list(
    clusters = clusters,
    centers = centers,
    scaled_profiles = scaled_profiles,
    mean_profiles = mean_profiles[valid_genes, ],
    time_points = time_points,
    pattern_types = pattern_types,
    genes_per_cluster = table(clusters)
  ))
}

# Classify temporal patterns
classify_patterns <- function(centers, time_points) {
  patterns <- character(nrow(centers))
  
  for (i in 1:nrow(centers)) {
    profile <- centers[i, ]
    
    # Find peak and trough
    peak_idx <- which.max(profile)
    trough_idx <- which.min(profile)
    
    # Calculate overall trend
    start_val <- profile[1]
    end_val <- profile[length(profile)]
    mid_val <- profile[ceiling(length(profile)/2)]
    
    # Classify pattern
    if (abs(profile[peak_idx] - profile[trough_idx]) < 0.5) {
      patterns[i] <- "Flat"
      
    } else if (end_val > start_val + 0.5) {
      if (peak_idx == length(profile)) {
        patterns[i] <- "Sustained Up"
      } else {
        patterns[i] <- "Transient Up"
      }
      
    } else if (end_val < start_val - 0.5) {
      if (trough_idx == length(profile)) {
        patterns[i] <- "Sustained Down"
      } else {
        patterns[i] <- "Transient Down"
      }
      
    } else if (peak_idx <= 2) {
      patterns[i] <- "Early Peak"
      
    } else if (peak_idx >= length(profile) - 1) {
      patterns[i] <- "Late Peak"
      
    } else if (peak_idx > 2 && peak_idx < length(profile) - 1) {
      patterns[i] <- "Mid Peak"
      
    } else {
      patterns[i] <- "Complex"
    }
  }
  
  return(patterns)
}

# Plot cluster profiles
plot_cluster_profiles <- function(cluster_result, output_file, method_name) {
  vcat("  Creating cluster profile plot...\n")
  
  centers <- cluster_result$centers
  time_points <- cluster_result$time_points
  pattern_types <- cluster_result$pattern_types
  genes_per_cluster <- cluster_result$genes_per_cluster
  
  # Prepare data
  plot_data <- data.frame(
    Cluster = rep(1:nrow(centers), each = ncol(centers)),
    Time = rep(time_points, nrow(centers)),
    Expression = as.vector(t(centers)),
    Pattern = rep(pattern_types, each = ncol(centers)),
    N_genes = rep(genes_per_cluster, each = ncol(centers))
  )
  
  plot_data$Cluster_label <- paste0("Cluster ", plot_data$Cluster, ": ", 
                                    plot_data$Pattern, " (n=", plot_data$N_genes, ")")
  
  # Plot
  p <- ggplot(plot_data, aes(x = Time, y = Expression, color = factor(Cluster))) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
    facet_wrap(~ Cluster_label, ncol = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste0("Temporal Pattern Clusters - ", method_name),
         x = "Time", y = "Scaled Expression (Z-score)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  
  ggsave(output_file, p, width = 12, height = ceiling(nrow(centers)/3) * 3)
  
  return(p)
}

# Plot cluster heatmap with all genes
plot_cluster_heatmap <- function(cluster_result, output_file, method_name) {
  vcat("  Creating cluster heatmap...\n")
  
  # Order genes by cluster
  genes_ordered <- names(sort(cluster_result$clusters))
  
  # Get scaled profiles
  mat <- cluster_result$scaled_profiles[genes_ordered, ]
  
  # Create annotation
  cluster_annotation <- data.frame(
    Cluster = factor(cluster_result$clusters[genes_ordered]),
    row.names = genes_ordered
  )
  
  pattern_annotation <- data.frame(
    Pattern = cluster_result$pattern_types[cluster_result$clusters[genes_ordered]],
    row.names = genes_ordered
  )
  
  annotation_row <- cbind(cluster_annotation, pattern_annotation)
  
  # Colors
  n_clusters <- length(unique(cluster_result$clusters))
  cluster_colors <- colorRampPalette(brewer.pal(min(9, n_clusters), "Set1"))(n_clusters)
  names(cluster_colors) <- as.character(1:n_clusters)
  
  pattern_types <- unique(cluster_result$pattern_types)
  pattern_colors <- colorRampPalette(brewer.pal(min(9, length(pattern_types)), "Set2"))(length(pattern_types))
  names(pattern_colors) <- pattern_types
  
  annotation_colors <- list(
    Cluster = cluster_colors,
    Pattern = pattern_colors
  )
  
  heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # Limit number of genes for display if too many
  if (nrow(mat) > 200) {
    cat(sprintf("  Note: Showing %d representative genes (total: %d)\n", 200, nrow(mat)))
    # Sample genes from each cluster proportionally
    genes_to_show <- unlist(lapply(1:n_clusters, function(k) {
      cluster_genes <- names(cluster_result$clusters)[cluster_result$clusters == k]
      n_show <- min(ceiling(200 * length(cluster_genes) / nrow(mat)), length(cluster_genes))
      sample(cluster_genes, n_show)
    }))
    mat <- mat[genes_to_show, ]
    annotation_row <- annotation_row[genes_to_show, , drop = FALSE]
  }
  
  # Plot
  pdf(output_file, width = 10, height = 14)
  pheatmap(mat,
           color = heatmap_colors,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main = paste0("Clustered Temporal Patterns - ", method_name),
           fontsize_col = 10,
           breaks = seq(-2, 2, length.out = 101))
  dev.off()
}

# Create cluster summary table
create_cluster_summary <- function(cluster_result, output_file) {
  vcat("  Creating cluster summary table...\n")
  
  summary_df <- data.frame(
    Cluster = 1:nrow(cluster_result$centers),
    Pattern = cluster_result$pattern_types,
    N_genes = as.vector(cluster_result$genes_per_cluster),
    Percent = round(100 * as.vector(cluster_result$genes_per_cluster) / 
                     sum(cluster_result$genes_per_cluster), 1)
  )
  
  # Add profile characteristics
  for (i in 1:nrow(cluster_result$centers)) {
    profile <- cluster_result$centers[i, ]
    summary_df$Peak_time[i] <- cluster_result$time_points[which.max(profile)]
    summary_df$Trough_time[i] <- cluster_result$time_points[which.min(profile)]
    summary_df$Range[i] <- round(max(profile) - min(profile), 2)
  }
  
  write.table(summary_df, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE)
  
  # Print to console
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("CLUSTER SUMMARY\n")
  cat(rep("=", 70), "\n", sep = "")
  print(summary_df, row.names = FALSE)
  cat("\n")
  
  return(summary_df)
}

# Export gene lists by cluster
export_cluster_gene_lists <- function(cluster_result, results_df, output_dir, method_name = "limma") {
  vcat("  Exporting gene lists by cluster...\n")
  
  cluster_dir <- file.path(output_dir, "clusters")
  if (!dir.exists(cluster_dir)) dir.create(cluster_dir, recursive = TRUE)
  
  # Determine which p-value column to use
  if ("adj.P.Val" %in% colnames(results_df)) {
    pval_col <- "adj.P.Val"
  } else if ("padj" %in% colnames(results_df)) {
    pval_col <- "padj"
  } else {
    pval_col <- NULL
  }
  
  for (i in 1:length(unique(cluster_result$clusters))) {
    cluster_genes <- names(cluster_result$clusters)[cluster_result$clusters == i]
    
    # Only keep genes that exist in results
    cluster_genes <- cluster_genes[cluster_genes %in% rownames(results_df)]
    
    if (length(cluster_genes) == 0) {
      vcat(sprintf("    Warning: No genes found for cluster %d\n", i))
      next
    }
    
    # Get results for these genes
    cluster_results <- results_df[cluster_genes, , drop = FALSE]
    
    # Add profile information first
    cluster_results$Cluster <- i
    cluster_results$Pattern <- cluster_result$pattern_types[i]
    
    # Order by p-value if available
    if (!is.null(pval_col) && pval_col %in% colnames(cluster_results)) {
      pvals <- cluster_results[[pval_col]]
      if (!all(is.na(pvals))) {
        order_idx <- order(pvals, na.last = TRUE)
        cluster_results <- cluster_results[order_idx, , drop = FALSE]
      }
    }
    
    # Save
    outfile <- file.path(cluster_dir, 
                        sprintf("cluster_%d_%s_genes.txt", i, 
                               gsub(" ", "_", cluster_result$pattern_types[i])))
    
    write.table(cluster_results, file = outfile, sep = "\t", quote = FALSE,
                row.names = TRUE, col.names = NA)
    
    vcat(sprintf("    Cluster %d (%s): %d genes\n", i, 
                cluster_result$pattern_types[i], length(cluster_genes)))
  }
  
  # Also create a master file with all genes and cluster assignments
  all_genes_clusters <- data.frame(
    Gene = names(cluster_result$clusters),
    Cluster = cluster_result$clusters,
    Pattern = cluster_result$pattern_types[cluster_result$clusters],
    stringsAsFactors = FALSE
  )
  
  # Only merge genes that exist in results
  genes_in_results <- all_genes_clusters$Gene[all_genes_clusters$Gene %in% rownames(results_df)]
  
  if (length(genes_in_results) > 0) {
    # Create a temporary dataframe with rownames as a column for merging
    results_with_genes <- results_df[genes_in_results, , drop = FALSE]
    results_with_genes$Gene <- rownames(results_with_genes)
    
    # Merge
    all_genes_clusters <- merge(all_genes_clusters, 
                                results_with_genes, 
                                by = "Gene",
                                all.x = TRUE)
    
    # Order by cluster and p-value
    if (!is.null(pval_col) && pval_col %in% colnames(all_genes_clusters)) {
      pvals <- all_genes_clusters[[pval_col]]
      order_idx <- order(all_genes_clusters$Cluster, pvals, na.last = TRUE)
      all_genes_clusters <- all_genes_clusters[order_idx, ]
    } else {
      all_genes_clusters <- all_genes_clusters[order(all_genes_clusters$Cluster), ]
    }
    
    write.table(all_genes_clusters, 
                file = file.path(cluster_dir, "all_genes_with_clusters.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  cat(sprintf("\nGene lists exported to: %s\n", cluster_dir))
}

# Plot top genes per cluster (profiles)
plot_cluster_top_genes_profiles <- function(cluster_result, expr_data, results_df, 
                                            design, time_column, output_dir, 
                                            method_name, n_top = 20) {
  vcat("  Creating per-cluster top genes profile plots...\n")
  
  # Determine p-value column
  if ("adj.P.Val" %in% colnames(results_df)) {
    pval_col <- "adj.P.Val"
  } else if ("padj" %in% colnames(results_df)) {
    pval_col <- "padj"
  } else {
    return()
  }
  
  # Create directory for cluster plots
  cluster_plot_dir <- file.path(output_dir, "figures", "cluster_top_genes")
  if (!dir.exists(cluster_plot_dir)) dir.create(cluster_plot_dir, recursive = TRUE)
  
  for (i in 1:length(unique(cluster_result$clusters))) {
    cluster_genes <- names(cluster_result$clusters)[cluster_result$clusters == i]
    
    # Get top genes from this cluster by significance
    cluster_genes_in_results <- cluster_genes[cluster_genes %in% rownames(results_df)]
    if (length(cluster_genes_in_results) == 0) next
    
    cluster_results <- results_df[cluster_genes_in_results, , drop = FALSE]
    cluster_results <- cluster_results[order(cluster_results[[pval_col]], na.last = TRUE), ]
    
    top_genes <- rownames(cluster_results)[1:min(n_top, nrow(cluster_results))]
    
    # Create profile plot
    outfile <- file.path(cluster_plot_dir,
                        sprintf("cluster_%d_%s_top20_profiles.pdf", i,
                               gsub(" ", "_", cluster_result$pattern_types[i])))
    
    plot_timecourse_profiles(
      expr_data, design, time_column, top_genes, outfile,
      sprintf("%s - Cluster %d (%s)", method_name, i, cluster_result$pattern_types[i])
    )
    
    vcat(sprintf("    Cluster %d (%s): %d genes\n", i, 
                cluster_result$pattern_types[i], length(top_genes)))
  }
  
  cat(sprintf("\nPer-cluster profile plots saved to: %s\n", cluster_plot_dir))
}

# Plot top genes per cluster (heatmaps)
plot_cluster_top_genes_heatmaps <- function(cluster_result, expr_data, results_df,
                                            design, time_column, output_dir,
                                            method_name, n_top = 50, group_column = NULL) {
  vcat("  Creating per-cluster top genes heatmaps...\n")
  
  # Determine p-value column
  if ("adj.P.Val" %in% colnames(results_df)) {
    pval_col <- "adj.P.Val"
  } else if ("padj" %in% colnames(results_df)) {
    pval_col <- "padj"
  } else {
    return()
  }
  
  # Create directory for cluster heatmaps
  cluster_heatmap_dir <- file.path(output_dir, "figures", "cluster_top_genes")
  if (!dir.exists(cluster_heatmap_dir)) dir.create(cluster_heatmap_dir, recursive = TRUE)
  
  for (i in 1:length(unique(cluster_result$clusters))) {
    cluster_genes <- names(cluster_result$clusters)[cluster_result$clusters == i]
    
    # Get top genes from this cluster by significance
    cluster_genes_in_results <- cluster_genes[cluster_genes %in% rownames(results_df)]
    if (length(cluster_genes_in_results) == 0) next
    
    cluster_results <- results_df[cluster_genes_in_results, , drop = FALSE]
    cluster_results <- cluster_results[order(cluster_results[[pval_col]], na.last = TRUE), ]
    
    top_genes <- rownames(cluster_results)[1:min(n_top, nrow(cluster_results))]
    
    # Create heatmap
    outfile <- file.path(cluster_heatmap_dir,
                        sprintf("cluster_%d_%s_top50_heatmap.pdf", i,
                               gsub(" ", "_", cluster_result$pattern_types[i])))
    
    plot_heatmap(
      expr_data, design, time_column, top_genes, outfile,
      sprintf("%s - Cluster %d (%s)", method_name, i, cluster_result$pattern_types[i]),
      group_column = group_column
    )
    
    vcat(sprintf("    Cluster %d (%s): %d genes\n", i,
                cluster_result$pattern_types[i], length(top_genes)))
  }
  
  cat(sprintf("\nPer-cluster heatmaps saved to: %s\n", cluster_heatmap_dir))
}

################################################################################
# PATHWAY ENRICHMENT FUNCTIONS
################################################################################

# Load GMT file
load_gmt_file <- function(gmt_file) {
  vcat("  Loading GMT file:", gmt_file, "\n")
  
  if (!file.exists(gmt_file)) {
    stop(sprintf("GMT file not found: %s", gmt_file))
  }
  
  # Read GMT file
  gmt_lines <- readLines(gmt_file)
  
  # Parse GMT format: pathway_name \t description \t gene1 \t gene2 \t ...
  pathways <- list()
  pathway_names <- character()
  pathway_descriptions <- character()
  
  for (line in gmt_lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 3) {
      pathway_name <- fields[1]
      description <- fields[2]
      genes <- fields[3:length(fields)]
      genes <- genes[genes != ""]  # Remove empty strings
      
      pathways[[pathway_name]] <- genes
      pathway_names <- c(pathway_names, pathway_name)
      pathway_descriptions <- c(pathway_descriptions, description)
    }
  }
  
  vcat(sprintf("  Loaded %d pathways\n", length(pathways)))
  
  return(list(
    pathways = pathways,
    names = pathway_names,
    descriptions = pathway_descriptions
  ))
}

# Detect gene ID type and convert if needed
detect_and_convert_gene_ids <- function(gene_list, organism, target_type = "SYMBOL") {
  vcat(sprintf("  Detecting gene ID type for %d genes...\n", length(gene_list)))
  
  # Sample first 100 genes to check type
  sample_genes <- head(gene_list, 100)
  
  # Check if IDs look like Entrez (all numeric)
  all_numeric <- all(grepl("^[0-9]+$", sample_genes))
  
  # Check if IDs look like Ensembl
  has_ensembl <- any(grepl("^ENS", sample_genes))
  
  # Check if IDs look like symbols (contain letters)
  has_letters <- any(grepl("[A-Za-z]", sample_genes))
  
  if (all_numeric) {
    id_type <- "ENTREZID"
    vcat("  Detected: Entrez IDs\n")
  } else if (has_ensembl) {
    id_type <- "ENSEMBL"
    vcat("  Detected: Ensembl IDs\n")
  } else if (has_letters) {
    id_type <- "SYMBOL"
    vcat("  Detected: Gene Symbols\n")
  } else {
    id_type <- "SYMBOL"
    vcat("  Assuming: Gene Symbols\n")
  }
  
  # If already the target type, return as is
  if (id_type == target_type) {
    vcat(sprintf("  IDs are already %s type\n", target_type))
    return(list(genes = as.character(gene_list), type = id_type))
  }
  
  # For Reactome with Entrez IDs, just ensure they're character
  if (target_type == "ENTREZID" && id_type == "ENTREZID") {
    return(list(genes = as.character(gene_list), type = id_type))
  }
  
  # Attempt conversion if needed
  if (id_type != target_type) {
    vcat(sprintf("  Converting from %s to %s...\n", id_type, target_type))
    
    # Set organism database
    if (organism == "human") {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        cat("  Warning: org.Hs.eg.db required for conversion. Install with:\n")
        cat("  BiocManager::install('org.Hs.eg.db')\n")
        cat("  Using IDs as-is (may cause issues).\n")
        return(list(genes = as.character(gene_list), type = id_type))
      }
      library(org.Hs.eg.db)
      orgdb <- org.Hs.eg.db
    } else {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        cat("  Warning: org.Mm.eg.db required for conversion. Install with:\n")
        cat("  BiocManager::install('org.Mm.eg.db')\n")
        cat("  Using IDs as-is (may cause issues).\n")
        return(list(genes = as.character(gene_list), type = id_type))
      }
      library(org.Mm.eg.db)
      orgdb <- org.Mm.eg.db
    }
    
    # Convert IDs
    if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
      library(AnnotationDbi)
      converted <- tryCatch({
        AnnotationDbi::mapIds(orgdb, 
                             keys = as.character(gene_list),
                             column = target_type,
                             keytype = id_type,
                             multiVals = "first")
      }, error = function(e) {
        cat("  Error converting IDs:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(converted)) {
        converted <- converted[!is.na(converted)]
        vcat(sprintf("  Successfully converted %d/%d genes\n", length(converted), length(gene_list)))
        if (length(converted) < length(gene_list) * 0.5) {
          cat(sprintf("  Warning: Only %d%% of genes converted successfully\n",
                     round(100 * length(converted) / length(gene_list))))
        }
        return(list(genes = as.character(converted), type = target_type))
      }
    } else {
      cat("  Warning: AnnotationDbi package required for ID conversion\n")
      cat("  Install with: BiocManager::install('AnnotationDbi')\n")
    }
  }
  
  # Return original if conversion failed
  return(list(genes = as.character(gene_list), type = id_type))
}

# Run enrichment analysis
run_enrichment_analysis <- function(gene_list, background_genes, params) {
  
  # Load required packages
  if (params$enrichment_db %in% c("reactome", "gobp")) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      cat("  Error: clusterProfiler package required for Reactome/GOBP enrichment\n")
      cat("  Install with: BiocManager::install('clusterProfiler')\n")
      return(NULL)
    }
    library(clusterProfiler)
  }
  
  if (params$enrichment_db == "reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      cat("  Error: ReactomePA package required for Reactome enrichment\n")
      cat("  Install with: BiocManager::install('ReactomePA')\n")
      return(NULL)
    }
    library(ReactomePA)
    
    # Reactome requires ENTREZ IDs
    gene_data <- detect_and_convert_gene_ids(gene_list, params$organism, "ENTREZID")
    bg_data <- detect_and_convert_gene_ids(background_genes, params$organism, "ENTREZID")
    
    if (length(gene_data$genes) == 0) {
      cat("  Error: No genes remain after ID conversion\n")
      return(NULL)
    }
    
    # Set organism
    org_code <- if (params$organism == "human") "human" else "mouse"
    
    vcat(sprintf("  Running Reactome enrichment (%s) with %d Entrez IDs...\n", 
                params$organism, length(gene_data$genes)))
    
    result <- tryCatch({
      enrichPathway(
        gene = gene_data$genes,
        organism = org_code,
        pvalueCutoff = 1.0,  # Get all pathways, will filter for display later
        pAdjustMethod = "BH",
        universe = bg_data$genes,
        readable = TRUE
      )
    }, error = function(e) {
      cat("  Error in Reactome enrichment:", e$message, "\n")
      return(NULL)
    })
    
  } else if (params$enrichment_db == "gobp") {
    
    # Set organism database
    if (params$organism == "human") {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        cat("  Error: org.Hs.eg.db package required for human GO enrichment\n")
        cat("  Install with: BiocManager::install('org.Hs.eg.db')\n")
        return(NULL)
      }
      library(org.Hs.eg.db)
      orgdb <- org.Hs.eg.db
    } else {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        cat("  Error: org.Mm.eg.db package required for mouse GO enrichment\n")
        cat("  Install with: BiocManager::install('org.Mm.eg.db')\n")
        return(NULL)
      }
      library(org.Mm.eg.db)
      orgdb <- org.Mm.eg.db
    }
    
    # Detect ID type and set keyType appropriately
    gene_data <- detect_and_convert_gene_ids(gene_list, params$organism, "SYMBOL")
    bg_data <- detect_and_convert_gene_ids(background_genes, params$organism, "SYMBOL")
    
    # Determine keyType based on what we have
    keytype_to_use <- gene_data$type
    
    vcat(sprintf("  Running GO BP enrichment (%s) with %d genes (%s)...\n", 
                params$organism, length(gene_data$genes), keytype_to_use))
    
    result <- tryCatch({
      enrichGO(
        gene = gene_data$genes,
        OrgDb = orgdb,
        keyType = keytype_to_use,
        ont = "BP",
        pvalueCutoff = 1.0,  # Get all pathways, will filter for display later
        pAdjustMethod = "BH",
        universe = bg_data$genes,
        readable = TRUE
      )
    }, error = function(e) {
      cat("  Error in GO enrichment:", e$message, "\n")
      return(NULL)
    })
    
  } else if (params$enrichment_db == "gmt") {
    
    vcat("  Running GMT file enrichment...\n")
    
    # Load GMT file
    gmt_data <- load_gmt_file(params$gmt_file)
    
    # Perform hypergeometric test for each pathway
    results_list <- list()
    
    n_background <- length(background_genes)
    n_query <- length(gene_list)
    
    for (pathway_name in gmt_data$names) {
      pathway_genes <- gmt_data$pathways[[pathway_name]]
      pathway_genes <- pathway_genes[pathway_genes %in% background_genes]
      
      if (length(pathway_genes) == 0) next
      
      # Genes in both query and pathway
      overlap <- intersect(gene_list, pathway_genes)
      n_overlap <- length(overlap)
      
      if (n_overlap == 0) next
      
      # Hypergeometric test
      # P(X >= k) where X ~ Hypergeometric(N, K, n)
      # N = total background genes
      # K = genes in pathway
      # n = query genes
      # k = overlap
      
      p_value <- phyper(
        q = n_overlap - 1,
        m = length(pathway_genes),
        n = n_background - length(pathway_genes),
        k = n_query,
        lower.tail = FALSE
      )
      
      results_list[[pathway_name]] <- data.frame(
        ID = pathway_name,
        Description = gmt_data$descriptions[which(gmt_data$names == pathway_name)],
        GeneRatio = paste0(n_overlap, "/", n_query),
        BgRatio = paste0(length(pathway_genes), "/", n_background),
        pvalue = p_value,
        Count = n_overlap,
        geneID = paste(overlap, collapse = "/"),
        stringsAsFactors = FALSE
      )
    }
    
    if (length(results_list) == 0) {
      cat("  No enriched pathways found\n")
      return(NULL)
    }
    
    # Combine results
    result_df <- do.call(rbind, results_list)
    rownames(result_df) <- NULL
    
    # Adjust p-values
    result_df$p.adjust <- p.adjust(result_df$pvalue, method = "BH")
    
    # Sort by p-value (keep all pathways, not just significant)
    result_df <- result_df[order(result_df$pvalue), ]
    
    # Convert to enrichResult-like object for compatibility
    result <- result_df
    
  }
  
  return(result)
}

# Perform enrichment for all clusters
perform_cluster_enrichment <- function(cluster_result, all_genes_results, params, output_dir) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("PATHWAY ENRICHMENT ANALYSIS\n")
  cat(rep("=", 70), "\n", sep = "")
  
  cat(sprintf("Database: %s\n", params$enrichment_db))
  if (params$enrichment_db != "gmt") {
    cat(sprintf("Organism: %s\n", params$organism))
  }
  cat(sprintf("P-value cutoff: %.3f\n", params$enrichment_pval))
  
  # Create enrichment directory
  enrich_dir <- file.path(output_dir, "enrichment")
  if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)
  
  # Get background genes (all genes tested)
  # Use top 500 per cluster OR all FDR < 0.05, whichever is less
  pval_col <- if ("adj.P.Val" %in% colnames(all_genes_results)) "adj.P.Val" else "padj"
  
  sig_genes <- rownames(all_genes_results)[all_genes_results[[pval_col]] < 0.05 & 
                                           !is.na(all_genes_results[[pval_col]])]
  
  background_genes <- rownames(all_genes_results)
  vcat(sprintf("Background: %d genes\n", length(background_genes)))
  vcat(sprintf("Significant genes: %d genes (FDR < 0.05)\n", length(sig_genes)))
  
  # Store all enrichment results
  all_enrichments <- list()
  enrichment_summaries <- list()
  
  # Perform enrichment for each cluster
  for (i in 1:length(unique(cluster_result$clusters))) {
    cluster_genes_all <- names(cluster_result$clusters)[cluster_result$clusters == i]
    
    # Get genes for this cluster (top 500 or all sig genes, whichever is less)
    cluster_genes_sig <- cluster_genes_all[cluster_genes_all %in% sig_genes]
    
    if (length(cluster_genes_sig) > 500) {
      # Take top 500 by significance
      cluster_results <- all_genes_results[cluster_genes_sig, , drop = FALSE]
      cluster_results <- cluster_results[order(cluster_results[[pval_col]]), ]
      cluster_genes <- rownames(cluster_results)[1:500]
    } else {
      cluster_genes <- cluster_genes_sig
    }
    
    cat(sprintf("\nCluster %d (%s): %d genes for enrichment\n", 
                i, cluster_result$pattern_types[i], length(cluster_genes)))
    
    if (length(cluster_genes) < 5) {
      cat("  Skipping: too few genes\n")
      next
    }
    
    # Run enrichment
    enrich_result <- run_enrichment_analysis(
      cluster_genes,
      background_genes,
      params
    )
    
    if (is.null(enrich_result)) {
      cat("  No enrichment results\n")
      next
    }
    
    # Convert to data frame if needed
    if (is(enrich_result, "enrichResult")) {
      enrich_df <- as.data.frame(enrich_result)
    } else {
      enrich_df <- enrich_result
    }
    
    if (nrow(enrich_df) == 0) {
      cat("  No pathways tested\n")
      next
    }
    
    # Add significance column based on user threshold
    # Check which adjusted p-value column exists
    if ("p.adjust" %in% colnames(enrich_df)) {
      enrich_df$Significant <- ifelse(enrich_df$p.adjust < params$enrichment_pval, "Yes", "No")
      n_sig <- sum(enrich_df$p.adjust < params$enrichment_pval, na.rm = TRUE)
    } else if ("padj" %in% colnames(enrich_df)) {
      enrich_df$Significant <- ifelse(enrich_df$padj < params$enrichment_pval, "Yes", "No")
      n_sig <- sum(enrich_df$padj < params$enrichment_pval, na.rm = TRUE)
    } else if ("qvalue" %in% colnames(enrich_df)) {
      enrich_df$Significant <- ifelse(enrich_df$qvalue < params$enrichment_pval, "Yes", "No")
      n_sig <- sum(enrich_df$qvalue < params$enrichment_pval, na.rm = TRUE)
    } else {
      enrich_df$Significant <- "Unknown"
      n_sig <- 0
    }
    
    cat(sprintf("  Tested %d pathways, %d significant (FDR < %.2f)\n", 
                nrow(enrich_df), n_sig, params$enrichment_pval))
    
    # Save results (all pathways, with significance column)
    outfile <- file.path(enrich_dir, 
                        sprintf("cluster_%d_%s_enrichment.txt", i,
                               gsub(" ", "_", cluster_result$pattern_types[i])))
    
    write.table(enrich_df, file = outfile, sep = "\t", quote = FALSE,
                row.names = FALSE)
    
    # Store for summary
    all_enrichments[[as.character(i)]] <- enrich_result
    
    # Create summary entry
    if (nrow(enrich_df) > 0) {
      enrichment_summaries[[as.character(i)]] <- data.frame(
        Cluster = i,
        Pattern = cluster_result$pattern_types[i],
        N_genes = length(cluster_genes),
        N_pathways_tested = nrow(enrich_df),
        N_pathways_sig = n_sig,
        Top_pathway = enrich_df$Description[1],
        Top_pvalue = enrich_df$pvalue[1]
      )
    }
  }
  
  # Create summary table
  if (length(enrichment_summaries) > 0) {
    summary_df <- do.call(rbind, enrichment_summaries)
    write.table(summary_df, 
                file = file.path(enrich_dir, "enrichment_summary.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("ENRICHMENT SUMMARY\n")
    cat(rep("=", 70), "\n", sep = "")
    print(summary_df, row.names = FALSE)
    cat("\n")
  }
  
  # Create visualizations if we have results
  if (length(all_enrichments) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_enrichment_summary(all_enrichments, cluster_result, enrich_dir, params)
  }
  
  cat(sprintf("\nEnrichment results saved to: %s\n", enrich_dir))
  
  return(all_enrichments)
}

# Plot enrichment summary
plot_enrichment_summary <- function(all_enrichments, cluster_result, output_dir, params) {
  vcat("  Creating enrichment visualization...\n")
  
  # Combine top pathways from each cluster
  plot_data_list <- list()
  
  for (cluster_id in names(all_enrichments)) {
    enrich_result <- all_enrichments[[cluster_id]]
    
    # Convert to data frame
    if (is(enrich_result, "enrichResult")) {
      enrich_df <- as.data.frame(enrich_result)
    } else {
      enrich_df <- enrich_result
    }
    
    if (nrow(enrich_df) == 0) next
    
    # Filter to significant pathways only for visualization
    if ("p.adjust" %in% colnames(enrich_df)) {
      sig_df <- enrich_df[enrich_df$p.adjust < params$enrichment_pval, ]
    } else if ("padj" %in% colnames(enrich_df)) {
      sig_df <- enrich_df[enrich_df$padj < params$enrichment_pval, ]
    } else if ("qvalue" %in% colnames(enrich_df)) {
      sig_df <- enrich_df[enrich_df$qvalue < params$enrichment_pval, ]
    } else {
      sig_df <- enrich_df  # If no adjusted p-value column, use all
    }
    
    if (nrow(sig_df) == 0) next
    
    # Take top 10 significant pathways
    top_pathways <- head(sig_df, 10)
    top_pathways$Cluster <- as.integer(cluster_id)
    top_pathways$Pattern <- cluster_result$pattern_types[as.integer(cluster_id)]
    
    plot_data_list[[cluster_id]] <- top_pathways
  }
  
  if (length(plot_data_list) == 0) return()
  
  plot_data <- do.call(rbind, plot_data_list)
  
  # Create cluster label
  plot_data$Cluster_label <- paste0("Cluster ", plot_data$Cluster, ": ", plot_data$Pattern)
  
  # Shorten pathway names if too long
  plot_data$Description_short <- ifelse(
    nchar(as.character(plot_data$Description)) > 50,
    paste0(substr(plot_data$Description, 1, 47), "..."),
    as.character(plot_data$Description)
  )
  
  # Create dotplot
  p <- ggplot(plot_data, aes(x = Cluster_label, y = Description_short)) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Top Enriched Pathways per Cluster",
         x = "Cluster", y = "Pathway",
         size = "Gene Count", color = "P-value")
  
  ggsave(file.path(output_dir, "enrichment_dotplot.pdf"), p, 
         width = 12, height = max(8, nrow(plot_data) * 0.3))
  
  # Create barplot for top pathway per cluster
  top_per_cluster <- do.call(rbind, lapply(split(plot_data, plot_data$Cluster), function(x) x[1, ]))
  
  p2 <- ggplot(top_per_cluster, aes(x = reorder(Cluster_label, -pvalue), 
                                     y = -log10(pvalue))) +
    geom_bar(stat = "identity", aes(fill = Pattern)) +
    coord_flip() +
    theme_bw() +
    geom_text(aes(label = Description_short), hjust = -0.1, size = 3) +
    labs(title = "Top Pathway per Cluster",
         x = "Cluster", y = "-log10(P-value)") +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "enrichment_top_per_cluster.pdf"), p2,
         width = 12, height = max(6, length(unique(plot_data$Cluster)) * 0.8))
}

################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

plot_timecourse_profiles <- function(counts, design, time_column, genes, 
                                     output_file, method_name) {
  vcat("  Creating timecourse profile plot...\n")
  
  # Prepare data for plotting
  time_var <- as.numeric(as.character(design[[time_column]]))
  
  plot_data <- data.frame(
    Gene = rep(genes, each = ncol(counts)),
    Sample = rep(colnames(counts), length(genes)),
    Time = rep(time_var, length(genes)),
    Expression = as.vector(t(counts[genes, ]))
  )
  
  # Plot
  p <- ggplot(plot_data, aes(x = Time, y = Expression, group = Gene)) +
    geom_line(alpha = 0.5) +
    geom_smooth(aes(group = 1), method = "loess", color = "red", se = TRUE) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    theme_bw() +
    labs(title = paste0("Top Genes - ", method_name),
         x = "Time", y = "Normalized Expression")
  
  ggsave(output_file, p, width = 16, height = 10)
  
  return(p)
}

# New function for comparing timecourse profiles between two groups
plot_group_comparison_profiles <- function(counts, design, time_column, group_column,
                                          compare_groups, genes, output_file, method_name) {
  vcat("  Creating group comparison profile plot...\n")
  
  # Prepare data for plotting
  time_var <- as.numeric(as.character(design[[time_column]]))
  group_var <- factor(design[[group_column]], levels = compare_groups)
  
  plot_data <- data.frame(
    Gene = rep(genes, each = ncol(counts)),
    Sample = rep(colnames(counts), length(genes)),
    Time = rep(time_var, length(genes)),
    Group = rep(group_var, length(genes)),
    Expression = as.vector(t(counts[genes, ]))
  )
  
  # Plot with groups colored differently
  p <- ggplot(plot_data, aes(x = Time, y = Expression, color = Group, group = Gene)) +
    geom_line(alpha = 0.3, aes(group = interaction(Gene, Group))) +
    geom_smooth(aes(group = Group), method = "loess", se = TRUE, size = 1.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("#E41A1C", "#377EB8")) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(title = paste0("Differential Timecourse - ", method_name),
         subtitle = paste0(compare_groups[1], " vs ", compare_groups[2]),
         x = "Time", y = "Normalized Expression")
  
  ggsave(output_file, p, width = 16, height = 12)
  
  return(p)
}

plot_heatmap <- function(counts, design, time_column, genes, output_file, 
                        method_name, group_column = NULL) {
  vcat("  Creating heatmap...\n")
  
  # Scale expression values
  mat <- counts[genes, ]
  mat <- t(scale(t(mat)))
  
  # Get time values
  time_values <- design[[time_column]]
  
  # Determine sample ordering
  if (!is.null(group_column) && group_column %in% colnames(design)) {
    # INTERACTION ANALYSIS: Sort by group first, then by time within each group
    vcat("  Sorting samples by group, then time within group\n")
    
    group_values <- design[[group_column]]
    
    # Convert time to numeric for sorting if needed
    if (is.numeric(time_values)) {
      time_numeric <- time_values
    } else {
      time_numeric <- suppressWarnings(as.numeric(as.character(time_values)))
      if (all(is.na(time_numeric))) {
        time_numeric <- as.numeric(factor(time_values))
      }
    }
    
    # Sort by group first, then by time within each group
    sample_order <- order(group_values, time_numeric)
    
  } else {
    # STANDARD ANALYSIS: Sort by time only
    if (is.numeric(time_values)) {
      sample_order <- order(time_values)
    } else {
      # For categorical, try to convert to numeric for sorting, otherwise keep as is
      time_numeric <- suppressWarnings(as.numeric(as.character(time_values)))
      if (all(!is.na(time_numeric))) {
        sample_order <- order(time_numeric)
      } else {
        # If can't convert to numeric, keep original order
        sample_order <- 1:length(time_values)
      }
    }
  }
  
  # Reorder matrix and design by the determined order
  mat <- mat[, sample_order, drop = FALSE]
  time_values <- time_values[sample_order]
  
  # Create annotation
  annotation_col <- data.frame(
    Time = factor(time_values, levels = sort(unique(time_values))),
    row.names = colnames(mat)
  )
  
  # Add group annotation if this is an interaction analysis
  if (!is.null(group_column) && group_column %in% colnames(design)) {
    group_values <- design[[group_column]][sample_order]
    annotation_col$Group <- factor(group_values)
  }
  
  # Create DISTINCT color palette for each unique timepoint
  unique_times <- sort(unique(time_values))
  n_times <- length(unique_times)
  
  # Choose color palette based on number of timepoints
  if (n_times <= 9) {
    # Use Set1 for up to 9 timepoints (bright, distinct colors)
    time_colors <- brewer.pal(max(3, n_times), "Set1")[1:n_times]
  } else if (n_times <= 12) {
    # Use Set3 for 10-12 timepoints (pastel, distinct colors)
    time_colors <- brewer.pal(n_times, "Set3")
  } else {
    # For many timepoints, create distinct colors
    time_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_times)
  }
  
  names(time_colors) <- as.character(unique_times)
  
  # Create annotation colors
  annotation_colors <- list(Time = time_colors)
  
  # Add group colors if interaction analysis
  if (!is.null(group_column) && group_column %in% colnames(design)) {
    group_values <- design[[group_column]][sample_order]
    unique_groups <- unique(group_values)
    n_groups <- length(unique_groups)
    
    # Use contrasting colors for groups
    if (n_groups == 2) {
      group_colors <- c("#E41A1C", "#377EB8")  # Red and Blue
    } else {
      group_colors <- brewer.pal(max(3, n_groups), "Dark2")[1:n_groups]
    }
    names(group_colors) <- as.character(unique_groups)
    annotation_colors$Group <- group_colors
  }
  
  # Colors for heatmap
  heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # Create title
  if (!is.null(group_column) && group_column %in% colnames(design)) {
    title <- paste0("Top ", length(genes), " Genes (grouped by treatment, ordered by time) - ", method_name)
  } else {
    title <- paste0("Top ", length(genes), " Genes (ordered by time) - ", method_name)
  }
  
  # Plot with legend, columns ordered appropriately
  pdf(output_file, width = 12, height = 12)
  pheatmap(mat,
           color = heatmap_colors,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           annotation_legend = TRUE,
           cluster_cols = FALSE,  # Don't cluster - already ordered
           cluster_rows = TRUE,   # Cluster genes by expression pattern
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = title,
           fontsize_row = max(6, 15 - length(genes)/10),
           fontsize_col = 8,
           breaks = seq(-2, 2, length.out = 101),
           legend_breaks = c(-2, -1, 0, 1, 2),
           legend_labels = c("-2", "-1", "0", "1", "2"))
  dev.off()
  
  vcat("  Samples ordered: ", paste(colnames(mat), collapse=", "), "\n")
  if (!is.null(group_column) && group_column %in% colnames(design)) {
    vcat("  Grouped by: ", group_column, "\n")
  }
  vcat("  Timepoint colors: ", paste(names(time_colors), "=", time_colors, collapse=", "), "\n")
}

################################################################################
# COMPARISON FUNCTIONS
################################################################################

compare_methods <- function(limma_results, deseq2_results, params) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("COMPARING LIMMA AND DESEQ2 RESULTS\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Get common genes
  common_genes <- intersect(rownames(limma_results), rownames(deseq2_results))
  
  # Extract p-values
  limma_pval <- limma_results[common_genes, "adj.P.Val"]
  deseq2_pval <- deseq2_results[common_genes, "padj"]
  
  # Remove NAs
  valid <- !is.na(limma_pval) & !is.na(deseq2_pval)
  limma_pval <- limma_pval[valid]
  deseq2_pval <- deseq2_pval[valid]
  genes <- common_genes[valid]
  
  # Correlation
  cor_pval <- cor(-log10(limma_pval), -log10(deseq2_pval), 
                  use = "complete.obs", method = "spearman")
  
  cat(sprintf("\nCorrelation of -log10(p-values): %.3f\n", cor_pval))
  
  # Overlap of significant genes
  limma_sig <- names(limma_pval)[limma_pval < params$fdr_threshold]
  deseq2_sig <- names(deseq2_pval)[deseq2_pval < params$fdr_threshold]
  
  overlap <- length(intersect(limma_sig, deseq2_sig))
  
  cat(sprintf("Significant genes (FDR < %.2f):\n", params$fdr_threshold))
  cat(sprintf("  limma: %d\n", length(limma_sig)))
  cat(sprintf("  DESeq2: %d\n", length(deseq2_sig)))
  cat(sprintf("  Overlap: %d (%.1f%% of union)\n", 
              overlap, 100 * overlap / length(union(limma_sig, deseq2_sig))))
  
  # Scatter plot
  if (params$generate_plots) {
    pdf(file.path(params$output_dir, "figures", "method_comparison_pvalues.pdf"), 
        width = 8, height = 8)
    plot(-log10(limma_pval), -log10(deseq2_pval),
         xlab = "-log10(limma FDR)", ylab = "-log10(DESeq2 FDR)",
         main = "Comparison of limma and DESeq2",
         pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.5)
    abline(0, 1, col = "red", lty = 2)
    abline(h = -log10(params$fdr_threshold), col = "blue", lty = 2)
    abline(v = -log10(params$fdr_threshold), col = "blue", lty = 2)
    legend("topleft", 
           legend = c(sprintf("Correlation: %.3f", cor_pval),
                     sprintf("Overlap: %d genes", overlap)),
           bty = "n")
    dev.off()
  }
  
  # Comparison data
  comparison_data <- data.frame(
    Gene = genes,
    limma_FDR = limma_pval,
    DESeq2_FDR = deseq2_pval,
    limma_sig = limma_pval < params$fdr_threshold,
    DESeq2_sig = deseq2_pval < params$fdr_threshold
  )
  
  write.table(comparison_data,
              file = file.path(params$output_dir, "comparison", "method_comparison.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(comparison_data)
}

################################################################################
# MAIN ANALYSIS
################################################################################

main <- function() {
  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("TIMECOURSE RNA-SEQ ANALYSIS\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Print configuration
  cat("Configuration:\n")
  cat(sprintf("  Counts file: %s\n", params$counts_file))
  cat(sprintf("  Design file: %s\n", params$design_file))
  cat(sprintf("  Time column: %s (%s)\n", 
              params$time_column, 
              if(params$time_as_numeric) "numeric" else "categorical"))
  if (params$time_as_numeric && !is.null(params$spline_df)) {
    cat(sprintf("  Spline df: %d\n", params$spline_df))
  }
  if (length(params$covariates) > 0) {
    cat(sprintf("  Covariates: %s\n", paste(params$covariates, collapse=", ")))
  }
  cat(sprintf("  Method(s): %s\n", opt$method))
  cat(sprintf("  Output directory: %s\n", params$output_dir))
  cat("\n")
  
  # Create output directories
  create_output_dirs(params$output_dir)
  
  # Load data
  data <- load_data(params$counts_file, params$design_file, params$duplicate_genes)
  counts <- data$counts
  design <- data$design
  
  # Filter low-expressed genes
  counts <- filter_genes(counts, params$min_count, params$min_samples)
  
  # Check time column exists
  if (!params$time_column %in% colnames(design)) {
    stop(sprintf("Time column '%s' not found in design file!\nAvailable columns: %s", 
                params$time_column, paste(colnames(design), collapse=", ")))
  }
  
  # Initialize results storage
  limma_res <- NULL
  deseq2_res <- NULL
  
  # Run limma analysis
  if (run_limma) {
    limma_res <- run_limma_analysis(counts, design, params)
    
    # Get top genes
    top_genes_limma <- rownames(limma_res$results)[1:min(params$top_n_genes, 
                                                          nrow(limma_res$results))]
    
    # Visualizations
    if (params$generate_plots && params$time_as_numeric) {
      if (params$test_interaction) {
        # Group comparison plots
        plot_group_comparison_profiles(
          limma_res$voom$E, limma_res$design, params$time_column, params$group_column,
          params$compare_groups,
          top_genes_limma[1:min(20, length(top_genes_limma))],
          file.path(params$output_dir, "figures", "limma_interaction_profiles.pdf"),
          "limma-voom"
        )
        
        plot_heatmap(
          limma_res$voom$E, limma_res$design, params$time_column, top_genes_limma,
          file.path(params$output_dir, "figures", "limma_interaction_heatmap.pdf"),
          "limma-voom (Interaction)",
          group_column = params$group_column
        )
      } else {
        # Standard time effect plots
        plot_timecourse_profiles(
          limma_res$voom$E, design, params$time_column, 
          top_genes_limma[1:min(20, length(top_genes_limma))],
          file.path(params$output_dir, "figures", "limma_top_genes_profiles.pdf"),
          "limma-voom"
        )
        
        plot_heatmap(
          limma_res$voom$E, design, params$time_column, top_genes_limma,
          file.path(params$output_dir, "figures", "limma_top_genes_heatmap.pdf"),
          "limma-voom"
        )
      }
    }
  }
  
  # Run DESeq2 analysis
  if (run_deseq2) {
    deseq2_res <- run_deseq2_analysis(counts, design, params)
    
    # Get top genes
    top_genes_deseq2 <- rownames(deseq2_res$results)[1:min(params$top_n_genes, 
                                                            nrow(deseq2_res$results))]
    
    # Visualizations
    if (params$generate_plots && params$time_as_numeric) {
      if (params$test_interaction) {
        # Group comparison plots
        plot_group_comparison_profiles(
          log2(deseq2_res$norm_counts + 1), deseq2_res$design, params$time_column,
          params$group_column, params$compare_groups,
          top_genes_deseq2[1:min(20, length(top_genes_deseq2))],
          file.path(params$output_dir, "figures", "deseq2_interaction_profiles.pdf"),
          "DESeq2"
        )
        
        plot_heatmap(
          log2(deseq2_res$norm_counts + 1), deseq2_res$design, params$time_column, 
          top_genes_deseq2,
          file.path(params$output_dir, "figures", "deseq2_interaction_heatmap.pdf"),
          "DESeq2 (Interaction)",
          group_column = params$group_column
        )
      } else {
        # Standard time effect plots
        plot_timecourse_profiles(
          log2(deseq2_res$norm_counts + 1), design, params$time_column,
          top_genes_deseq2[1:min(20, length(top_genes_deseq2))],
          file.path(params$output_dir, "figures", "deseq2_top_genes_profiles.pdf"),
          "DESeq2"
        )
        
        plot_heatmap(
          log2(deseq2_res$norm_counts + 1), design, params$time_column, top_genes_deseq2,
          file.path(params$output_dir, "figures", "deseq2_top_genes_heatmap.pdf"),
          "DESeq2"
        )
      }
    }
  }
  
  # Perform clustering if requested
  if (params$cluster_genes && params$time_as_numeric && params$generate_plots) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("TEMPORAL PATTERN CLUSTERING\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Use limma results if available, otherwise DESeq2
    if (run_limma) {
      cluster_method_name <- "limma-voom"
      sig_genes <- rownames(limma_res$results)[limma_res$results$adj.P.Val < params$fdr_threshold]
      expr_data <- limma_res$voom$E
      results_for_clustering <- limma_res$results
      analysis_design <- if (params$test_interaction) limma_res$design else design
    } else if (run_deseq2) {
      cluster_method_name <- "DESeq2"
      sig_genes <- rownames(deseq2_res$results)[deseq2_res$results$padj < params$fdr_threshold & 
                                                 !is.na(deseq2_res$results$padj)]
      expr_data <- log2(deseq2_res$norm_counts + 1)
      results_for_clustering <- deseq2_res$results
      analysis_design <- if (params$test_interaction) deseq2_res$design else design
    }
    
    if (length(sig_genes) < params$n_clusters) {
      cat(sprintf("  Warning: Only %d significant genes found, need at least %d for clustering\n",
                 length(sig_genes), params$n_clusters))
      cat("  Skipping clustering analysis\n")
    } else if (length(sig_genes) < 10) {
      cat(sprintf("  Warning: Only %d significant genes found, clustering may not be meaningful\n",
                 length(sig_genes)))
      cat("  Skipping clustering analysis\n")
    } else {
      cat(sprintf("Clustering %d significant genes (FDR < %.2f)\n", 
                 length(sig_genes), params$fdr_threshold))
      
      # Subset to significant genes
      expr_sig <- expr_data[sig_genes, , drop = FALSE]
      
      # Perform clustering
      cluster_result <- cluster_temporal_patterns(
        expr_sig, 
        analysis_design, 
        params$time_column,
        params$n_clusters,
        params$cluster_method
      )
      
      # Create visualizations
      plot_cluster_profiles(
        cluster_result,
        file.path(params$output_dir, "figures", 
                 paste0(tolower(gsub("-", "_", cluster_method_name)), "_cluster_profiles.pdf")),
        cluster_method_name
      )
      
      plot_cluster_heatmap(
        cluster_result,
        file.path(params$output_dir, "figures",
                 paste0(tolower(gsub("-", "_", cluster_method_name)), "_cluster_heatmap.pdf")),
        cluster_method_name
      )
      
      # Create summary
      create_cluster_summary(
        cluster_result,
        file.path(params$output_dir, 
                 tolower(gsub("-", "_", cluster_method_name)),
                 "cluster_summary.txt")
      )
      
      # Export gene lists
      export_cluster_gene_lists(
        cluster_result,
        results_for_clustering[sig_genes, , drop = FALSE],
        file.path(params$output_dir, tolower(gsub("-", "_", cluster_method_name)))
      )
      
      # Create per-cluster top genes plots
      cat("\n", rep("=", 70), "\n", sep = "")
      cat("PER-CLUSTER TOP GENES VISUALIZATIONS\n")
      cat(rep("=", 70), "\n", sep = "")
      
      plot_cluster_top_genes_profiles(
        cluster_result,
        expr_data,
        results_for_clustering,
        analysis_design,
        params$time_column,
        file.path(params$output_dir, tolower(gsub("-", "_", cluster_method_name))),
        cluster_method_name,
        n_top = 20
      )
      
      plot_cluster_top_genes_heatmaps(
        cluster_result,
        expr_data,
        results_for_clustering,
        analysis_design,
        params$time_column,
        file.path(params$output_dir, tolower(gsub("-", "_", cluster_method_name))),
        cluster_method_name,
        n_top = 50,
        group_column = if (params$test_interaction) params$group_column else NULL
      )
      
      # Run enrichment analysis if requested
      if (params$run_enrichment) {
        perform_cluster_enrichment(
          cluster_result,
          results_for_clustering,
          params,
          file.path(params$output_dir, tolower(gsub("-", "_", cluster_method_name)))
        )
      }
    }
  }
  
  # Compare methods if both were run
  if (run_limma && run_deseq2) {
    compare_methods(limma_res$results, deseq2_res$results, params)
  }
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS COMPLETE!\n")
  cat(rep("=", 70), "\n", sep = "")
  cat(sprintf("\nResults saved to: %s\n\n", params$output_dir))
}

# Run the analysis
main()
