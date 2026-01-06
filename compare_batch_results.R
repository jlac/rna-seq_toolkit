#!/usr/bin/env Rscript

# Cross-Group Comparison Script for Batch Timecourse Analysis
# Compares genes and pathways across multiple group analyses

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input-dir"), type="character", default=NULL,
              help="Input directory containing batch results", metavar="DIR"),
  make_option(c("-o", "--output-dir"), type="character", default="batch_comparison",
              help="Output directory for comparison results [default: %default]", metavar="DIR"),
  make_option(c("--method"), type="character", default="deseq2",
              help="Analysis method used: limma or deseq2 [default: %default]", metavar="STRING"),
  make_option(c("--fdr-cutoff"), type="double", default=0.05,
              help="FDR cutoff for significant genes [default: %default]", metavar="FLOAT"),
  make_option(c("--top-n"), type="integer", default=20,
              help="Number of top pathways to show [default: %default]", metavar="INT"),
  make_option(c("--annotation"), type="character", default=NULL,
              help="Gene annotation file with Gene_ID and Gene_Symbol columns", metavar="FILE"),
  make_option(c("--organism"), type="character", default="human",
              help="Organism: human or mouse [default: %default]", metavar="STRING"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$`input-dir`)) {
  stop("Error: --input-dir is required\n")
}

params <- list(
  input_dir = opt$`input-dir`,
  output_dir = opt$`output-dir`,
  method = opt$method,
  fdr_cutoff = opt$`fdr-cutoff`,
  top_n = opt$`top-n`,
  annotation = opt$annotation,
  organism = opt$organism,
  verbose = opt$verbose
)

vcat <- function(...) {
  if (params$verbose) cat(...)
}

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("CROSS-GROUP COMPARISON ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Configuration:\n")
cat(sprintf("  Input directory: %s\n", params$input_dir))
cat(sprintf("  Output directory: %s\n", params$output_dir))
cat(sprintf("  Method: %s\n", params$method))
cat(sprintf("  FDR cutoff: %.3f\n", params$fdr_cutoff))
cat(sprintf("  Organism: %s\n", params$organism))
if (!is.null(params$annotation)) {
  cat(sprintf("  Annotation file: %s\n", params$annotation))
}
cat("\n")

# Create output directories
dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(params$output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(params$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Find all group directories
group_dirs <- list.dirs(params$input_dir, recursive = FALSE, full.names = FALSE)
group_dirs <- group_dirs[!grepl("_design\\.txt$", group_dirs)]
group_dirs <- group_dirs[group_dirs != "cross_group_comparison"]  # Exclude comparison dir

if (length(group_dirs) == 0) {
  stop("No group directories found in ", params$input_dir)
}

cat(sprintf("Found %d groups to compare:\n", length(group_dirs)))
for (g in group_dirs) {
  cat(sprintf("  - %s\n", g))
}
cat("\n")

################################################################################
# LOAD RESULTS FROM ALL GROUPS
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("LOADING RESULTS\n")
cat(rep("=", 80), "\n\n", sep = "")

all_results <- list()
all_enrichments <- list()

for (group in group_dirs) {
  vcat(sprintf("Loading results for %s...\n", group))
  
  # Load DE results - try both time effect and interaction files
  results_file <- file.path(params$input_dir, group, params$method,
                           paste0(params$method, "_time_effect_results.txt"))
  
  # If time effect file doesn't exist, try interaction file
  if (!file.exists(results_file)) {
    results_file <- file.path(params$input_dir, group, params$method,
                             paste0(params$method, "_interaction_results.txt"))
  }
  
  if (file.exists(results_file)) {
    res <- read.table(results_file, header = TRUE, sep = "\t", 
                     stringsAsFactors = FALSE, check.names = FALSE)
    all_results[[group]] <- res
    vcat(sprintf("  Loaded %d genes\n", nrow(res)))
  } else {
    warning(sprintf("Results file not found for %s (tried both time_effect and interaction files)", group))
  }
  
  # Load cluster assignments (if available)
  cluster_file <- file.path(params$input_dir, group, params$method,
                            "cluster_assignments.txt")
  if (file.exists(cluster_file)) {
    tryCatch({
      clusters <- read.table(cluster_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)
      if (group %in% names(all_results) && nrow(all_results[[group]]) > 0) {
        all_results[[group]]$Cluster <- clusters$Cluster[match(rownames(all_results[[group]]), 
                                                                clusters$Gene)]
        all_results[[group]]$Pattern <- clusters$Pattern[match(rownames(all_results[[group]]), 
                                                                clusters$Gene)]
        vcat(sprintf("  Loaded cluster assignments\n"))
      }
    }, error = function(e) {
      warning(sprintf("Could not load cluster assignments for %s: %s", group, e$message))
    })
  }
  
  # Load enrichment results (if available)
  enrich_dir <- file.path(params$input_dir, group, params$method, "enrichment")
  if (dir.exists(enrich_dir)) {
    enrich_files <- list.files(enrich_dir, pattern = "cluster_.*_enrichment\\.txt$",
                               full.names = TRUE)
    
    if (length(enrich_files) > 0) {
      all_enrichments[[group]] <- list()
      
      for (efile in enrich_files) {
        cluster_name <- sub(".*cluster_(\\d+)_(.+)_enrichment\\.txt", "\\1_\\2", 
                           basename(efile))
        enrich <- read.table(efile, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)
        all_enrichments[[group]][[cluster_name]] <- enrich
        vcat(sprintf("  Loaded enrichment for cluster %s: %d pathways\n",
                    cluster_name, nrow(enrich)))
      }
    }
  }
}

if (length(all_results) == 0) {
  stop("No results could be loaded from any group")
}

cat(sprintf("\nSuccessfully loaded results from %d groups\n\n", length(all_results)))

################################################################################
# HELPER FUNCTIONS
################################################################################

# Load gene annotation from file or create mapping
gene_annotation <- NULL

if (!is.null(params$annotation) && file.exists(params$annotation)) {
  vcat(sprintf("Loading gene annotation from %s...\n", params$annotation))
  gene_annotation <- read.table(params$annotation, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
  
  # Check for required columns
  if (!("Gene_ID" %in% colnames(gene_annotation)) || 
      !("Gene_Symbol" %in% colnames(gene_annotation))) {
    warning("Annotation file must have Gene_ID and Gene_Symbol columns")
    gene_annotation <- NULL
  } else {
    vcat(sprintf("  Loaded %d gene annotations\n", nrow(gene_annotation)))
  }
}

# Function to convert Entrez IDs to gene symbols
convert_entrez_to_symbol <- function(entrez_ids, organism = "human") {
  symbols <- rep(NA, length(entrez_ids))
  names(symbols) <- entrez_ids
  
  # Try using annotation database
  if (organism == "human") {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      vcat("  Converting Entrez IDs to symbols using org.Hs.eg.db...\n")
      tryCatch({
        library(org.Hs.eg.db)
        mapped <- mapIds(org.Hs.eg.db, 
                        keys = as.character(entrez_ids),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
        symbols <- as.character(mapped)
        names(symbols) <- entrez_ids
        vcat(sprintf("    Converted %d/%d IDs\n", sum(!is.na(symbols)), length(symbols)))
      }, error = function(e) {
        vcat(sprintf("    Error: %s\n", e$message))
      })
    }
  } else if (organism == "mouse") {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      vcat("  Converting Entrez IDs to symbols using org.Mm.eg.db...\n")
      tryCatch({
        library(org.Mm.eg.db)
        mapped <- mapIds(org.Mm.eg.db, 
                        keys = as.character(entrez_ids),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
        symbols <- as.character(mapped)
        names(symbols) <- entrez_ids
        vcat(sprintf("    Converted %d/%d IDs\n", sum(!is.na(symbols)), length(symbols)))
      }, error = function(e) {
        vcat(sprintf("    Error: %s\n", e$message))
      })
    }
  }
  
  return(symbols)
}

# Function to convert Ensembl IDs to gene symbols
convert_ensembl_to_symbol <- function(ensembl_ids, organism = "human") {
  symbols <- rep(NA, length(ensembl_ids))
  names(symbols) <- ensembl_ids
  
  if (organism == "human") {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      vcat("  Converting Ensembl IDs to symbols using org.Hs.eg.db...\n")
      tryCatch({
        library(org.Hs.eg.db)
        mapped <- mapIds(org.Hs.eg.db, 
                        keys = as.character(ensembl_ids),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
        symbols <- as.character(mapped)
        names(symbols) <- ensembl_ids
        vcat(sprintf("    Converted %d/%d IDs\n", sum(!is.na(symbols)), length(symbols)))
      }, error = function(e) {
        vcat(sprintf("    Error: %s\n", e$message))
      })
    }
  } else if (organism == "mouse") {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      vcat("  Converting Ensembl IDs to symbols using org.Mm.eg.db...\n")
      tryCatch({
        library(org.Mm.eg.db)
        mapped <- mapIds(org.Mm.eg.db, 
                        keys = as.character(ensembl_ids),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
        symbols <- as.character(mapped)
        names(symbols) <- ensembl_ids
        vcat(sprintf("    Converted %d/%d IDs\n", sum(!is.na(symbols)), length(symbols)))
      }, error = function(e) {
        vcat(sprintf("    Error: %s\n", e$message))
      })
    }
  }
  
  return(symbols)
}

# Function to detect gene ID type
detect_id_type <- function(gene_ids) {
  sample_ids <- head(gene_ids, 100)
  
  # Check for Ensembl (starts with ENS)
  if (sum(grepl("^ENS", sample_ids)) > length(sample_ids) * 0.8) {
    return("ensembl")
  }
  
  # Check for numeric Entrez IDs
  if (sum(grepl("^[0-9]+$", sample_ids)) > length(sample_ids) * 0.8) {
    return("entrez")
  }
  
  # Check for gene symbols (letters)
  if (sum(grepl("^[A-Z][A-Z0-9-]+$", sample_ids, ignore.case = TRUE)) > length(sample_ids) * 0.5) {
    return("symbol")
  }
  
  return("unknown")
}

# Function to get gene symbols from results
get_gene_symbols <- function(gene_ids, results_list) {
  vcat(sprintf("Converting %d gene IDs to symbols...\n", length(gene_ids)))
  
  symbols <- rep(NA, length(gene_ids))
  names(symbols) <- gene_ids
  
  # Method 1: Use provided annotation file
  if (!is.null(gene_annotation)) {
    vcat("  Using provided annotation file...\n")
    matched <- match(gene_ids, gene_annotation$Gene_ID)
    valid <- !is.na(matched)
    symbols[valid] <- gene_annotation$Gene_Symbol[matched[valid]]
    vcat(sprintf("    Matched %d/%d IDs\n", sum(!is.na(symbols)), length(symbols)))
  }
  
  # Method 2: Check results files for symbol columns
  if (sum(is.na(symbols)) > 0) {
    vcat("  Checking results files for symbol columns...\n")
    remaining_ids <- gene_ids[is.na(symbols)]
    
    for (group in names(results_list)) {
      res <- results_list[[group]]
      
      # Common column names for gene symbols
      symbol_cols <- c("gene_name", "symbol", "Symbol", "Gene", "SYMBOL", 
                      "gene_symbol", "GeneSymbol", "hgnc_symbol", "external_gene_name")
      
      for (col in symbol_cols) {
        if (col %in% colnames(res)) {
          matched <- match(remaining_ids, rownames(res))
          valid <- !is.na(matched)
          if (sum(valid) > 0) {
            symbols[gene_ids %in% remaining_ids][valid] <- res[matched[valid], col]
            vcat(sprintf("    Found %s column in %s: matched %d IDs\n", 
                        col, group, sum(valid)))
            remaining_ids <- gene_ids[is.na(symbols)]
            break
          }
        }
      }
      
      if (sum(is.na(symbols)) == 0) break
    }
  }
  
  # Method 3: Auto-detect ID type and convert
  if (sum(is.na(symbols)) > 0) {
    remaining_ids <- gene_ids[is.na(symbols)]
    id_type <- detect_id_type(remaining_ids)
    vcat(sprintf("  Detected ID type: %s\n", id_type))
    
    converted <- rep(NA, length(remaining_ids))
    
    if (id_type == "entrez") {
      converted <- convert_entrez_to_symbol(remaining_ids, params$organism)
    } else if (id_type == "ensembl") {
      converted <- convert_ensembl_to_symbol(remaining_ids, params$organism)
    } else if (id_type == "symbol") {
      vcat("  IDs appear to already be gene symbols\n")
      converted <- remaining_ids
    }
    
    if (!all(is.na(converted))) {
      symbols[gene_ids %in% remaining_ids] <- converted
    }
  }
  
  # Method 4: Fall back to gene IDs if still not resolved
  still_missing <- is.na(symbols)
  if (sum(still_missing) > 0) {
    vcat(sprintf("  Warning: Could not convert %d gene IDs, using IDs as symbols\n", 
                sum(still_missing)))
    symbols[still_missing] <- gene_ids[still_missing]
  }
  
  vcat(sprintf("  Final: %d gene symbols obtained\n", sum(!is.na(symbols))))
  
  return(symbols)
}

# Create gene annotation data frame
create_gene_annotation <- function(gene_ids, results_list) {
  symbols <- get_gene_symbols(gene_ids, results_list)
  
  data.frame(
    Gene_ID = gene_ids,
    Gene_Symbol = symbols,
    stringsAsFactors = FALSE
  )
}

################################################################################
# COMPARE SIGNIFICANT GENES
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("GENE OVERLAP ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Determine p-value column
pval_col <- if (params$method == "limma") "adj.P.Val" else "padj"

# Extract significant genes for each group
sig_genes_list <- lapply(all_results, function(res) {
  if (pval_col %in% colnames(res)) {
    rownames(res)[which(res[[pval_col]] < params$fdr_cutoff)]
  } else {
    character(0)
  }
})

# Count significant genes per group
sig_counts <- data.frame(
  Group = names(sig_genes_list),
  N_significant = sapply(sig_genes_list, length),
  stringsAsFactors = FALSE
)
sig_counts <- sig_counts[order(sig_counts$N_significant, decreasing = TRUE), ]

cat("Significant genes per group:\n")
print(sig_counts, row.names = FALSE)
cat("\n")

# Write to file
write.table(sig_counts,
           file = file.path(params$output_dir, "tables", "significant_genes_per_group.txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)

# Create Venn diagram data (for 2-5 groups)
if (length(sig_genes_list) >= 2 && length(sig_genes_list) <= 5) {
  vcat("Creating Venn diagram...\n")
  
  # Calculate overlaps
  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    library(VennDiagram)
    
    pdf(file.path(params$output_dir, "figures", "gene_overlap_venn.pdf"),
        width = 10, height = 10)
    
    if (length(sig_genes_list) == 2) {
      grid.draw(venn.diagram(
        x = sig_genes_list,
        category.names = names(sig_genes_list),
        filename = NULL,
        fill = brewer.pal(2, "Set2"),
        alpha = 0.5
      ))
    } else if (length(sig_genes_list) == 3) {
      grid.draw(venn.diagram(
        x = sig_genes_list,
        category.names = names(sig_genes_list),
        filename = NULL,
        fill = brewer.pal(3, "Set2"),
        alpha = 0.5
      ))
    } else if (length(sig_genes_list) == 4) {
      grid.draw(venn.diagram(
        x = sig_genes_list,
        category.names = names(sig_genes_list),
        filename = NULL,
        fill = brewer.pal(4, "Set2"),
        alpha = 0.5
      ))
    } else {
      grid.draw(venn.diagram(
        x = sig_genes_list,
        category.names = names(sig_genes_list),
        filename = NULL,
        fill = brewer.pal(5, "Set2"),
        alpha = 0.5
      ))
    }
    
    dev.off()
  }
}

# Create UpSet plot for complex overlaps (if available)
if (length(sig_genes_list) > 2 && requireNamespace("UpSetR", quietly = TRUE)) {
  vcat("Creating UpSet plot...\n")
  library(UpSetR)
  
  # Create presence/absence matrix
  all_genes <- unique(unlist(sig_genes_list))
  upset_data <- data.frame(
    Gene = all_genes,
    stringsAsFactors = FALSE
  )
  
  for (group in names(sig_genes_list)) {
    upset_data[[group]] <- as.integer(upset_data$Gene %in% sig_genes_list[[group]])
  }
  
  pdf(file.path(params$output_dir, "figures", "gene_overlap_upset.pdf"),
      width = 12, height = 8)
  print(upset(upset_data, sets = names(sig_genes_list), 
             order.by = "freq", keep.order = FALSE))
  dev.off()
}

# Find shared and unique genes
all_sig_genes <- unique(unlist(sig_genes_list))
shared_in_all <- Reduce(intersect, sig_genes_list)
unique_per_group <- lapply(names(sig_genes_list), function(g) {
  setdiff(sig_genes_list[[g]], unlist(sig_genes_list[names(sig_genes_list) != g]))
})
names(unique_per_group) <- names(sig_genes_list)

overlap_summary <- data.frame(
  Category = c("Total unique genes across all groups",
              "Genes significant in ALL groups",
              sapply(names(sig_genes_list), function(g) 
                paste0("Genes unique to ", g))),
  Count = c(length(all_sig_genes),
           length(shared_in_all),
           sapply(unique_per_group, length)),
  stringsAsFactors = FALSE
)

cat("\nGene overlap summary:\n")
print(overlap_summary, row.names = FALSE)
cat("\n")

write.table(overlap_summary,
           file = file.path(params$output_dir, "tables", "gene_overlap_summary.txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)

# Save shared genes with symbols
if (length(shared_in_all) > 0) {
  shared_genes_annot <- create_gene_annotation(shared_in_all, all_results)
  write.table(shared_genes_annot,
             file = file.path(params$output_dir, "tables", "genes_shared_in_all_groups.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved %d shared genes to genes_shared_in_all_groups.txt\n", 
             length(shared_in_all)))
}

# Save unique genes per group with symbols
for (group in names(unique_per_group)) {
  if (length(unique_per_group[[group]]) > 0) {
    unique_genes_annot <- create_gene_annotation(unique_per_group[[group]], all_results)
    write.table(unique_genes_annot,
               file = file.path(params$output_dir, "tables", 
                              paste0("genes_unique_to_", group, ".txt")),
               sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

################################################################################
# DIRECTION OF CHANGE COMPARISON
################################################################################

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("DIRECTION OF CHANGE ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Get log fold changes for shared genes
lfc_col <- if (params$method == "limma") "logFC" else "log2FoldChange"

if (length(shared_in_all) > 0) {
  lfc_matrix <- matrix(NA, nrow = length(shared_in_all), 
                      ncol = length(all_results))
  rownames(lfc_matrix) <- shared_in_all
  colnames(lfc_matrix) <- names(all_results)
  
  for (group in names(all_results)) {
    res <- all_results[[group]]
    if (lfc_col %in% colnames(res)) {
      matched_genes <- intersect(shared_in_all, rownames(res))
      lfc_matrix[matched_genes, group] <- res[matched_genes, lfc_col]
    }
  }
  
  # Create heatmap of log fold changes
  vcat("Creating log fold change heatmap for shared genes...\n")
  
  # Filter to genes with non-NA values
  complete_genes <- rownames(lfc_matrix)[complete.cases(lfc_matrix)]
  
  if (length(complete_genes) > 0) {
    lfc_subset <- lfc_matrix[complete_genes, , drop = FALSE]
    
    # Limit to top N genes by variance
    if (nrow(lfc_subset) > 100) {
      gene_vars <- apply(lfc_subset, 1, var, na.rm = TRUE)
      top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:100])
      lfc_subset <- lfc_subset[top_genes, , drop = FALSE]
    }
    
    # Replace rownames with gene symbols
    gene_symbols <- get_gene_symbols(rownames(lfc_subset), all_results)
    rownames(lfc_subset) <- gene_symbols
    
    # Adjust font size based on number of genes
    row_fontsize <- if (nrow(lfc_subset) <= 30) 8 else if (nrow(lfc_subset) <= 50) 6 else 5
    
    pdf(file.path(params$output_dir, "figures", "shared_genes_lfc_heatmap.pdf"),
        width = 10, height = max(10, nrow(lfc_subset) * 0.15))
    pheatmap(lfc_subset,
            color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
            breaks = seq(-3, 3, length.out = 101),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,  # Always show gene names
            main = "Log Fold Changes - Shared Genes",
            fontsize_row = row_fontsize)
    dev.off()
  }
  
  # Count concordant vs discordant changes
  direction_matrix <- sign(lfc_matrix)
  
  concordant <- apply(direction_matrix, 1, function(x) {
    all(x == x[1], na.rm = TRUE) && all(!is.na(x))
  })
  
  direction_summary <- data.frame(
    Category = c("Concordant (all same direction)",
                "Discordant (different directions)",
                "Up in all groups",
                "Down in all groups"),
    Count = c(sum(concordant, na.rm = TRUE),
             sum(!concordant, na.rm = TRUE),
             sum(apply(direction_matrix, 1, function(x) 
               all(x > 0, na.rm = TRUE) && all(!is.na(x))), na.rm = TRUE),
             sum(apply(direction_matrix, 1, function(x) 
               all(x < 0, na.rm = TRUE) && all(!is.na(x))), na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  cat("Direction of change for shared genes:\n")
  print(direction_summary, row.names = FALSE)
  cat("\n")
  
  write.table(direction_summary,
             file = file.path(params$output_dir, "tables", "direction_of_change_summary.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################################################
# CLUSTER PATTERN COMPARISON
################################################################################

# Check if any groups have cluster assignments
has_clusters <- sapply(all_results, function(res) "Pattern" %in% colnames(res))

if (any(has_clusters) && length(shared_in_all) > 0) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("TEMPORAL PATTERN COMPARISON\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  cat("Comparing temporal patterns across groups for shared genes...\n\n")
  
  # Create matrix of patterns for shared genes
  pattern_matrix <- matrix(NA, nrow = length(shared_in_all),
                          ncol = length(all_results))
  rownames(pattern_matrix) <- shared_in_all
  colnames(pattern_matrix) <- names(all_results)
  
  for (group in names(all_results)) {
    if ("Pattern" %in% colnames(all_results[[group]])) {
      res <- all_results[[group]]
      matched_genes <- intersect(shared_in_all, rownames(res))
      pattern_matrix[matched_genes, group] <- res[matched_genes, "Pattern"]
    }
  }
  
  # Get groups with cluster data
  groups_with_clusters <- names(all_results)[has_clusters]
  
  if (length(groups_with_clusters) >= 2) {
    # Filter to genes with pattern in at least 2 groups
    pattern_subset <- pattern_matrix[, groups_with_clusters, drop = FALSE]
    complete_patterns <- rowSums(!is.na(pattern_subset)) >= 2
    pattern_subset <- pattern_subset[complete_patterns, , drop = FALSE]
    
    if (nrow(pattern_subset) > 0) {
      cat(sprintf("Found %d shared genes with pattern assignments\n\n", nrow(pattern_subset)))
      
      # Count pattern concordance
      concordant_patterns <- apply(pattern_subset, 1, function(x) {
        x_clean <- x[!is.na(x)]
        if (length(x_clean) >= 2) {
          all(x_clean == x_clean[1])
        } else {
          FALSE
        }
      })
      
      # Get gene symbols for pattern matrix
      pattern_gene_symbols <- get_gene_symbols(rownames(pattern_subset), all_results)
      
      # Create pattern comparison table with gene symbols
      pattern_comparison <- data.frame(
        Gene_ID = rownames(pattern_subset),
        Gene_Symbol = pattern_gene_symbols,
        pattern_subset,
        Concordant = concordant_patterns,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      # Summary statistics
      pattern_summary <- data.frame(
        Category = c("Genes with same pattern in all groups",
                    "Genes with different patterns across groups"),
        Count = c(sum(concordant_patterns),
                 sum(!concordant_patterns)),
        stringsAsFactors = FALSE
      )
      
      cat("Pattern concordance:\n")
      print(pattern_summary, row.names = FALSE)
      cat("\n")
      
      # Count each pattern combination
      pattern_combos <- apply(pattern_subset, 1, function(x) {
        paste(x[!is.na(x)], collapse = " | ")
      })
      
      pattern_freq <- sort(table(pattern_combos), decreasing = TRUE)
      
      cat(sprintf("Top pattern combinations (showing up to 20):\n"))
      print(head(pattern_freq, 20))
      cat("\n")
      
      # Save pattern comparison table
      write.table(pattern_comparison,
                 file = file.path(params$output_dir, "tables", "cluster_pattern_comparison.txt"),
                 sep = "\t", quote = FALSE, row.names = FALSE)
      
      # Save concordant patterns
      concordant_genes <- pattern_comparison[pattern_comparison$Concordant, ]
      if (nrow(concordant_genes) > 0) {
        write.table(concordant_genes,
                   file = file.path(params$output_dir, "tables", "genes_with_concordant_patterns.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("Saved %d genes with concordant patterns\n", nrow(concordant_genes)))
      }
      
      # Save discordant patterns
      discordant_genes <- pattern_comparison[!pattern_comparison$Concordant, ]
      if (nrow(discordant_genes) > 0) {
        write.table(discordant_genes,
                   file = file.path(params$output_dir, "tables", "genes_with_discordant_patterns.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("Saved %d genes with discordant patterns\n", nrow(discordant_genes)))
      }
      
      # Create heatmap of patterns if there are distinct patterns
      unique_patterns <- unique(as.vector(pattern_subset))
      unique_patterns <- unique_patterns[!is.na(unique_patterns)]
      
      if (length(unique_patterns) >= 2) {
        vcat("\nCreating pattern comparison heatmap...\n")
        
        # Convert patterns to numeric for heatmap
        pattern_levels <- sort(unique(unique_patterns))
        pattern_numeric <- matrix(NA, nrow = nrow(pattern_subset), 
                                 ncol = ncol(pattern_subset))
        rownames(pattern_numeric) <- pattern_gene_symbols  # Use symbols
        colnames(pattern_numeric) <- colnames(pattern_subset)
        
        for (i in 1:nrow(pattern_subset)) {
          for (j in 1:ncol(pattern_subset)) {
            if (!is.na(pattern_subset[i, j])) {
              pattern_numeric[i, j] <- match(pattern_subset[i, j], pattern_levels)
            }
          }
        }
        
        # Limit to reasonable number for visualization
        if (nrow(pattern_numeric) > 100) {
          # Prioritize concordant genes
          concordant_indices <- which(concordant_patterns)
          discordant_indices <- which(!concordant_patterns)
          
          selected_indices <- c(
            head(concordant_indices, 50),
            head(discordant_indices, 50)
          )
          pattern_numeric <- pattern_numeric[selected_indices, , drop = FALSE]
        }
        
        # Create color palette
        n_patterns <- length(pattern_levels)
        pattern_colors <- colorRampPalette(brewer.pal(min(12, n_patterns), "Set3"))(n_patterns)
        
        # Adjust font size
        row_fontsize <- if (nrow(pattern_numeric) <= 30) 8 else if (nrow(pattern_numeric) <= 50) 6 else 5
        
        pdf(file.path(params$output_dir, "figures", "cluster_pattern_comparison_heatmap.pdf"),
            width = 10, height = max(10, nrow(pattern_numeric) * 0.15))
        
        tryCatch({
          pheatmap(pattern_numeric,
                  color = pattern_colors,
                  breaks = seq(0.5, n_patterns + 0.5, length.out = n_patterns + 1),
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  show_rownames = TRUE,
                  fontsize_row = row_fontsize,
                  main = "Temporal Pattern Comparison - Shared Genes",
                  legend_breaks = 1:n_patterns,
                  legend_labels = pattern_levels,
                  na_col = "grey90")
        }, error = function(e) {
          warning("Could not create pattern heatmap: ", e$message)
        })
        
        dev.off()
      }
    }
  }
}

################################################################################
# PATHWAY COMPARISON
################################################################################

if (length(all_enrichments) > 0) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("PATHWAY ENRICHMENT COMPARISON\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Combine all pathways across all groups and clusters
  all_pathways <- list()
  
  for (group in names(all_enrichments)) {
    for (cluster in names(all_enrichments[[group]])) {
      enrich <- all_enrichments[[group]][[cluster]]
      
      # Get significant pathways
      if ("p.adjust" %in% colnames(enrich)) {
        sig_pathways <- enrich$Description[enrich$p.adjust < params$fdr_cutoff]
      } else if ("padj" %in% colnames(enrich)) {
        sig_pathways <- enrich$Description[enrich$padj < params$fdr_cutoff]
      } else {
        sig_pathways <- enrich$Description[enrich$pvalue < params$fdr_cutoff]
      }
      
      # Filter NA and empty values
      sig_pathways <- sig_pathways[!is.na(sig_pathways) & nchar(sig_pathways) > 0]
      
      if (length(sig_pathways) > 0) {
        key <- paste(group, cluster, sep = "_")
        all_pathways[[key]] <- sig_pathways
      }
    }
  }
  
  if (length(all_pathways) > 0) {
    # Get unique pathways (filter NA values)
    all_pathway_values <- unlist(all_pathways)
    all_pathway_values <- all_pathway_values[!is.na(all_pathway_values)]
    unique_pathways <- unique(all_pathway_values)
    
    cat(sprintf("Found %d unique enriched pathways across all groups\n\n",
               length(unique_pathways)))
    
    # Create pathway presence/absence matrix
    pathway_matrix <- matrix(0, nrow = length(unique_pathways),
                            ncol = length(names(all_enrichments)))
    rownames(pathway_matrix) <- unique_pathways
    colnames(pathway_matrix) <- names(all_enrichments)
    
    for (group in names(all_enrichments)) {
      group_pathways <- unique(unlist(lapply(all_enrichments[[group]], function(x) {
        if ("p.adjust" %in% colnames(x)) {
          x$Description[x$p.adjust < params$fdr_cutoff]
        } else if ("padj" %in% colnames(x)) {
          x$Description[x$padj < params$fdr_cutoff]
        } else {
          x$Description[x$pvalue < params$fdr_cutoff]
        }
      })))
      
      # Only use pathways that exist in the matrix (handle NA, empty, and mismatches)
      group_pathways <- group_pathways[!is.na(group_pathways) & nchar(group_pathways) > 0]
      valid_pathways <- intersect(group_pathways, rownames(pathway_matrix))
      
      if (length(valid_pathways) > 0) {
        pathway_matrix[valid_pathways, group] <- 1
      }
    }
    
    # Filter to pathways present in at least 2 groups
    pathway_counts <- rowSums(pathway_matrix)
    shared_pathways <- names(pathway_counts[pathway_counts >= 2])
    
    if (length(shared_pathways) > 0) {
      vcat(sprintf("Creating pathway comparison heatmap for %d shared pathways...\n",
                  length(shared_pathways)))
      
      pathway_subset <- pathway_matrix[shared_pathways, , drop = FALSE]
      
      # Limit to top N by frequency
      if (nrow(pathway_subset) > params$top_n) {
        top_pathways <- names(sort(rowSums(pathway_subset), 
                                  decreasing = TRUE)[1:params$top_n])
        pathway_subset <- pathway_subset[top_pathways, , drop = FALSE]
      }
      
      # Shorten pathway names if needed
      rownames(pathway_subset) <- ifelse(nchar(rownames(pathway_subset)) > 50,
                                         paste0(substr(rownames(pathway_subset), 1, 47), "..."),
                                         rownames(pathway_subset))
      
      pdf(file.path(params$output_dir, "figures", "pathway_comparison_heatmap.pdf"),
          width = 10, height = max(8, nrow(pathway_subset) * 0.2))
      
      # For binary data, use specific settings
      tryCatch({
        pheatmap(pathway_subset,
                color = c("white", "steelblue"),
                breaks = c(-0.1, 0.5, 1.1),  # Explicit breaks for 0 and 1
                scale = "none",  # Don't scale binary data
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                clustering_distance_rows = "binary",  # Binary distance for 0/1 data
                clustering_distance_cols = "binary",
                show_rownames = TRUE,
                fontsize_row = 8,
                main = "Pathway Enrichment Across Groups",
                legend_breaks = c(0, 1),
                legend_labels = c("Not enriched", "Enriched"))
      }, error = function(e) {
        # If clustering fails, try without clustering
        warning("Clustering failed, creating heatmap without clustering: ", e$message)
        pheatmap(pathway_subset,
                color = c("white", "steelblue"),
                breaks = c(-0.1, 0.5, 1.1),
                scale = "none",
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                fontsize_row = 8,
                main = "Pathway Enrichment Across Groups",
                legend_breaks = c(0, 1),
                legend_labels = c("Not enriched", "Enriched"))
      })
      dev.off()
    }
    
    # Create summary table
    pathway_summary <- data.frame(
      Pathway = rownames(pathway_matrix),
      N_groups = rowSums(pathway_matrix),
      Groups = apply(pathway_matrix, 1, function(x) {
        paste(colnames(pathway_matrix)[x == 1], collapse = ", ")
      }),
      stringsAsFactors = FALSE
    )
    pathway_summary <- pathway_summary[order(pathway_summary$N_groups, 
                                             decreasing = TRUE), ]
    
    # Save top shared pathways
    top_shared <- head(pathway_summary[pathway_summary$N_groups >= 2, ], params$top_n)
    
    cat(sprintf("\nTop %d shared pathways:\n", nrow(top_shared)))
    print(top_shared, row.names = FALSE)
    cat("\n")
    
    write.table(pathway_summary,
               file = file.path(params$output_dir, "tables", "pathway_comparison.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Pathways unique to each group
    unique_pathways_per_group <- list()
    for (group in colnames(pathway_matrix)) {
      unique_pw <- rownames(pathway_matrix)[pathway_matrix[, group] == 1 & 
                                            rowSums(pathway_matrix) == 1]
      if (length(unique_pw) > 0) {
        unique_pathways_per_group[[group]] <- unique_pw
      }
    }
    
    if (length(unique_pathways_per_group) > 0) {
      cat("Group-specific pathways:\n")
      for (group in names(unique_pathways_per_group)) {
        cat(sprintf("  %s: %d unique pathways\n", 
                   group, length(unique_pathways_per_group[[group]])))
      }
      cat("\n")
    }
  }
}

################################################################################
# ENHANCED PATHWAY ENRICHMENT WITH RATIOS
################################################################################

if (length(all_enrichments) > 0) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("ENHANCED PATHWAY ANALYSIS\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Create enrichment ratio-based heatmap
  pathway_data_detailed <- list()
  
  for (group in names(all_enrichments)) {
    pathway_data_detailed[[group]] <- list()
    
    for (cluster in names(all_enrichments[[group]])) {
      enrich <- all_enrichments[[group]][[cluster]]
      
      if ("p.adjust" %in% colnames(enrich)) {
        sig_idx <- which(enrich$p.adjust < params$fdr_cutoff)
      } else if ("padj" %in% colnames(enrich)) {
        sig_idx <- which(enrich$padj < params$fdr_cutoff)
      } else {
        sig_idx <- which(enrich$pvalue < params$fdr_cutoff)
      }
      
      if (length(sig_idx) > 0) {
        for (i in sig_idx) {
          pathway <- enrich$Description[i]
          if (is.na(pathway) || nchar(pathway) == 0) next
          
          # Calculate enrichment ratio
          if ("GeneRatio" %in% colnames(enrich) && "BgRatio" %in% colnames(enrich)) {
            gene_ratio <- as.numeric(strsplit(as.character(enrich$GeneRatio[i]), "/")[[1]])
            bg_ratio <- as.numeric(strsplit(as.character(enrich$BgRatio[i]), "/")[[1]])
            if (length(gene_ratio) == 2 && length(bg_ratio) == 2) {
              enrich_ratio <- (gene_ratio[1] / gene_ratio[2]) / (bg_ratio[1] / bg_ratio[2])
            } else {
              enrich_ratio <- 1
            }
          } else if ("Count" %in% colnames(enrich)) {
            enrich_ratio <- enrich$Count[i]
          } else {
            enrich_ratio <- 1
          }
          
          # Get p-value
          if ("p.adjust" %in% colnames(enrich)) {
            pval <- enrich$p.adjust[i]
          } else if ("padj" %in% colnames(enrich)) {
            pval <- enrich$padj[i]
          } else {
            pval <- enrich$pvalue[i]
          }
          
          if (!pathway %in% names(pathway_data_detailed[[group]])) {
            pathway_data_detailed[[group]][[pathway]] <- list()
          }
          
          # Extract cluster pattern
          cluster_pattern <- sub("^[0-9]+_", "", cluster)
          
          pathway_data_detailed[[group]][[pathway]][[cluster]] <- list(
            enrichment_ratio = enrich_ratio,
            pvalue = pval,
            cluster = cluster,
            pattern = cluster_pattern
          )
        }
      }
    }
  }
  
  # Create enrichment ratio heatmap
  all_pw_names <- unique(unlist(lapply(pathway_data_detailed, names)))
  all_pw_names <- all_pw_names[!is.na(all_pw_names)]
  
  if (length(all_pw_names) > 0) {
    enrich_ratio_matrix <- matrix(NA, nrow = length(all_pw_names),
                                  ncol = length(names(pathway_data_detailed)))
    rownames(enrich_ratio_matrix) <- all_pw_names
    colnames(enrich_ratio_matrix) <- names(pathway_data_detailed)
    
    for (group in names(pathway_data_detailed)) {
      for (pathway in names(pathway_data_detailed[[group]])) {
        cluster_data <- pathway_data_detailed[[group]][[pathway]]
        avg_ratio <- mean(sapply(cluster_data, function(x) x$enrichment_ratio), na.rm = TRUE)
        enrich_ratio_matrix[pathway, group] <- avg_ratio
      }
    }
    
    # Filter and create heatmap
    pathway_present <- rowSums(!is.na(enrich_ratio_matrix))
    shared_pw <- rownames(enrich_ratio_matrix)[pathway_present >= 2]
    
    if (length(shared_pw) > 0) {
      enrich_subset <- enrich_ratio_matrix[shared_pw, , drop = FALSE]
      
      if (nrow(enrich_subset) > params$top_n * 2) {
        pw_scores <- rowSums(enrich_subset, na.rm = TRUE)
        top_pw <- names(sort(pw_scores, decreasing = TRUE)[1:(params$top_n * 2)])
        enrich_subset <- enrich_subset[top_pw, , drop = FALSE]
      }
      
      display_names <- rownames(enrich_subset)
      display_names <- ifelse(nchar(display_names) > 60,
                             paste0(substr(display_names, 1, 57), "..."),
                             display_names)
      rownames(enrich_subset) <- display_names
      
      color_palette <- colorRampPalette(c("white", "lightyellow", "yellow", "orange", "red", "darkred"))(100)
      
      pdf(file.path(params$output_dir, "figures", "pathway_enrichment_ratios_heatmap.pdf"),
          width = 12, height = max(10, nrow(enrich_subset) * 0.15))
      
      tryCatch({
        pheatmap(enrich_subset,
                color = color_palette,
                scale = "none",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = TRUE,
                fontsize_row = 7,
                main = "Pathway Enrichment Ratios (White = Not Significant)",
                na_col = "gray90",
                border_color = "gray80",
                cellwidth = 25,
                cellheight = 12)
      }, error = function(e) {
        warning("Enrichment ratio heatmap failed: ", e$message)
      })
      
      dev.off()
      
      cat(sprintf("Created enrichment ratio heatmap for %d pathways\n\n", nrow(enrich_subset)))
    }
    
    # Cluster-level pathway comparison
    cat("Comparing pathways by cluster patterns...\n")
    
    pathway_cluster_patterns <- list()
    
    for (pathway in all_pw_names) {
      group_patterns <- character()
      
      for (group in names(pathway_data_detailed)) {
        if (pathway %in% names(pathway_data_detailed[[group]])) {
          patterns <- sapply(pathway_data_detailed[[group]][[pathway]], function(x) x$pattern)
          group_patterns[group] <- paste(unique(patterns), collapse = ", ")
        }
      }
      
      if (length(group_patterns) >= 2) {
        pathway_cluster_patterns[[pathway]] <- data.frame(
          Pathway = pathway,
          t(group_patterns),
          Concordant = length(unique(group_patterns)) == 1,
          Pattern_Summary = paste(unique(group_patterns), collapse = " | "),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
    
    if (length(pathway_cluster_patterns) > 0) {
      pathway_cluster_comparison <- do.call(rbind, pathway_cluster_patterns)
      rownames(pathway_cluster_comparison) <- NULL
      
      write.table(pathway_cluster_comparison,
                 file = file.path(params$output_dir, "tables", "pathway_cluster_patterns.txt"),
                 sep = "\t", quote = FALSE, row.names = FALSE)
      
      n_concordant <- sum(pathway_cluster_comparison$Concordant)
      n_discordant <- sum(!pathway_cluster_comparison$Concordant)
      
      cat(sprintf("  Pathways in multiple groups: %d\n", nrow(pathway_cluster_comparison)))
      cat(sprintf("  Concordant (same cluster pattern): %d\n", n_concordant))
      cat(sprintf("  Discordant (different cluster patterns): %d\n\n", n_discordant))
    }
  }
  
  # Gene-level regulation comparison
  cat("Comparing gene regulation patterns...\n")
  
  gene_regulation <- list()
  lfc_col <- if (params$method == "limma") "logFC" else "log2FoldChange"
  pval_col <- if (params$method == "limma") "adj.P.Val" else "padj"
  
  for (group in names(all_results)) {
    res <- all_results[[group]]
    
    if (!lfc_col %in% colnames(res) || !pval_col %in% colnames(res)) next
    
    sig_genes <- rownames(res)[which(res[[pval_col]] < params$fdr_cutoff)]
    
    for (gene in sig_genes) {
      lfc <- res[gene, lfc_col]
      
      if ("Pattern" %in% colnames(res)) {
        pattern <- res[gene, "Pattern"]
      } else {
        pattern <- NA
      }
      
      if (!is.na(lfc)) {
        regulation <- if (lfc > 0) "Up" else "Down"
      } else {
        regulation <- NA
      }
      
      if (!gene %in% names(gene_regulation)) {
        gene_regulation[[gene]] <- list()
      }
      
      gene_regulation[[gene]][[group]] <- list(
        lfc = lfc,
        regulation = regulation,
        pattern = pattern
      )
    }
  }
  
  if (length(gene_regulation) > 0) {
    gene_counts <- sapply(gene_regulation, length)
    multi_group_genes <- names(gene_counts[gene_counts >= 2])
    
    if (length(multi_group_genes) > 0) {
      gene_reg_data <- list()
      
      for (gene in multi_group_genes) {
        group_data <- gene_regulation[[gene]]
        regulations <- sapply(group_data, function(x) x$regulation)
        patterns <- sapply(group_data, function(x) x$pattern)
        lfcs <- sapply(group_data, function(x) x$lfc)
        
        reg_concordant <- length(unique(regulations[!is.na(regulations)])) == 1
        pattern_concordant <- length(unique(patterns[!is.na(patterns)])) == 1
        
        gene_reg_data[[gene]] <- data.frame(
          Gene_ID = gene,
          Regulation_Concordant = reg_concordant,
          Pattern_Concordant = pattern_concordant,
          Regulation_Summary = paste(unique(regulations), collapse = " | "),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
        for (group in names(group_data)) {
          gene_reg_data[[gene]][[paste0(group, "_Regulation")]] <- group_data[[group]]$regulation
          gene_reg_data[[gene]][[paste0(group, "_LFC")]] <- round(group_data[[group]]$lfc, 3)
        }
      }
      
      gene_reg_comparison <- do.call(rbind, gene_reg_data)
      rownames(gene_reg_comparison) <- NULL
      
      # Add gene symbols
      gene_symbols <- get_gene_symbols(gene_reg_comparison$Gene_ID, all_results)
      gene_reg_comparison <- cbind(Gene_Symbol = gene_symbols, gene_reg_comparison)
      
      write.table(gene_reg_comparison,
                 file = file.path(params$output_dir, "tables", "gene_regulation_patterns.txt"),
                 sep = "\t", quote = FALSE, row.names = FALSE)
      
      n_reg_concordant <- sum(gene_reg_comparison$Regulation_Concordant)
      n_reg_discordant <- sum(!gene_reg_comparison$Regulation_Concordant)
      
      cat(sprintf("  Genes in multiple groups: %d\n", nrow(gene_reg_comparison)))
      cat(sprintf("  Same regulation direction: %d\n", n_reg_concordant))
      cat(sprintf("  Opposite regulation direction: %d\n\n", n_reg_discordant))
      
      # Find opposite regulation genes
      opposite_reg <- gene_reg_comparison[grepl("Up.*Down|Down.*Up", gene_reg_comparison$Regulation_Summary), ]
      
      if (nrow(opposite_reg) > 0) {
        write.table(opposite_reg,
                   file = file.path(params$output_dir, "tables", "genes_opposite_regulation.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
        
        cat(sprintf("  Genes with opposite regulation: %d (saved separately)\n\n", nrow(opposite_reg)))
      }
    }
  }
}

################################################################################
# CREATE SUMMARY REPORT
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("GENERATING SUMMARY REPORT\n")
cat(rep("=", 80), "\n\n", sep = "")

summary_file <- file.path(params$output_dir, "COMPARISON_SUMMARY.txt")

sink(summary_file)

cat(rep("=", 80), "\n", sep = "")
cat("CROSS-GROUP COMPARISON SUMMARY\n")
cat(rep("=", 80), "\n\n", sep = "")

cat(sprintf("Analysis Date: %s\n", Sys.Date()))
cat(sprintf("Input Directory: %s\n", params$input_dir))
cat(sprintf("Number of Groups: %d\n", length(group_dirs)))
cat(sprintf("Groups Analyzed: %s\n", paste(group_dirs, collapse = ", ")))
cat(sprintf("Method: %s\n", params$method))
cat(sprintf("FDR Cutoff: %.3f\n\n", params$fdr_cutoff))

cat(rep("-", 80), "\n", sep = "")
cat("GENE ANALYSIS\n")
cat(rep("-", 80), "\n\n", sep = "")

cat("Significant genes per group:\n")
print(sig_counts, row.names = FALSE)
cat("\n")

cat("Gene overlap:\n")
print(overlap_summary, row.names = FALSE)
cat("\n")

if (exists("direction_summary")) {
  cat("Direction of change (shared genes):\n")
  print(direction_summary, row.names = FALSE)
  cat("\n")
}

if (exists("pattern_summary")) {
  cat("Temporal pattern concordance:\n")
  print(pattern_summary, row.names = FALSE)
  cat("\n")
  
  if (exists("pattern_freq")) {
    cat("Top pattern combinations:\n")
    print(head(pattern_freq, 10))
    cat("\n")
  }
}

if (exists("pathway_summary") && nrow(pathway_summary) > 0) {
  cat(rep("-", 80), "\n", sep = "")
  cat("PATHWAY ANALYSIS\n")
  cat(rep("-", 80), "\n\n", sep = "")
  
  cat(sprintf("Total unique pathways: %d\n", nrow(pathway_summary)))
  cat(sprintf("Pathways in multiple groups: %d\n", 
             sum(pathway_summary$N_groups >= 2)))
  cat(sprintf("Pathways in all groups: %d\n\n", 
             sum(pathway_summary$N_groups == length(all_enrichments))))
  
  if (nrow(top_shared) > 0) {
    cat(sprintf("Top %d shared pathways:\n", nrow(top_shared)))
    print(top_shared, row.names = FALSE)
    cat("\n")
  }
}

cat(rep("-", 80), "\n", sep = "")
cat("OUTPUT FILES\n")
cat(rep("-", 80), "\n\n", sep = "")

cat("Tables:\n")
cat("  - significant_genes_per_group.txt\n")
cat("  - gene_overlap_summary.txt\n")
if (length(shared_in_all) > 0) {
  cat("  - genes_shared_in_all_groups.txt\n")
}
if (exists("direction_summary")) {
  cat("  - direction_of_change_summary.txt\n")
}
if (exists("pattern_comparison")) {
  cat("  - cluster_pattern_comparison.txt\n")
  cat("  - genes_with_concordant_patterns.txt\n")
  cat("  - genes_with_discordant_patterns.txt\n")
}
if (exists("pathway_summary")) {
  cat("  - pathway_comparison.txt\n")
}
if (file.exists(file.path(params$output_dir, "tables", "pathway_cluster_patterns.txt"))) {
  cat("  - pathway_cluster_patterns.txt\n")
}
if (file.exists(file.path(params$output_dir, "tables", "gene_regulation_patterns.txt"))) {
  cat("  - gene_regulation_patterns.txt\n")
}
if (file.exists(file.path(params$output_dir, "tables", "genes_opposite_regulation.txt"))) {
  cat("  - genes_opposite_regulation.txt\n")
}
cat("\nFigures:\n")
if (length(sig_genes_list) >= 2 && length(sig_genes_list) <= 5) {
  cat("  - gene_overlap_venn.pdf\n")
}
if (length(sig_genes_list) > 2) {
  cat("  - gene_overlap_upset.pdf\n")
}
if (length(complete_genes) > 0) {
  cat("  - shared_genes_lfc_heatmap.pdf (with gene symbols)\n")
}
if (exists("pattern_numeric")) {
  cat("  - cluster_pattern_comparison_heatmap.pdf\n")
}
if (exists("pathway_subset")) {
  cat("  - pathway_comparison_heatmap.pdf (binary)\n")
}
if (file.exists(file.path(params$output_dir, "figures", "pathway_enrichment_ratios_heatmap.pdf"))) {
  cat("  - pathway_enrichment_ratios_heatmap.pdf (NEW: with enrichment ratios)\n")
}

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("COMPARISON ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n", sep = "")

sink()

cat(sprintf("\nSummary report saved to: %s\n", summary_file))
cat(sprintf("All results saved to: %s\n\n", params$output_dir))
