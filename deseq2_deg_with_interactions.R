# DESeq2 differential expression (CLI) with CPM filtering, PDFs, and a SELF-CONTAINED HTML report
# NOW WITH DIFFERENCE-OF-DIFFERENCES (INTERACTION) CONTRASTS!
# ---------------------------------------------------------------------------------
# USAGE (examples):
#   Rscript deseq2_cpm_cli.R \
#     --counts counts.tsv \
#     --design design.tsv \
#     --contrasts contrasts.txt \
#     --outdir results_dir \
#     --sample-col SampleID \
#     --gene-col gene \
#     --min-cpm 1 --min-samples 2 \
#     --sig-alpha 0.05 \
#     --top-heatmap-genes 100 \
#     --pca-color Group,Batch,Treatment \
#     --html report.html
#
# REQUIRED FILES
#   --counts:    Raw counts matrix (rows=genes, columns=samples). Can be .txt/.tsv/.csv or a .zip containing one.
#   --design:    Design with columns: Group, Patient_ID, and a sample ID column.
#                The sample ID column defaults to 'SampleID' (override with --sample-col). Must match counts column names.
#   --contrasts: Two-column text/CSV/TSV file listing comparisons. Either header 'treat,ctrl' (case-insensitive)
#                or just two columns with treat in col1 and ctrl in col2.
#                
#                NEW: Also supports difference-of-differences (interaction) contrasts!
#                Format for diff-of-diff: (GroupA-GroupB)-(GroupC-GroupD) in the 'treat' column, leave 'ctrl' empty or write "."
#                Example contrasts.txt:
#                  treat                                      ctrl
#                  Healthy_ADP-h_1hr                          Healthy_un-stimulate
#                  (Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)    .
#
# OPTIONAL PARAMETERS
#   --pca-color: Comma-separated list of design column names to color PCA plots by. Default: "Group"
#                Example: --pca-color Group,Batch,Treatment will generate 3 PCA plots
#
# NOTES
# * Input to DESeq2 must be RAW integer counts. CPM is used ONLY for per-contrast filtering.
# * Design for simple contrasts: ~ Patient_ID + Group **when estimable**; otherwise falls back to ~ Group.
# * Design for interaction contrasts: ~ Patient_ID + Group (if Patient_ID available and estimable across all 4 groups)
# * Duplicate gene symbols are summed prior to analysis.
# * Outputs:
#     - <outdir>/DESeq2_all_contrasts_by_gene.csv
#     - <outdir>/volcano_plots/volcano_<treat>_vs_<ctrl>.pdf
#     - <outdir>/contrast_tables/deseq2_<treat>_vs_<ctrl>.tsv (full per-contrast results)
#     - <outdir>/report.html (SELF-CONTAINED HTML: PCA(s), volcanoes, UpSet, top-100 heatmaps, summary table)

# ---------------------------------------------------------------------------------
# Simple CLI arg parser (no extra packages required)
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  kv <- list()
  i <- 1
  while (i <= length(a)) {
    ai <- a[i]
    if (startsWith(ai, "--")) {
      k <- substring(ai, 3)
      if (grepl("=", k, fixed = TRUE)) {
        parts <- strsplit(k, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        val <- paste(parts[-1], collapse = "=")
        kv[[key]] <- val
        i <- i + 1
      } else {
        if (i + 1 <= length(a) && !startsWith(a[i + 1], "--")) {
          kv[[k]] <- a[i + 1]
          i <- i + 2
        } else {
          kv[[k]] <- TRUE
          i <- i + 1
        }
      }
    } else {
      i <- i + 1
    }
  }
  return(kv)
}

args <- parse_args()

# Required args
counts_file   <- if (!is.null(args$counts)) args$counts else stop("--counts is required")
design_file   <- if (!is.null(args$design)) args$design else stop("--design is required")
contr_file    <- if (!is.null(args$contrasts)) args$contrasts else stop("--contrasts is required")
outdir        <- if (!is.null(args$outdir)) args$outdir else "deseq2_out"

# Optional args
sample_col    <- if (!is.null(args$`sample-col`)) args$`sample-col` else "SampleID"
gene_col_opt  <- if (!is.null(args$`gene-col`)) args$`gene-col` else NA
min_cpm       <- if (!is.null(args$`min-cpm`)) as.numeric(args$`min-cpm`) else 1
min_samples   <- if (!is.null(args$`min-samples`)) as.integer(args$`min-samples`) else 2
sig_alpha     <- if (!is.null(args$`sig-alpha`)) as.numeric(args$`sig-alpha`) else 0.05
top_heatmap   <- if (!is.null(args$`top-heatmap-genes`)) as.integer(args$`top-heatmap-genes`) else 100
html_report   <- if (!is.null(args$html)) args$html else file.path(outdir, "report.html")

# PCA coloring: can provide multiple comma-separated column names
pca_color_arg <- if (!is.null(args$`pca-color`)) args$`pca-color` else "Group"
pca_color_cols <- unlist(strsplit(pca_color_arg, ",", fixed = TRUE))
pca_color_cols <- trimws(pca_color_cols)  # Remove any whitespace

# ---------------------------------------------------------------------------------
# Libraries (plotly is used for interactive widgets & fallback HTML)
suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(htmltools)
})

# Avoid any attempt to open temporary previews in non-interactive environments
options(viewer = function(...) invisible(NULL))
options(browser = function(...) invisible(NULL))

# ---------------------------------------------------------------------------------
# Helpers
`%||%` <- function(a,b) { if (!is.null(a)) return(a); return(b) }

make_outdir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE) }

read_table_auto <- function(f) {
  is_zip <- endsWith(tolower(f), ".zip")
  dt <- NULL
  if (is_zip) {
    td <- tempdir()
    lst <- tryCatch(unzip(f, list = TRUE), error = function(e) NULL)
    if (is.null(lst) || nrow(lst) == 0) stop("Cannot list contents of zip: ", f)
    pick <- lst$Name[tolower(tools::file_ext(lst$Name)) %in% c("txt","tsv","csv")]
    if (length(pick) == 0) stop("No .txt/.tsv/.csv file found inside zip: ", f)
    extracted <- unzip(f, files = pick[1], exdir = td, overwrite = TRUE)
    dt <- data.table::fread(extracted)
  } else {
    dt <- data.table::fread(f)
  }
  return(as.data.frame(dt))
}

cpm_matrix <- function(mat) {
  if (requireNamespace("edgeR", quietly = TRUE)) return(edgeR::cpm(mat, prior.count = 0))
  lib <- colSums(mat)
  return(sweep(mat, 2, lib, FUN = "/") * 1e6)
}

sanitize_for_filename <- function(x) { x <- gsub("[^A-Za-z0-9_-]+", "_", x); x <- gsub("_+", "_", x); return(x) }

aggregate_duplicate_genes <- function(counts_df, gene_col) {
  if (!(gene_col %in% colnames(counts_df))) stop("counts table must have a gene column named '", gene_col, "'")
  sample_cols <- setdiff(colnames(counts_df), gene_col)
  for (cn in sample_cols) counts_df[[cn]] <- suppressWarnings(as.numeric(as.character(counts_df[[cn]])))
  agg <- aggregate(counts_df[sample_cols], by = list(gene = counts_df[[gene_col]]), FUN = function(x) sum(x, na.rm = TRUE))
  return(agg)
}

make_discrete_palette <- function(n) { grDevices::hcl.colors(n, palette = "Dynamic") }

plotly_volcano <- function(res_df, contrast_key) {
  df <- res_df
  df$padj_use <- ifelse(is.na(df$padj), df$pvalue, df$padj)
  df$neglog10 <- -log10(df$padj_use)
  df <- df[is.finite(df$neglog10) & is.finite(df$log2FoldChange), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  hover <- paste0("<b>", df$gene, "</b><br>log2FC=", round(df$log2FoldChange, 3),
                  "<br>padj=", signif(df$padj, 3), "<br>p=", signif(df$pvalue, 3))
  p <- plot_ly(df, x = ~log2FoldChange, y = ~neglog10, type = "scatter", mode = "markers",
               text = hover, hoverinfo = "text") %>%
    layout(title = paste0("Volcano: ", contrast_key), xaxis = list(title = "log2 fold change"),
           yaxis = list(title = "-log10(padj)"))
  return(p)
}

heatmap_widget <- function(mat, row_labels, title = "Heatmap", col_annotations = NULL) {
  if (requireNamespace("heatmaply", quietly = TRUE)) {
    if (!is.null(col_annotations)) {
      # Ensure annotations are factors for proper legend display
      for (col_name in colnames(col_annotations)) {
        if (!is.factor(col_annotations[[col_name]])) {
          col_annotations[[col_name]] <- factor(col_annotations[[col_name]])
        }
      }
      return(heatmaply::heatmaply(mat, Rowv = TRUE, Colv = TRUE, labRow = row_labels, 
                                  col_side_colors = col_annotations, 
                                  col_side_palette = NULL,  # Use default palette
                                  showticklabels = c(TRUE, TRUE),
                                  main = title))
    } else {
      return(heatmaply::heatmaply(mat, Rowv = TRUE, Colv = TRUE, labRow = row_labels, main = title))
    }
  } else {
    return(plot_ly(z = mat, type = "heatmap") %>% layout(title = title))
  }
}

volcano_plot_pdf <- function(res_df, contrast_key, outdir) {
  df <- res_df
  df$padj_use <- ifelse(is.na(df$padj), df$pvalue, df$padj)
  df$padj_use[df$padj_use <= 0] <- NA
  df$neglog10 <- -log10(df$padj_use)
  df$neglog10[is.infinite(df$neglog10)] <- NA
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10)) +
    geom_point(alpha = 0.6, size = 1) +
    labs(title = paste0("Volcano: ", contrast_key), x = "log2 fold change", y = expression(-log[10](padj))) +
    theme_bw(base_size = 12)
  pdf_file <- file.path(outdir, paste0("volcano_", sanitize_for_filename(contrast_key), ".pdf"))
  grDevices::pdf(pdf_file, width = 6, height = 5)
  print(p)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------------------
# NEW: Parse interaction contrast
# Input: "(GroupA-GroupB)-(GroupC-GroupD)"
# Output: list(type="interaction", groupA="GroupA", groupB="GroupB", groupC="GroupC", groupD="GroupD")
# For simple contrasts: list(type="simple", treat="GroupA", ctrl="GroupB")
parse_contrast_spec <- function(treat_str, ctrl_str) {
  # Clean up inputs
  treat_str <- trimws(treat_str)
  ctrl_str <- trimws(ctrl_str)
  
  # Check if it's an interaction contrast (has parentheses)
  if (grepl("\\(.*\\).*-.*\\(.*\\)", treat_str)) {
    # Pattern: (GroupA-GroupB)-(GroupC-GroupD)
    # Extract the two differences
    pattern <- "\\(([^)]+)\\)-\\(([^)]+)\\)"
    matches <- regmatches(treat_str, regexec(pattern, treat_str))[[1]]
    
    if (length(matches) != 3) {
      stop("Invalid interaction contrast format: ", treat_str, 
           "\nExpected format: (GroupA-GroupB)-(GroupC-GroupD)")
    }
    
    diff1 <- matches[2]  # GroupA-GroupB
    diff2 <- matches[3]  # GroupC-GroupD
    
    # Split each difference
    parts1 <- strsplit(diff1, "-", fixed = TRUE)[[1]]
    parts2 <- strsplit(diff2, "-", fixed = TRUE)[[1]]
    
    if (length(parts1) != 2 || length(parts2) != 2) {
      stop("Invalid interaction contrast format: ", treat_str,
           "\nEach part should have exactly 2 groups separated by '-'")
    }
    
    return(list(
      type = "interaction",
      groupA = trimws(parts1[1]),
      groupB = trimws(parts1[2]),
      groupC = trimws(parts2[1]),
      groupD = trimws(parts2[2])
    ))
  } else {
    # Simple contrast
    return(list(
      type = "simple",
      treat = treat_str,
      ctrl = ctrl_str
    ))
  }
}

# ---------------------------------------------------------------------------------
# Read inputs & dirs
make_outdir(outdir)
plots_dir <- file.path(outdir, "volcano_plots"); make_outdir(plots_dir)
tsv_dir   <- file.path(outdir, "contrast_tables"); make_outdir(tsv_dir)

message("Reading counts: ", counts_file)
counts_raw <- read_table_auto(counts_file)

# Determine gene column name
if (!is.na(gene_col_opt)) {
  gene_col <- gene_col_opt
} else {
  potential_gene_cols <- c("gene", "Gene", "gene_id", "gene_name", "GeneID", "GeneName", "GENE", "symbol", "Symbol")
  found_gene_col <- potential_gene_cols[potential_gene_cols %in% colnames(counts_raw)]
  if (length(found_gene_col) == 0) {
    guess <- colnames(counts_raw)[1]
    message("No gene column found among common names. Using first column as gene column: ", guess)
    gene_col <- guess
  } else {
    gene_col <- found_gene_col[1]
    message("Detected gene column: ", gene_col)
  }
}

# Aggregate duplicates
message("Aggregating duplicate genes...")
counts_agg <- aggregate_duplicate_genes(counts_raw, gene_col)
rownames(counts_agg) <- counts_agg$gene
counts_agg$gene <- NULL

message("Reading design: ", design_file)
design_df <- read_table_auto(design_file)
if (!(sample_col %in% colnames(design_df))) stop("Design file must have a column named '", sample_col, "'")
if (!"Group" %in% colnames(design_df)) stop("Design file must have a column named 'Group'")
has_patient <- "Patient_ID" %in% colnames(design_df)
if (!has_patient) message("No Patient_ID found. Will use ~ Group design.")

message("Reading contrasts: ", contr_file)
contrasts_raw <- read_table_auto(contr_file)
if (ncol(contrasts_raw) < 2) stop("Contrasts file must have at least 2 columns (treat, control).")
if (any(tolower(colnames(contrasts_raw)) == "treat")) {
  colnames(contrasts_raw)[tolower(colnames(contrasts_raw)) == "treat"] <- "treat"
}
if (any(tolower(colnames(contrasts_raw)) == "ctrl")) {
  colnames(contrasts_raw)[tolower(colnames(contrasts_raw)) == "ctrl"] <- "ctrl"
}
if (!"treat" %in% colnames(contrasts_raw)) colnames(contrasts_raw)[1] <- "treat"
if (!"ctrl" %in% colnames(contrasts_raw)) colnames(contrasts_raw)[2] <- "ctrl"

# Align samples
common_samples <- intersect(design_df[[sample_col]], colnames(counts_agg))
if (length(common_samples) < 2) stop("Fewer than 2 common samples between counts and design.")
design_df <- design_df[design_df[[sample_col]] %in% common_samples, , drop = FALSE]
counts_agg <- counts_agg[, common_samples, drop = FALSE]
rownames(design_df) <- design_df[[sample_col]]
design_df <- design_df[colnames(counts_agg), , drop = FALSE]

# Validate PCA color columns
valid_pca_cols <- character()
for (col in pca_color_cols) {
  if (col %in% colnames(design_df)) {
    valid_pca_cols <- c(valid_pca_cols, col)
    message("Will generate PCA colored by: ", col)
  } else {
    message("Warning: Column '", col, "' not found in design file. Skipping.")
  }
}

if (length(valid_pca_cols) == 0) {
  message("No valid PCA color columns found. Defaulting to 'Group'.")
  if ("Group" %in% colnames(design_df)) {
    valid_pca_cols <- "Group"
  } else {
    stop("Cannot find 'Group' column in design file and no valid --pca-color columns provided.")
  }
}

counts_mat <- as.matrix(counts_agg)
mode(counts_mat) <- "integer"

cpm_mat <- cpm_matrix(counts_mat)

# ---------------------------------------------------------------------------------
# Contrast loop
message("Running DESeq2 for ", nrow(contrasts_raw), " contrasts...")

results_list <- list()
sig_summary_df <- data.frame(contrast = character(), sig_genes = integer(), stringsAsFactors = FALSE)

for (r in seq_len(nrow(contrasts_raw))) {
  treat_str <- as.character(contrasts_raw[r, "treat"])
  ctrl_str <- as.character(contrasts_raw[r, "ctrl"])
  
  # Parse contrast specification
  contrast_spec <- tryCatch(
    parse_contrast_spec(treat_str, ctrl_str),
    error = function(e) {
      message(sprintf("  [%d/%d] Error parsing contrast: %s", r, nrow(contrasts_raw), e$message))
      return(NULL)
    }
  )
  
  if (is.null(contrast_spec)) next
  
  # ===== SIMPLE CONTRAST =====
  if (contrast_spec$type == "simple") {
    treat_g <- contrast_spec$treat
    ctrl_g <- contrast_spec$ctrl
    contrast_key <- paste(treat_g, "vs", ctrl_g, sep = "_")
    
    idx_treat <- which(design_df$Group == treat_g)
    idx_ctrl <- which(design_df$Group == ctrl_g)
    if (length(idx_treat) == 0 || length(idx_ctrl) == 0) {
      message(sprintf("  [%d/%d] Skipping %s: one or both groups empty.", r, nrow(contrasts_raw), contrast_key))
      next
    }
    
    idx <- c(idx_treat, idx_ctrl)
    sub_counts <- counts_mat[, idx, drop = FALSE]
    sub_cpm <- cpm_mat[, idx, drop = FALSE]
    sub_design <- design_df[idx, , drop = FALSE]
    sub_design$Group <- droplevels(factor(sub_design$Group, levels = c(ctrl_g, treat_g)))
    
    keep <- rowSums(sub_cpm >= min_cpm) >= min_samples
    sub_counts <- sub_counts[keep, , drop = FALSE]
    if (nrow(sub_counts) < 10) {
      message(sprintf("  [%d/%d] Skipping %s: too few genes after filtering (%d).", r, nrow(contrasts_raw), contrast_key, nrow(sub_counts)))
      next
    }
    
    use_patient <- FALSE
    if (has_patient) {
      tbl <- table(sub_design$Patient_ID, sub_design$Group)
      ok <- all(rowSums(tbl > 0) == 2)
      if (ok) use_patient <- TRUE
    }
    
    # Try with Patient_ID first if applicable, fall back to ~ Group if model is not full rank
    dds <- NULL
    if (use_patient) {
      sub_design$Patient_ID <- factor(sub_design$Patient_ID)
      dds <- tryCatch({
        DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_design, design = ~ Patient_ID + Group)
      }, error = function(e) {
        if (grepl("not full rank|linear combination", e$message, ignore.case = TRUE)) {
          message("    Note: Patient_ID blocking not possible (confounded with Group). Using ~ Group instead.")
          return(NULL)
        } else {
          stop(e)
        }
      })
    }
    
    # If Patient_ID model failed or wasn't attempted, use simple model
    if (is.null(dds)) {
      dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_design, design = ~ Group)
    }
    
    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) { message("DESeq failed: ", e$message); NULL })
    if (is.null(dds)) next
    
    res <- tryCatch(results(dds, contrast = c("Group", treat_g, ctrl_g)), error = function(e) { message("results() failed: ", e$message); NULL })
    if (is.null(res)) next
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$contrast <- contrast_key
    results_list[[contrast_key]] <- res_df
    
    nsig <- sum(!is.na(res_df$padj) & res_df$padj <= sig_alpha)
    sig_summary_df <- rbind(sig_summary_df, data.frame(contrast = contrast_key, sig_genes = nsig, stringsAsFactors = FALSE))
    
    volcano_plot_pdf(res_df, contrast_key, plots_dir)
    tsv_path <- file.path(tsv_dir, paste0("deseq2_", sanitize_for_filename(contrast_key), ".tsv"))
    data.table::fwrite(res_df, tsv_path, sep = "\t")
    
    message(sprintf("  [%d/%d] %s: %d sig genes (padj <= %.3f)", r, nrow(contrasts_raw), contrast_key, nsig, sig_alpha))
    
  # ===== INTERACTION CONTRAST (DIFFERENCE-OF-DIFFERENCES) =====
  } else if (contrast_spec$type == "interaction") {
    groupA <- contrast_spec$groupA
    groupB <- contrast_spec$groupB
    groupC <- contrast_spec$groupC
    groupD <- contrast_spec$groupD
    
    contrast_key <- sprintf("(%s-%s)-(%s-%s)", groupA, groupB, groupC, groupD)
    
    # Validate and match group names to actual design groups
    # Helper function to find closest matching group
    find_matching_group <- function(name, available_groups) {
      # Try exact match first
      if (name %in% available_groups) return(name)
      
      # Try case-insensitive match
      idx <- match(tolower(name), tolower(available_groups))
      if (!is.na(idx)) return(available_groups[idx])
      
      # Try with different separators (-, _, space)
      name_normalized <- gsub("[-_ ]+", "", tolower(name))
      for (g in available_groups) {
        g_normalized <- gsub("[-_ ]+", "", tolower(g))
        if (name_normalized == g_normalized) return(g)
      }
      
      return(NULL)
    }
    
    available_groups <- unique(as.character(design_df$Group))
    
    groupA_matched <- find_matching_group(groupA, available_groups)
    groupB_matched <- find_matching_group(groupB, available_groups)
    groupC_matched <- find_matching_group(groupC, available_groups)
    groupD_matched <- find_matching_group(groupD, available_groups)
    
    if (is.null(groupA_matched) || is.null(groupB_matched) || is.null(groupC_matched) || is.null(groupD_matched)) {
      message(sprintf("  [%d/%d] Skipping %s: cannot find all groups in design.", r, nrow(contrasts_raw), contrast_key))
      message("    Specified groups: A='", groupA, "', B='", groupB, "', C='", groupC, "', D='", groupD, "'")
      message("    Matched to: A=", ifelse(is.null(groupA_matched), "NOT FOUND", groupA_matched),
              ", B=", ifelse(is.null(groupB_matched), "NOT FOUND", groupB_matched),
              ", C=", ifelse(is.null(groupC_matched), "NOT FOUND", groupC_matched),
              ", D=", ifelse(is.null(groupD_matched), "NOT FOUND", groupD_matched))
      message("    Available groups in design: ", paste(head(available_groups, 10), collapse = ", "), 
              if(length(available_groups) > 10) "..." else "")
      next
    }
    
    # Use matched names
    groupA <- groupA_matched
    groupB <- groupB_matched
    groupC <- groupC_matched
    groupD <- groupD_matched
    
    # Update contrast key with matched names
    contrast_key <- sprintf("(%s-%s)-(%s-%s)", groupA, groupB, groupC, groupD)
    
    # Find samples for all 4 groups
    idx_A <- which(design_df$Group == groupA)
    idx_B <- which(design_df$Group == groupB)
    idx_C <- which(design_df$Group == groupC)
    idx_D <- which(design_df$Group == groupD)
    
    if (length(idx_A) == 0 || length(idx_B) == 0 || length(idx_C) == 0 || length(idx_D) == 0) {
      message(sprintf("  [%d/%d] Skipping %s: one or more groups empty.", r, nrow(contrasts_raw), contrast_key))
      next
    }
    
    # Combine all samples
    idx <- c(idx_A, idx_B, idx_C, idx_D)
    sub_counts <- counts_mat[, idx, drop = FALSE]
    sub_cpm <- cpm_mat[, idx, drop = FALSE]
    sub_design <- design_df[idx, , drop = FALSE]
    
    # Set Group factor levels
    sub_design$Group <- factor(sub_design$Group, levels = c(groupD, groupC, groupB, groupA))
    
    # CPM filtering across all 4 groups
    keep <- rowSums(sub_cpm >= min_cpm) >= min_samples
    sub_counts <- sub_counts[keep, , drop = FALSE]
    if (nrow(sub_counts) < 10) {
      message(sprintf("  [%d/%d] Skipping %s: too few genes after filtering (%d).", r, nrow(contrasts_raw), contrast_key, nrow(sub_counts)))
      next
    }
    
    # Check if Patient_ID blocking is feasible across all 4 groups
    use_patient <- FALSE
    if (has_patient) {
      tbl <- table(sub_design$Patient_ID, sub_design$Group)
      # Check if each patient has samples in all 4 groups or at least a balanced design
      # For true interaction, ideally patients should span groups, but we'll be flexible
      ok <- nrow(tbl) > 1  # At least 2 patients
      if (ok) use_patient <- TRUE
    }
    
    # Try with Patient_ID first if applicable, fall back to ~ Group if model is not full rank
    dds <- NULL
    if (use_patient) {
      sub_design$Patient_ID <- factor(sub_design$Patient_ID)
      dds <- tryCatch({
        DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_design, design = ~ Patient_ID + Group)
      }, error = function(e) {
        if (grepl("not full rank|linear combination", e$message, ignore.case = TRUE)) {
          message("    Note: Patient_ID blocking not possible (confounded with Group). Using ~ Group instead.")
          return(NULL)
        } else {
          stop(e)
        }
      })
    }
    
    # If Patient_ID model failed or wasn't attempted, use simple model
    if (is.null(dds)) {
      dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_design, design = ~ Group)
    }
    
    dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) { message("DESeq failed: ", e$message); NULL })
    if (is.null(dds)) next
    
    # Build custom contrast vector for (A-B)-(C-D)
    # Get resultsNames to see available coefficients
    rn <- resultsNames(dds)
    
    # In the model ~ Group, the first level is the reference
    # We set levels as c(groupD, groupC, groupB, groupA), so groupD is reference
    # The coefficients are:
    #   GroupC vs GroupD (coef for C)
    #   GroupB vs GroupD (coef for B)  
    #   GroupA vs GroupD (coef for A)
    # 
    # We want: (A-B) - (C-D) 
    #        = (A - groupD) - (B - groupD) - (C - groupD) + (groupD - groupD)
    #        = coef_A - coef_B - coef_C
    
    # Helper function to find coefficient for a group
    # DESeq2 sanitizes names with make.names(), converting -, +, space to .
    find_group_coef <- function(group_name, coef_names, ref_name) {
      # Sanitize the group name the same way DESeq2 does
      group_safe <- make.names(group_name)
      ref_safe <- make.names(ref_name)
      
      # Try to find coefficient with pattern: Group_{group}_vs_{ref}
      # or just Group{group} (older DESeq2 versions)
      
      # Pattern 1: Group_{safe_name}_vs_{ref}
      pattern <- paste0("^Group_", group_safe, "_vs_", ref_safe, "$")
      matches <- grep(pattern, coef_names, value = TRUE, fixed = FALSE)
      
      if (length(matches) == 0) {
        # Pattern 2: Group{safe_name} (no underscore, no vs)
        pattern <- paste0("^Group", group_safe, "$")
        matches <- grep(pattern, coef_names, value = TRUE, fixed = FALSE)
      }
      
      if (length(matches) == 0) {
        # Pattern 3: More flexible - just look for the group name anywhere
        pattern <- paste0("Group.*", group_safe)
        matches <- grep(pattern, coef_names, value = TRUE, fixed = FALSE)
      }
      
      if (length(matches) > 0) {
        if (length(matches) > 1) {
          # Multiple matches - try to find the best one
          # Prefer ones that include _vs_
          vs_matches <- grep("_vs_", matches, value = TRUE)
          if (length(vs_matches) > 0) return(vs_matches[1])
        }
        return(matches[1])
      }
      
      return(NULL)
    }
    
    group_coefs <- grep("^Group", rn, value = TRUE)
    
    coef_A <- find_group_coef(groupA, rn, groupD)
    coef_B <- find_group_coef(groupB, rn, groupD)
    coef_C <- find_group_coef(groupC, rn, groupD)
    
    if (is.null(coef_A) || is.null(coef_B) || is.null(coef_C)) {
      message(sprintf("  [%d/%d] Skipping %s: cannot find all required coefficients in model.", r, nrow(contrasts_raw), contrast_key))
      message("    Looking for groups: A='", groupA, "', B='", groupB, "', C='", groupC, "' (ref: '", groupD, "')")
      message("    Sanitized names: A='", make.names(groupA), "', B='", make.names(groupB), 
              "', C='", make.names(groupC), "', D='", make.names(groupD), "'")
      message("    Found coefficients: A=", ifelse(is.null(coef_A), "NULL", coef_A), 
              ", B=", ifelse(is.null(coef_B), "NULL", coef_B),
              ", C=", ifelse(is.null(coef_C), "NULL", coef_C))
      message("    Available coefficients: ", paste(rn, collapse = ", "))
      next
    }
    
    message("    Using coefficients: +1*(", coef_A, ") -1*(", coef_B, ") -1*(", coef_C, ")")
    
    # Create numeric contrast vector
    contrast_vec <- numeric(length(rn))
    names(contrast_vec) <- rn
    contrast_vec[coef_A] <- 1
    contrast_vec[coef_B] <- -1
    contrast_vec[coef_C] <- -1
    
    res <- tryCatch(
      results(dds, contrast = contrast_vec),
      error = function(e) { 
        message("results() failed for interaction contrast: ", e$message)
        NULL 
      }
    )
    
    if (is.null(res)) next
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$contrast <- contrast_key
    results_list[[contrast_key]] <- res_df
    
    nsig <- sum(!is.na(res_df$padj) & res_df$padj <= sig_alpha)
    sig_summary_df <- rbind(sig_summary_df, data.frame(contrast = contrast_key, sig_genes = nsig, stringsAsFactors = FALSE))
    
    volcano_plot_pdf(res_df, contrast_key, plots_dir)
    tsv_path <- file.path(tsv_dir, paste0("deseq2_", sanitize_for_filename(contrast_key), ".tsv"))
    data.table::fwrite(res_df, tsv_path, sep = "\t")
    
    message(sprintf("  [%d/%d] %s: %d sig genes (padj <= %.3f)", r, nrow(contrasts_raw), contrast_key, nsig, sig_alpha))
  }
}

if (length(results_list) == 0) stop("No contrasts succeeded; nothing to report.")

# Create wide-format table with genes as rows and contrasts as columns
message("Creating combined wide-format results table...")

# Get all unique genes across all contrasts
all_genes <- unique(unlist(lapply(results_list, function(x) x$gene)))

# Initialize wide table with genes
wide_res <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

# Add columns from each contrast
for (contrast_name in names(results_list)) {
  res_df <- results_list[[contrast_name]]
  
  # Select columns to include (exclude the contrast column if it exists)
  cols_to_add <- setdiff(colnames(res_df), c("gene", "contrast"))
  
  # Rename columns with contrast suffix
  contrast_suffix <- paste0("_", contrast_name)
  for (col in cols_to_add) {
    new_col_name <- paste0(col, contrast_suffix)
    wide_res[[new_col_name]] <- res_df[[col]][match(wide_res$gene, res_df$gene)]
  }
}

# Save wide format table
out_csv <- file.path(outdir, "DESeq2_all_contrasts_by_gene_wide.csv")
data.table::fwrite(wide_res, out_csv)
message("Wide-format table saved to: ", out_csv)

# Also save the long format for backwards compatibility
all_res_long <- do.call(rbind, results_list)
out_csv_long <- file.path(outdir, "DESeq2_all_contrasts_by_gene_long.csv")
data.table::fwrite(all_res_long, out_csv_long)
message("Long-format table saved to: ", out_csv_long)

# ---------------------------------------------------------------------------------
# PCA
message("Generating PCA...")
dds_all <- DESeqDataSetFromMatrix(countData = counts_mat, colData = design_df, design = ~ Group)
dds_all <- estimateSizeFactors(dds_all)
vsd <- tryCatch(vst(dds_all, blind = FALSE), error = function(e) NULL)
if (is.null(vsd)) vsd <- tryCatch(rlog(dds_all, blind = FALSE), error = function(e) NULL)
if (is.null(vsd)) stop("vst and rlog both failed.")

vsd_mat <- assay(vsd)
pca_data <- prcomp(t(vsd_mat))
pc_df <- as.data.frame(pca_data$x[, 1:2])
pc_df$Sample <- rownames(pc_df)
var_exp <- (pca_data$sdev^2 / sum(pca_data$sdev^2)) * 100

# Add all PCA color columns to pc_df and convert to character for safe serialization
for (col in valid_pca_cols) {
  pc_df[[col]] <- as.character(design_df[[col]])
}

# ---------------------------------------------------------------------------------
# Save data for report (NOT widgets, just data)
report_data_rds <- file.path(outdir, "report_data.rds")

saveRDS(list(
  sig_summary_df = sig_summary_df,
  vsd_mat = vsd_mat,
  pc_df = pc_df,
  var_exp = var_exp,
  pca_color_cols = valid_pca_cols,
  results_list = results_list  # Include results_list for heatmap generation
), report_data_rds)
message("Saved report data to: ", report_data_rds)

# ---------------------------------------------------------------------------------
# HTML report generation
# We try to render with RMarkdown first. If that fails, fall back to htmlwidgets.
message("Generating HTML report...")

did_rmd <- FALSE

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  # Build RMarkdown content
  rmd_file <- file.path(outdir, "report.Rmd")
  
  # Create PCA sections for each color column
  pca_sections <- paste(sapply(valid_pca_cols, function(col) {
    paste0('
### PCA colored by ', col, '

```{r}
pal <- make_discrete_palette(length(unique(pc_df[["', col, '"]])))
plot_ly(pc_df, x = ~PC1, y = ~PC2, 
        color = ~', col, ', 
        text = ~Sample,
        type = "scatter", mode = "markers",
        colors = pal) %>%
  layout(xaxis = list(title = sprintf("PC1 (%.1f%%)", var_exp[1])),
         yaxis = list(title = sprintf("PC2 (%.1f%%)", var_exp[2])),
         legend = list(title = list(text = "', col, '")))
```
')
  }), collapse = "\n")
  
  # Create volcano sections
  volcano_sections <- paste(sapply(names(results_list), function(nm) {
    paste0('
### ', nm, '

```{r}
tryCatch({
  df <- results_list[["', nm, '"]]
  if (is.null(df) || nrow(df) == 0) {
    htmltools::div("No results available for this contrast.")
  } else {
    df$padj_use <- ifelse(is.na(df$padj), df$pvalue, df$padj)
    df$neglog10 <- -log10(df$padj_use)
    df <- df[is.finite(df$neglog10) & is.finite(df$log2FoldChange), ]
    if (nrow(df) > 0) {
      hover <- paste0("<b>", df$gene, "</b><br>log2FC=", round(df$log2FoldChange, 3),
                      "<br>padj=", signif(df$padj, 3), "<br>p=", signif(df$pvalue, 3))
      plot_ly(df, x = ~log2FoldChange, y = ~neglog10, type = "scatter", mode = "markers",
              text = hover, hoverinfo = "text") %>%
        layout(xaxis = list(title = "log2 fold change"),
               yaxis = list(title = "-log10(padj)"))
    } else {
      htmltools::div("No plottable data after filtering.")
    }
  }
}, error = function(e) {
  htmltools::div("Error generating volcano plot: ", e$message)
})
```
')
  }), collapse = "\n")
  
  # Create heatmap sections
  heatmap_sections <- paste(sapply(names(results_list), function(nm) {
    paste0('
### ', nm, '

```{r}
tryCatch({
  res_df <- results_list[["', nm, '"]]
  if (is.null(res_df) || nrow(res_df) == 0) {
    htmltools::div("No results available for this contrast.")
  } else {
    ord <- res_df$gene[order(res_df$padj, res_df$pvalue, na.last=TRUE)]
    ord <- ord[!is.na(ord)]
    topg <- unique(ord)[seq_len(min(top_n, length(unique(ord))))]
    common <- intersect(topg, rownames(vsd_mat))
    if (length(common) >= 2) {
      mat <- vsd_mat[common, , drop=FALSE]
      mat <- t(scale(t(mat)))
      mat[!is.finite(mat)] <- 0
      
      # Create sample annotations for column colors with proper legend
      sample_names <- colnames(mat)
      sample_groups <- pc_df[sample_names, "Group", drop=TRUE]
      
      # Create a data frame for annotations (this ensures a legend is created)
      col_annotations <- data.frame(
        Group = factor(sample_groups),
        row.names = sample_names
      )
      
      if (requireNamespace("heatmaply", quietly=TRUE)) {
        # Use heatmaply with column annotations and explicit legend
        heatmaply::heatmaply(mat, 
                            Rowv=TRUE, Colv=TRUE, 
                            labRow=common,
                            col_side_colors=col_annotations,
                            col_side_palette=NULL,  # Use default color palette
                            showticklabels=c(TRUE, TRUE),
                            main=sprintf("Top %d genes (%s)", length(common), "', nm, '"))
      } else {
        # Fallback to plotly (no annotations)
        plot_ly(z=mat, type="heatmap") %>% 
          layout(title=sprintf("Top %d genes (%s)", length(common), "', nm, '"))
      }
    } else {
      htmltools::div("Not enough significant genes for heatmap (", length(common), " genes found).")
    }
  }
}, error = function(e) {
  htmltools::div("Error generating heatmap: ", e$message)
})
```
')
  }), collapse = "\n")
  
  # Paths need to be absolute for the RMarkdown document
  rds_path_str <- normalizePath(report_data_rds, mustWork = FALSE)
  tsv_dir_str <- normalizePath(tsv_dir, mustWork = FALSE)
  
  rmd_content <- paste0('---
title: "DESeq2 Differential Expression Report"
output: 
  html_document:
    self_contained: true
    toc: true
    toc_float: true
    code_folding: hide
    theme: flatly
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)

# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(ggplot2)
  library(htmltools)
  if (requireNamespace("heatmaply", quietly=TRUE)) library(heatmaply)
  if (requireNamespace("upsetjs", quietly=TRUE)) library(upsetjs)
  if (requireNamespace("UpSetR", quietly=TRUE)) library(UpSetR)
})

# Set paths (absolute)
rds_path <- "', rds_path_str, '"
tsv_dir <- "', tsv_dir_str, '"
alpha <- ', sig_alpha, '
top_n <- ', top_heatmap, '

# Load data
if (!file.exists(rds_path)) stop(paste("RDS not found:", rds_path))
report_data <- readRDS(rds_path)
sig_summary_df <- report_data$sig_summary_df
vsd_mat <- report_data$vsd_mat
pc_df <- report_data$pc_df
var_exp <- report_data$var_exp
pca_color_cols <- report_data$pca_color_cols
results_list <- report_data$results_list  # Load results_list from RDS

# Convert all PCA color columns to factors
for (col in pca_color_cols) {
  if (col %in% colnames(pc_df)) {
    pc_df[[col]] <- factor(pc_df[[col]])
  }
}

# Helper functions
find_tsvs <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  list.files(dir, pattern="deseq2_.*\\\\.tsv$", full.names=TRUE)
}

make_discrete_palette <- function(n) { 
  grDevices::hcl.colors(n, palette = "Dynamic") 
}
```

## Summary

- **Counts file:** `', counts_file, '`
- **Design file:** `', design_file, '`
- **Contrasts file:** `', contr_file, '`
- **Significance threshold (padj):** ', sig_alpha, '

### Significant genes per contrast

```{r}
if (nrow(sig_summary_df) > 0) {
  knitr::kable(sig_summary_df, caption=sprintf("Significant genes at padj <= %.3f", alpha))
} else {
  htmltools::div("No significant genes at this threshold.")
}
```

## PCA (vst)
', pca_sections, '

## UpSet: overlap of significant genes

```{r}
tsvs <- find_tsvs(tsv_dir)
if (length(tsvs) >= 2) {
  sig_sets <- list()
  for (f in tsvs) {
    df <- tryCatch(data.table::fread(f), error=function(e) NULL)
    if (is.null(df)) next
    bn <- basename(f)
    lab <- sub("deseq2_", "", tools::file_path_sans_ext(bn), fixed=TRUE)
    sig <- df$gene[!is.na(df$padj) & df$padj <= alpha]
    sig_sets[[lab]] <- unique(sig)
  }
  if (length(sig_sets) >= 2) {
    if (requireNamespace("upsetjs", quietly=TRUE)) {
      uj <- upsetjs::upsetjs()
      upsetjs::fromList(uj, sig_sets)
    } else if (requireNamespace("UpSetR", quietly=TRUE)) {
      UpSetR::upset(UpSetR::fromList(sig_sets), nsets=length(sig_sets), order.by="freq")
    } else {
      htmltools::div("Install upsetjs or UpSetR for UpSet plot.")
    }
  } else {
    htmltools::div("Not enough contrasts for UpSet.")
  }
} else {
  htmltools::div("Not enough TSV files for UpSet.")
}
```

## Interactive Volcano Plots
', volcano_sections, '

## Heatmaps: top DE genes per contrast (all samples)
', heatmap_sections, '
')
  
  writeLines(rmd_content, rmd_file)
  message("Created RMarkdown file: ", rmd_file)
  
  tryCatch({
    rmarkdown::render(
      input = rmd_file,
      output_file = basename(html_report),
      output_dir = dirname(html_report),
      quiet = FALSE
    )
    did_rmd <- TRUE
    message("Successfully created HTML report: ", html_report)
  }, error = function(e) {
    message("RMarkdown rendering failed: ", e$message)
  })
}

# Fallback to htmlwidgets if RMarkdown failed
if (!did_rmd) {
  message("Using htmlwidgets fallback method...")
  
  # Convert color columns back to factors for plotting
  pc_df_plot <- pc_df
  for (col in valid_pca_cols) {
    if (col %in% colnames(pc_df_plot)) {
      pc_df_plot[[col]] <- factor(pc_df_plot[[col]])
    }
  }
  
  # Generate PCA widgets for each color column
  pca_widgets <- list()
  for (color_col in valid_pca_cols) {
    if (color_col %in% colnames(pc_df_plot)) {
      pca_widget <- plot_ly(pc_df_plot, x = ~PC1, y = ~PC2, 
                           color = pc_df_plot[[color_col]], 
                           text = ~Sample,
                           type = "scatter", mode = "markers", 
                           colors = make_discrete_palette(length(unique(pc_df_plot[[color_col]])))) %>%
        layout(xaxis = list(title = sprintf("PC1 (%.1f%%)", var_exp[1])),
               yaxis = list(title = sprintf("PC2 (%.1f%%)", var_exp[2])),
               legend = list(title = list(text = color_col)))
      pca_widgets[[color_col]] <- pca_widget
    }
  }
  
  volcano_widgets <- list()
  heatmap_widgets <- list()
  
  for (nm in names(results_list)) {
    res_df <- results_list[[nm]]
    vw <- plotly_volcano(res_df, nm)
    if (!is.null(vw)) volcano_widgets[[nm]] <- vw
    
    ord <- res_df$gene[order(res_df$padj, res_df$pvalue, na.last = TRUE)]
    ord <- ord[!is.na(ord)]
    topg <- unique(ord)[seq_len(min(top_heatmap, length(unique(ord))))]
    common <- intersect(topg, rownames(vsd_mat))
    if (length(common) >= 2) {
      mat <- vsd_mat[common, , drop = FALSE]
      mat <- t(scale(t(mat)))
      mat[!is.finite(mat)] <- 0
      
      # Create sample annotations for heatmap with proper formatting
      sample_names <- colnames(mat)
      sample_groups <- pc_df[sample_names, "Group", drop = TRUE]
      
      # Create data frame with row names for proper legend display
      col_annot <- data.frame(
        Group = factor(sample_groups),
        row.names = sample_names
      )
      
      hw <- heatmap_widget(mat, common, sprintf("Top %d genes (%s)", length(common), nm), 
                          col_annotations = col_annot)
      heatmap_widgets[[nm]] <- hw
    }
  }
  
  # Build summary table
  if (nrow(sig_summary_df) > 0) {
    header <- tags$tr(lapply(colnames(sig_summary_df), tags$th))
    body_rows <- lapply(seq_len(nrow(sig_summary_df)), function(i) {
      tags$tr(lapply(colnames(sig_summary_df), function(nm) tags$td(as.character(sig_summary_df[i, nm]))))
    })
    summary_widget <- tags$table(style = "width:100%;border-collapse:collapse;", tags$thead(header), tags$tbody(body_rows))
  } else {
    summary_widget <- tags$p(sprintf("No significant genes at padj ≤ %.3f", sig_alpha))
  }
  
  page <- browsable(tagList(list(
    tags$h1("DESeq2 Report"),
    tags$p(sprintf("Counts: %s", counts_file)),
    tags$p(sprintf("Design: %s", design_file)),
    tags$p(sprintf("Contrasts: %s", contr_file)),
    tags$h2(sprintf("Summary: significant genes (padj ≤ %.3f)", sig_alpha)),
    summary_widget,
    tags$h2("PCA (vst)"), 
    tagList(lapply(names(pca_widgets), function(nm) tagList(tags$h3(paste("PCA colored by", nm)), pca_widgets[[nm]]))),
    tags$h2("Interactive Volcano Plots"), 
    tagList(lapply(names(volcano_widgets), function(nm) tagList(tags$h3(nm), volcano_widgets[[nm]]))),
    tags$h2(sprintf("Heatmaps: top %d genes per contrast (all samples)", top_heatmap)),
    tagList(lapply(names(heatmap_widgets), function(nm) tagList(tags$h3(nm), heatmap_widgets[[nm]])))
  )))
  
  tryCatch({
    htmlwidgets::saveWidget(page, file = html_report, selfcontained = TRUE)
    message("Created HTML report using htmlwidgets: ", html_report)
  }, error = function(e) {
    message("htmlwidgets saveWidget failed: ", e$message)
    try(htmltools::save_html(page, file = html_report))
  })
}

message("\nAll done!\n- Combined results (wide format): ", file.path(outdir, "DESeq2_all_contrasts_by_gene_wide.csv"),
        "\n- Combined results (long format): ", file.path(outdir, "DESeq2_all_contrasts_by_gene_long.csv"),
        "\n- Volcano plots (PDF): ", plots_dir,
        "\n- Per-contrast TSVs: ", tsv_dir,
        "\n- HTML report: ", html_report)
