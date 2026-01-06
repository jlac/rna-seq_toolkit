#!/usr/bin/env Rscript

# Pre-ranked GSEA Analysis using ClusterProfiler
# Processes all differential expression TSV files in a directory
# Using Gene Ontology: Biological Processes gene sets

# Load required libraries
library(clusterProfiler, quietly = TRUE, warn.conflicts = FALSE)
library(enrichplot, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(DOSE, quietly = TRUE, warn.conflicts = FALSE)
library(ComplexUpset, quietly = TRUE, warn.conflicts = FALSE)
library(ReactomePA, quietly = TRUE, warn.conflicts = FALSE)
library(plotly, quietly = TRUE, warn.conflicts = FALSE)
library(htmltools, quietly = TRUE, warn.conflicts = FALSE)
library(base64enc, quietly = TRUE, warn.conflicts = FALSE)
library(htmlwidgets, quietly = TRUE, warn.conflicts = FALSE)

# Function to perform pre-ranked GSEA
perform_gsea <- function(de_file, 
                        organism = "human",
                        rank_by = "stat",
                        output_dir = "gsea_results",
                        padj_cutoff = 0.25,
                        min_gset_size = 10,
                        max_gset_size = 2000,
                        pathway_database = "GO") {
  
  # Set organism database
  if (organism == "human") {
    library(org.Hs.eg.db, quietly = TRUE, warn.conflicts = FALSE)
    orgdb <- org.Hs.eg.db
    organism_code <- "hsa"
    message("Using human (Homo sapiens) database")
  } else if (organism == "mouse") {
    library(org.Mm.eg.db, quietly = TRUE, warn.conflicts = FALSE)
    orgdb <- org.Mm.eg.db
    organism_code <- "mmu"
    message("Using mouse (Mus musculus) database")
  } else {
    stop("Organism must be either 'human' or 'mouse'")
  }
  
  # Read differential expression data
  message(paste("Reading file:", de_file))
  de_data <- read.delim(de_file, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("gene", rank_by)
  if (!all(required_cols %in% colnames(de_data))) {
    stop(paste("Missing required columns. Need:", paste(required_cols, collapse = ", ")))
  }
  
  # Remove NAs and prepare ranked gene list
  de_data <- de_data[!is.na(de_data[[rank_by]]), ]
  de_data <- de_data[!is.na(de_data$gene), ]
  
  # Create ranked list
  gene_list <- de_data[[rank_by]]
  names(gene_list) <- de_data$gene
  
  # Sort in decreasing order (most upregulated to most downregulated)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  message(paste("Total genes in ranked list:", length(gene_list)))
  
  # Perform GSEA based on chosen pathway database
  if (pathway_database == "GO") {
    message("Running GSEA with GO:BP gene sets...")
    gsea_results <- gseGO(
      geneList = gene_list,
      OrgDb = orgdb,
      ont = "BP",  # Biological Processes
      keyType = "SYMBOL",
      minGSSize = min_gset_size,
      maxGSSize = max_gset_size,
      pvalueCutoff = 1,  # Set to 1 to return all gene sets tested
      pAdjustMethod = "BH",
      verbose = TRUE,
      by = "fgsea"
    )
  } else if (pathway_database == "Reactome") {
    message("Running GSEA with Reactome pathways...")
    
    # Convert gene symbols to Entrez IDs for Reactome
    if (organism == "human") {
      gene_entrez <- mapIds(orgdb, names(gene_list), 'ENTREZID', 'SYMBOL')
    } else if (organism == "mouse") {
      gene_entrez <- mapIds(orgdb, names(gene_list), 'ENTREZID', 'SYMBOL')
    }
    
    # Remove NAs and create new ranked list with Entrez IDs
    gene_entrez <- gene_entrez[!is.na(gene_entrez)]
    gene_list_entrez <- gene_list[names(gene_list) %in% names(gene_entrez)]
    names(gene_list_entrez) <- gene_entrez[names(gene_list_entrez)]
    
    # Remove duplicates (keep highest absolute value)
    gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
    
    message(paste("Converted to", length(gene_list_entrez), "Entrez IDs"))
    
    gsea_results <- gsePathway(
      geneList = gene_list_entrez,
      organism = organism,
      minGSSize = min_gset_size,
      maxGSSize = max_gset_size,
      pvalueCutoff = 1,  # Set to 1 to return all gene sets tested
      pAdjustMethod = "BH",
      verbose = TRUE
    )
    
    # Convert Entrez IDs back to gene symbols in results
    if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
      message("Converting Entrez IDs back to gene symbols in results...")
      
      # Create reverse mapping: Entrez ID -> Symbol
      entrez_to_symbol <- setNames(names(gene_entrez), gene_entrez)
      
      # Convert core_enrichment column from Entrez IDs to Symbols
      gsea_results@result$core_enrichment <- sapply(gsea_results@result$core_enrichment, function(genes_str) {
        if (is.na(genes_str) || genes_str == "") {
          return(genes_str)
        }
        # Split by "/" to get individual Entrez IDs
        entrez_ids <- unlist(strsplit(genes_str, "/"))
        # Convert to symbols
        symbols <- entrez_to_symbol[entrez_ids]
        # Replace any NAs with original Entrez ID
        symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
        # Rejoin with "/"
        paste(symbols, collapse = "/")
      })
      
      message("Gene symbol conversion complete")
    }
  } else if (grepl("\\.gmt$", pathway_database, ignore.case = TRUE)) {
    # GMT file provided
    message(paste("Running GSEA with custom GMT file:", pathway_database))
    
    # Read GMT file
    message("Reading GMT file...")
    gmt_data <- clusterProfiler::read.gmt(pathway_database)
    
    message(paste("Loaded", length(unique(gmt_data$term)), "gene sets from GMT file"))
    message(paste("Total gene-set associations:", nrow(gmt_data)))
    
    # Run GSEA with custom gene sets
    gsea_results <- GSEA(
      geneList = gene_list,
      TERM2GENE = gmt_data,
      minGSSize = min_gset_size,
      maxGSSize = max_gset_size,
      pvalueCutoff = 1,  # Set to 1 to return all gene sets tested
      pAdjustMethod = "BH",
      verbose = TRUE,
      by = "fgsea"
    )
  } else {
    stop("pathway_database must be either 'GO', 'Reactome', or a path to a GMT file (ending in .gmt)")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate base filename from input file
  base_name <- tools::file_path_sans_ext(basename(de_file))
  
  # Save results
  if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
    
    # Save full results table (all gene sets tested)
    results_file_all <- file.path(output_dir, paste0(base_name, "_gsea_all_results.tsv"))
    write.table(as.data.frame(gsea_results), 
                results_file_all, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    message(paste("All results saved to:", results_file_all))
    message(paste("Total gene sets tested:", nrow(gsea_results@result)))
    
    # Filter and save significant results (FDR < 0.1)
    sig_results <- gsea_results@result[gsea_results@result$p.adjust < 0.1, ]
    
    if (nrow(sig_results) > 0) {
      results_file_sig <- file.path(output_dir, paste0(base_name, "_gsea_significant_FDR0.1.tsv"))
      write.table(sig_results, 
                  results_file_sig, 
                  sep = "\t", 
                  row.names = FALSE, 
                  quote = FALSE)
      message(paste("Significant results (FDR < 0.1) saved to:", results_file_sig))
      message(paste("Significant pathways found (FDR < 0.1):", nrow(sig_results)))
    } else {
      message("No significant pathways found at FDR < 0.1")
    }
    
    # Generate plots if significant results exist
    plot_sig_results <- gsea_results@result[gsea_results@result$p.adjust < padj_cutoff, ]
    
    if (nrow(plot_sig_results) > 0) {
      # Dotplot
      tryCatch({
        p1 <- dotplot(gsea_results, showCategory = 20, font.size = 10) +
          ggtitle(paste(pathway_database, "GSEA -", base_name))
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_dotplot.pdf")), 
               p1, width = 10, height = 8)
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_dotplot.png")), 
               p1, width = 10, height = 8, dpi = 300)
        message("Dotplot saved")
      }, error = function(e) {
        message(paste("Could not create dotplot:", e$message))
      })
      
      # Enrichment plots - split by direction
      # Positive NES (upregulated pathways)
      pos_pathways <- plot_sig_results[plot_sig_results$NES > 0, ]
      if (nrow(pos_pathways) > 0) {
        tryCatch({
          # Sort by p.adjust and take top 20 (or fewer)
          pos_pathways <- pos_pathways[order(pos_pathways$p.adjust), ]
          top_pos <- head(pos_pathways, 20)
          
          p2_up <- gseaplot2(gsea_results, 
                         geneSetID = top_pos$ID,
                         pvalue_table = TRUE,
                         title = paste("Top Upregulated Pathways (Positive NES) -", base_name))
          ggsave(file.path(output_dir, paste0(base_name, "_gsea_enrichment_plot_upregulated.pdf")), 
                 p2_up, width = 12, height = 10)
          message(paste("Upregulated enrichment plot saved (", nrow(top_pos), "of", nrow(pos_pathways), "pathways)"))
        }, error = function(e) {
          message(paste("Could not create upregulated enrichment plot:", e$message))
        })
      } else {
        message("No upregulated pathways to plot")
      }
      
      # Negative NES (downregulated pathways)
      neg_pathways <- plot_sig_results[plot_sig_results$NES < 0, ]
      if (nrow(neg_pathways) > 0) {
        tryCatch({
          # Sort by p.adjust and take top 20 (or fewer)
          neg_pathways <- neg_pathways[order(neg_pathways$p.adjust), ]
          top_neg <- head(neg_pathways, 20)
          
          p2_down <- gseaplot2(gsea_results, 
                           geneSetID = top_neg$ID,
                           pvalue_table = TRUE,
                           title = paste("Top Downregulated Pathways (Negative NES) -", base_name))
          ggsave(file.path(output_dir, paste0(base_name, "_gsea_enrichment_plot_downregulated.pdf")), 
                 p2_down, width = 12, height = 10)
          message(paste("Downregulated enrichment plot saved (", nrow(top_neg), "of", nrow(neg_pathways), "pathways)"))
        }, error = function(e) {
          message(paste("Could not create downregulated enrichment plot:", e$message))
        })
      } else {
        message("No downregulated pathways to plot")
      }
      
      # Ridge plot (showing expression distributions)
      tryCatch({
        p3 <- ridgeplot(gsea_results, showCategory = 20) +
          ggtitle(paste(pathway_database, "GSEA Ridge Plot -", base_name))
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_ridgeplot.pdf")), 
               p3, width = 10, height = 8)
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_ridgeplot.png")), 
               p3, width = 10, height = 8, dpi = 300)
        message("Ridge plot saved")
      }, error = function(e) {
        message(paste("Could not create ridge plot:", e$message))
      })
      
      # Treeplot (showing hierarchical clustering of GO terms)
      tryCatch({
        # pairwise_termsim is needed for treeplot
        gsea_results_sim <- pairwise_termsim(gsea_results)
        p4 <- treeplot(gsea_results_sim, showCategory = 30) +
          ggtitle(paste(pathway_database, "GSEA Tree Plot -", base_name))
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_treeplot.pdf")), 
               p4, width = 12, height = 10)
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_treeplot.png")), 
               p4, width = 12, height = 10, dpi = 300)
        message("Tree plot saved")
      }, error = function(e) {
        message(paste("Could not create tree plot:", e$message))
      })
      
      # Enrichment map plot (network of enriched terms)
      tryCatch({
        # Use the same similarity calculation
        if (!exists("gsea_results_sim")) {
          gsea_results_sim <- pairwise_termsim(gsea_results)
        }
        p5 <- emapplot(gsea_results_sim, showCategory = 30) +
          ggtitle(paste(pathway_database, "GSEA Enrichment Map -", base_name))
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_emapplot.pdf")), 
               p5, width = 12, height = 10)
        ggsave(file.path(output_dir, paste0(base_name, "_gsea_emapplot.png")), 
               p5, width = 12, height = 10, dpi = 300)
        message("Enrichment map plot saved")
      }, error = function(e) {
        message(paste("Could not create enrichment map plot:", e$message))
      })
    } else {
      message("No significant pathways to plot")
    }
    
  } else {
    message("No significant GSEA results found")
  }
  
  return(gsea_results)
}

# Function to merge all GSEA results into wide format
merge_gsea_results <- function(results_list, file_names, output_file) {
  
  message("\nMerging all GSEA results into wide format...")
  
  # Extract result data frames and add file identifier
  all_results <- list()
  
  for (i in seq_along(results_list)) {
    if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
      df <- results_list[[i]]@result
      file_name <- tools::file_path_sans_ext(basename(file_names[i]))
      
      # Select key columns
      df_subset <- df[, c("ID", "Description", "NES", "pvalue", "p.adjust")]
      
      # Rename columns to include file name
      colnames(df_subset)[3:5] <- paste0(file_name, "_", c("NES", "pvalue", "padj"))
      
      all_results[[i]] <- df_subset
    }
  }
  
  if (length(all_results) == 0) {
    message("No results to merge")
    return(NULL)
  }
  
  # Merge all results by ID and Description
  merged <- all_results[[1]]
  
  if (length(all_results) > 1) {
    for (i in 2:length(all_results)) {
      merged <- merge(merged, all_results[[i]], 
                     by = c("ID", "Description"), 
                     all = TRUE)
    }
  }
  
  # Sort by ID
  merged <- merged[order(merged$ID), ]
  
  # Write merged results
  write.table(merged, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Merged results saved to:", output_file))
  message(paste("Total unique gene sets:", nrow(merged)))
  
  return(merged)
}

# Function to create comparative dotplots across all analyses
create_comparative_dotplots <- function(results_list, file_names, output_dir) {
  
  message("\nCreating comparative dotplots across all analyses...")
  
  # Extract data from each result
  all_data <- list()
  comparison_names <- c()
  
  for (i in seq_along(results_list)) {
    if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
      df <- results_list[[i]]@result
      # Filter for FDR < 0.1
      df <- df[df$p.adjust < 0.1, ]
      
      if (nrow(df) > 0) {
        file_name <- tools::file_path_sans_ext(basename(file_names[i]))
        comparison_names <- c(comparison_names, file_name)
        all_data[[file_name]] <- df
      }
    }
  }
  
  if (length(all_data) == 0) {
    message("No significant results to create comparative plots")
    return(NULL)
  }
  
  # UPREGULATED PATHWAYS
  tryCatch({
    # Get top 10 upregulated from each comparison
    top_up_pathways <- c()
    for (comp_name in names(all_data)) {
      df <- all_data[[comp_name]]
      up_df <- df[df$NES > 0, ]
      if (nrow(up_df) > 0) {
        up_df <- up_df[order(up_df$p.adjust), ]
        top_up <- head(up_df$ID, 10)
        top_up_pathways <- c(top_up_pathways, top_up)
      }
    }
    top_up_pathways <- unique(top_up_pathways)
    
    if (length(top_up_pathways) > 0) {
      # Create matrix for NES and -log10(padj)
      nes_matrix <- matrix(NA, nrow = length(top_up_pathways), ncol = length(all_data))
      padj_matrix <- matrix(NA, nrow = length(top_up_pathways), ncol = length(all_data))
      rownames(nes_matrix) <- top_up_pathways
      colnames(nes_matrix) <- names(all_data)
      rownames(padj_matrix) <- top_up_pathways
      colnames(padj_matrix) <- names(all_data)
      
      # Get descriptions for pathways - use a list instead of named vector
      desc_map <- list()
      for (comp_name in names(all_data)) {
        df <- all_data[[comp_name]]
        for (pathway in top_up_pathways) {
          if (is.null(desc_map[[pathway]])) {
            idx <- which(df$ID == pathway)
            if (length(idx) > 0) {
              desc_map[[pathway]] <- df$Description[idx[1]]
            }
          }
        }
      }
      
      # Fill matrices
      for (comp_name in names(all_data)) {
        df <- all_data[[comp_name]]
        for (pathway in top_up_pathways) {
          idx <- which(df$ID == pathway)
          if (length(idx) > 0) {
            nes_matrix[pathway, comp_name] <- df$NES[idx[1]]
            padj_matrix[pathway, comp_name] <- -log10(df$p.adjust[idx[1]])
          }
        }
      }
      
      # Convert to long format for ggplot
      plot_data_up <- data.frame()
      for (i in 1:nrow(nes_matrix)) {
        for (j in 1:ncol(nes_matrix)) {
          if (!is.na(nes_matrix[i,j])) {
            pathway_id <- rownames(nes_matrix)[i]
            pathway_desc <- desc_map[[pathway_id]]
            if (is.null(pathway_desc)) pathway_desc <- pathway_id
            
            plot_data_up <- rbind(plot_data_up, data.frame(
              Pathway = pathway_desc,
              Comparison = colnames(nes_matrix)[j],
              NES = nes_matrix[i,j],
              NegLog10Padj = padj_matrix[i,j],
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      if (nrow(plot_data_up) > 0) {
        # Create dotplot
        p_up <- ggplot(plot_data_up, aes(x = Comparison, y = Pathway)) +
          geom_point(aes(size = NegLog10Padj, color = NES)) +
          scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          scale_size_continuous(range = c(2, 10), name = "-log10(padj)") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(size = 9),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "right") +
          labs(title = "Top 10 Upregulated Pathways Across All Comparisons (FDR<0.1)",
               x = "Comparison", y = "Pathway", color = "NES")
        
        ggsave(file.path(output_dir, "comparative_dotplot_upregulated_top10.pdf"), 
               p_up, width = max(8, length(all_data) * 1.5), height = max(8, length(top_up_pathways) * 0.3))
        ggsave(file.path(output_dir, "comparative_dotplot_upregulated_top10.png"), 
               p_up, width = max(8, length(all_data) * 1.5), height = max(8, length(top_up_pathways) * 0.3), dpi = 300)
        message(paste("Comparative upregulated dotplot saved (", length(top_up_pathways), "unique pathways)"))
      }
    }
  }, error = function(e) {
    message(paste("Could not create comparative upregulated dotplot:", e$message))
  })
  
  # DOWNREGULATED PATHWAYS
  tryCatch({
    # Get top 10 downregulated from each comparison
    top_down_pathways <- c()
    for (comp_name in names(all_data)) {
      df <- all_data[[comp_name]]
      down_df <- df[df$NES < 0, ]
      if (nrow(down_df) > 0) {
        down_df <- down_df[order(down_df$p.adjust), ]
        top_down <- head(down_df$ID, 10)
        top_down_pathways <- c(top_down_pathways, top_down)
      }
    }
    top_down_pathways <- unique(top_down_pathways)
    
    if (length(top_down_pathways) > 0) {
      # Create matrix for NES and -log10(padj)
      nes_matrix <- matrix(NA, nrow = length(top_down_pathways), ncol = length(all_data))
      padj_matrix <- matrix(NA, nrow = length(top_down_pathways), ncol = length(all_data))
      rownames(nes_matrix) <- top_down_pathways
      colnames(nes_matrix) <- names(all_data)
      rownames(padj_matrix) <- top_down_pathways
      colnames(padj_matrix) <- names(all_data)
      
      # Get descriptions for pathways - use a list instead of named vector
      desc_map <- list()
      for (comp_name in names(all_data)) {
        df <- all_data[[comp_name]]
        for (pathway in top_down_pathways) {
          if (is.null(desc_map[[pathway]])) {
            idx <- which(df$ID == pathway)
            if (length(idx) > 0) {
              desc_map[[pathway]] <- df$Description[idx[1]]
            }
          }
        }
      }
      
      # Fill matrices
      for (comp_name in names(all_data)) {
        df <- all_data[[comp_name]]
        for (pathway in top_down_pathways) {
          idx <- which(df$ID == pathway)
          if (length(idx) > 0) {
            nes_matrix[pathway, comp_name] <- df$NES[idx[1]]
            padj_matrix[pathway, comp_name] <- -log10(df$p.adjust[idx[1]])
          }
        }
      }
      
      # Convert to long format for ggplot
      plot_data_down <- data.frame()
      for (i in 1:nrow(nes_matrix)) {
        for (j in 1:ncol(nes_matrix)) {
          if (!is.na(nes_matrix[i,j])) {
            pathway_id <- rownames(nes_matrix)[i]
            pathway_desc <- desc_map[[pathway_id]]
            if (is.null(pathway_desc)) pathway_desc <- pathway_id
            
            plot_data_down <- rbind(plot_data_down, data.frame(
              Pathway = pathway_desc,
              Comparison = colnames(nes_matrix)[j],
              NES = nes_matrix[i,j],
              NegLog10Padj = padj_matrix[i,j]
            ))
          }
        }
      }
      
      if (nrow(plot_data_down) > 0) {
        # Create dotplot
        p_down <- ggplot(plot_data_down, aes(x = Comparison, y = Pathway)) +
          geom_point(aes(size = NegLog10Padj, color = NES)) +
          scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          scale_size_continuous(range = c(2, 10), name = "-log10(padj)") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(size = 9),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "right") +
          labs(title = "Top 10 Downregulated Pathways Across All Comparisons (FDR<0.1)",
               x = "Comparison", y = "Pathway", color = "NES")
        
        ggsave(file.path(output_dir, "comparative_dotplot_downregulated_top10.pdf"), 
               p_down, width = max(8, length(all_data) * 1.5), height = max(8, length(top_down_pathways) * 0.3))
        ggsave(file.path(output_dir, "comparative_dotplot_downregulated_top10.png"), 
               p_down, width = max(8, length(all_data) * 1.5), height = max(8, length(top_down_pathways) * 0.3), dpi = 300)
        message(paste("Comparative downregulated dotplot saved (", length(top_down_pathways), "unique pathways)"))
      }
    }
  }, error = function(e) {
    message(paste("Could not create comparative downregulated dotplot:", e$message))
  })
}

# ComplexUpset-based upset plot function
# Replace the create_upset_plot function in preranked_gsea_analysis.R with this

create_upset_plot <- function(results_list, file_names, output_dir, padj_cutoff) {
  
  message("\nCreating upset plots for pathway overlap...")
  
  # Extract significant pathways from each comparison, split by direction
  pathway_sets_up <- list()
  pathway_sets_down <- list()
  
  for (i in seq_along(results_list)) {
    if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
      df <- results_list[[i]]@result
      sig_df <- df[df$p.adjust < padj_cutoff, ]
      
      if (nrow(sig_df) > 0) {
        file_name <- tools::file_path_sans_ext(basename(file_names[i]))
        up_pathways <- sig_df$ID[sig_df$NES > 0]
        down_pathways <- sig_df$ID[sig_df$NES < 0]
        
        if (length(up_pathways) > 0) {
          pathway_sets_up[[file_name]] <- up_pathways
        }
        if (length(down_pathways) > 0) {
          pathway_sets_down[[file_name]] <- down_pathways
        }
      }
    }
  }
  
  # UPREGULATED PATHWAYS UPSET PLOT
  if (length(pathway_sets_up) >= 2) {
    tryCatch({
      message(paste("Creating ComplexUpset plot with", length(unique(unlist(pathway_sets_up))), "upregulated pathways"))
      
      all_pathways_up <- unique(unlist(pathway_sets_up))
      upset_data <- data.frame(Pathway = all_pathways_up, stringsAsFactors = FALSE)
      
      # Use actual comparison names
      comp_names <- names(pathway_sets_up)
      for (i in seq_along(pathway_sets_up)) {
        upset_data[[comp_names[i]]] <- upset_data$Pathway %in% pathway_sets_up[[i]]
      }
      
      message("ComplexUpset data structure:")
      print(head(upset_data))
      
      p <- upset(
        upset_data,
        intersect = comp_names,
        name = "Pathway Overlap",
        width_ratio = 0.1,
        min_size = 0,
        sort_sets = FALSE,
        sort_intersections_by = "cardinality",
        base_annotations = list(
          'Intersection size' = intersection_size(text = list(size = 3))
        ),
        set_sizes = (
          upset_set_size()
          + geom_text(aes(label = ..count..), hjust = 1.1, stat = "count", size = 3.5)
          + theme(axis.text.x = element_text(angle = 90))
        ),
        themes = upset_default_themes(text = element_text(size = 9))
      ) + ggtitle("Upregulated Pathways Overlap") +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      ggsave(file.path(output_dir, "upset_plot_upregulated_pathways.pdf"), p, width = 16, height = 8, device = "pdf")
      ggsave(file.path(output_dir, "upset_plot_upregulated_pathways.png"), p, width = 16, height = 8, dpi = 300, device = "png")
      
      message("Upregulated upset plot completed successfully")
      
    }, error = function(e) {
      message(paste("ERROR creating upregulated upset plot:", e$message))
      print(e)
    })
  } else {
    message("Need at least 2 comparisons with upregulated pathways for upset plot")
  }
  
  # DOWNREGULATED PATHWAYS UPSET PLOT
  if (length(pathway_sets_down) >= 2) {
    tryCatch({
      message(paste("Creating ComplexUpset plot with", length(unique(unlist(pathway_sets_down))), "downregulated pathways"))
      
      all_pathways_down <- unique(unlist(pathway_sets_down))
      upset_data <- data.frame(Pathway = all_pathways_down, stringsAsFactors = FALSE)
      
      # Use actual comparison names
      comp_names <- names(pathway_sets_down)
      for (i in seq_along(pathway_sets_down)) {
        upset_data[[comp_names[i]]] <- upset_data$Pathway %in% pathway_sets_down[[i]]
      }
      
      p <- upset(
        upset_data,
        intersect = comp_names,
        name = "Pathway Overlap",
        width_ratio = 0.1,
        min_size = 0,
        sort_sets = FALSE,
        sort_intersections_by = "cardinality",
        base_annotations = list(
          'Intersection size' = intersection_size(text = list(size = 3))
        ),
        set_sizes = (
          upset_set_size()
          + geom_text(aes(label = ..count..), hjust = 1.1, stat = "count", size = 3.5)
          + theme(axis.text.x = element_text(angle = 90))
        ),
        themes = upset_default_themes(text = element_text(size = 9))
      ) + ggtitle("Downregulated Pathways Overlap") +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      ggsave(file.path(output_dir, "upset_plot_downregulated_pathways.pdf"), p, width = 16, height = 8, device = "pdf")
      ggsave(file.path(output_dir, "upset_plot_downregulated_pathways.png"), p, width = 16, height = 8, dpi = 300, device = "png")
      
      message("Downregulated upset plot completed successfully")
      
    }, error = function(e) {
      message(paste("ERROR creating downregulated upset plot:", e$message))
      print(e)
    })
  } else {
    message("Need at least 2 comparisons with downregulated pathways for upset plot")
  }
  
  # Save summary statistics
  if (length(pathway_sets_up) > 0 || length(pathway_sets_down) > 0) {
    overlap_summary <- data.frame(
      Comparison = unique(c(names(pathway_sets_up), names(pathway_sets_down))),
      NumUpregulated = 0,
      NumDownregulated = 0,
      stringsAsFactors = FALSE
    )
    
    for (comp_name in overlap_summary$Comparison) {
      if (comp_name %in% names(pathway_sets_up)) {
        overlap_summary$NumUpregulated[overlap_summary$Comparison == comp_name] <- length(pathway_sets_up[[comp_name]])
      }
      if (comp_name %in% names(pathway_sets_down)) {
        overlap_summary$NumDownregulated[overlap_summary$Comparison == comp_name] <- length(pathway_sets_down[[comp_name]])
      }
    }
    
    write.table(overlap_summary, file.path(output_dir, "pathway_overlap_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    message("Overlap summary statistics saved")
  }
}

# Function to perform PCA on NES values across comparisons
create_nes_pca <- function(results_list, file_names, output_dir, padj_cutoff = 0.25) {
  
  message("\nPerforming PCA on NES values across comparisons...")
  
  # Extract all pathway IDs and their NES values from each comparison
  nes_data_list <- list()
  
  for (i in seq_along(results_list)) {
    if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
      comp_name <- tools::file_path_sans_ext(basename(file_names[i]))
      df <- results_list[[i]]@result
      
      # Get NES values for all pathways (not just significant ones for PCA)
      nes_data_list[[comp_name]] <- data.frame(
        ID = df$ID,
        NES = df$NES,
        padj = df$p.adjust,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(nes_data_list) < 2) {
    message("Need at least 2 comparisons for PCA")
    return(NULL)
  }
  
  # Get all unique pathways across all comparisons
  all_pathways <- unique(unlist(lapply(nes_data_list, function(x) x$ID)))
  
  # Create matrix: rows = pathways, columns = comparisons
  nes_matrix <- matrix(0, nrow = length(all_pathways), ncol = length(nes_data_list))
  rownames(nes_matrix) <- all_pathways
  colnames(nes_matrix) <- names(nes_data_list)
  
  # Fill matrix with NES values
  for (comp_name in names(nes_data_list)) {
    df <- nes_data_list[[comp_name]]
    nes_matrix[df$ID, comp_name] <- df$NES
  }
  
  # Remove pathways with all zeros (not found in any comparison)
  nes_matrix <- nes_matrix[rowSums(abs(nes_matrix)) > 0, , drop = FALSE]
  
  message(paste("PCA matrix dimensions:", nrow(nes_matrix), "pathways ×", ncol(nes_matrix), "comparisons"))
  
  # Perform PCA on transposed matrix (samples = comparisons)
  pca_result <- prcomp(t(nes_matrix), scale. = TRUE, center = TRUE)
  
  # Get variance explained
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  # Create PCA plot of comparisons (scores plot)
  pca_scores <- as.data.frame(pca_result$x[, 1:min(5, ncol(pca_result$x))])
  pca_scores$Comparison <- rownames(pca_scores)
  
  # Static PCA biplot - comparisons (for PDF)
  p1_static <- ggplot(pca_scores, aes(x = PC1, y = PC2, label = Comparison)) +
    geom_point(size = 4, color = "steelblue") +
    geom_text(hjust = -0.1, vjust = -0.5, size = 3) +
    labs(
      title = "PCA of Comparisons Based on Pathway NES",
      x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(file.path(output_dir, "pca_comparisons_biplot.pdf"), p1_static, width = 10, height = 8)
  ggsave(file.path(output_dir, "pca_comparisons_biplot.png"), p1_static, width = 10, height = 8, dpi = 300)
  
  # Interactive PCA biplot using plotly
  p1_interactive <- plot_ly(data = pca_scores, 
                            x = ~PC1, 
                            y = ~PC2,
                            type = 'scatter',
                            mode = 'markers+text',
                            text = ~Comparison,
                            textposition = 'top right',
                            marker = list(size = 12, color = 'steelblue', opacity = 0.8),
                            hovertemplate = paste(
                              '<b>%{text}</b><br>',
                              'PC1: %{x:.2f}<br>',
                              'PC2: %{y:.2f}<br>',
                              '<extra></extra>'
                            )) %>%
    layout(
      title = list(text = "PCA of Comparisons Based on Pathway NES", 
                   font = list(size = 16, face = "bold")),
      xaxis = list(title = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
                   zeroline = TRUE, zerolinewidth = 1, zerolinecolor = 'lightgray'),
      yaxis = list(title = paste0("PC2 (", round(var_explained[2], 1), "% variance)"),
                   zeroline = TRUE, zerolinewidth = 1, zerolinecolor = 'lightgray'),
      hovermode = 'closest',
      plot_bgcolor = 'white',
      paper_bgcolor = 'white'
    )
  
  # PCA loadings - top contributing pathways
  loadings <- as.data.frame(pca_result$rotation[, 1:2])
  loadings$Pathway <- rownames(loadings)
  loadings$Contribution <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  
  # Get top 20 pathways by contribution
  top_pathways <- loadings[order(-loadings$Contribution), ][1:min(20, nrow(loadings)), ]
  
  # Loadings plot
  p2 <- ggplot(top_pathways, aes(x = PC1, y = PC2)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkred", alpha = 0.7) +
    geom_text(aes(label = Pathway), size = 2.5, hjust = -0.1, vjust = 0) +
    labs(
      title = "Top 20 Pathways Driving PC1 and PC2",
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12)
    )
  
  ggsave(file.path(output_dir, "pca_pathway_loadings.pdf"), p2, width = 12, height = 10)
  ggsave(file.path(output_dir, "pca_pathway_loadings.png"), p2, width = 12, height = 10, dpi = 300)
  
  # Interactive loadings plot using plotly
  p2_interactive <- plot_ly(data = top_pathways,
                            x = ~PC1,
                            y = ~PC2,
                            type = 'scatter',
                            mode = 'markers+text',
                            text = ~Pathway,
                            textposition = 'top right',
                            marker = list(size = 8, color = 'darkred', opacity = 0.7),
                            hovertemplate = paste(
                              '<b>%{text}</b><br>',
                              'PC1 loading: %{x:.3f}<br>',
                              'PC2 loading: %{y:.3f}<br>',
                              'Contribution: %{customdata:.3f}<br>',
                              '<extra></extra>'
                            ),
                            customdata = ~Contribution) %>%
    layout(
      title = list(text = "Top 20 Pathways Driving PC1 and PC2",
                   font = list(size = 16, face = "bold")),
      xaxis = list(title = paste0("PC1 (", round(var_explained[1], 1), "%)"),
                   zeroline = TRUE, zerolinewidth = 1, zerolinecolor = 'lightgray'),
      yaxis = list(title = paste0("PC2 (", round(var_explained[2], 1), "%)"),
                   zeroline = TRUE, zerolinewidth = 1, zerolinecolor = 'lightgray'),
      hovermode = 'closest',
      plot_bgcolor = 'white',
      paper_bgcolor = 'white'
    )
  
  # Scree plot - static version
  scree_data <- data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    Variance = var_explained
  )
  scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)
  
  p3_static <- ggplot(scree_data, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = paste0(round(Variance, 1), "%")), 
              vjust = -0.5, size = 3) +
    labs(
      title = "PCA Scree Plot - Variance Explained",
      x = "Principal Component",
      y = "% Variance Explained"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(output_dir, "pca_scree_plot.pdf"), p3_static, width = 10, height = 6)
  ggsave(file.path(output_dir, "pca_scree_plot.png"), p3_static, width = 10, height = 6, dpi = 300)
  
  # Interactive scree plot using plotly
  p3_interactive <- plot_ly(data = scree_data,
                            x = ~PC,
                            y = ~Variance,
                            type = 'bar',
                            marker = list(color = 'steelblue'),
                            text = ~paste0(round(Variance, 1), "%"),
                            textposition = 'outside',
                            hovertemplate = paste(
                              '<b>%{x}</b><br>',
                              'Variance: %{y:.1f}%<br>',
                              '<extra></extra>'
                            )) %>%
    layout(
      title = list(text = "PCA Scree Plot - Variance Explained",
                   font = list(size = 16, face = "bold")),
      xaxis = list(title = "Principal Component"),
      yaxis = list(title = "% Variance Explained"),
      plot_bgcolor = 'white',
      paper_bgcolor = 'white'
    )
  
  # Save PCA results
  write.table(pca_scores, 
              file.path(output_dir, "pca_comparison_scores.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(loadings[order(-loadings$Contribution), ], 
              file.path(output_dir, "pca_pathway_loadings.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save variance explained
  var_summary <- data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    VarianceExplained = var_explained,
    CumulativeVariance = cumsum(var_explained)
  )
  write.table(var_summary,
              file.path(output_dir, "pca_variance_explained.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("PCA analysis complete")
  message(paste("  - Comparison biplot saved"))
  message(paste("  - Top pathway loadings saved"))
  message(paste("  - Scree plot saved"))
  message(paste("  - PC1 explains", round(var_explained[1], 1), "% of variance"))
  message(paste("  - PC2 explains", round(var_explained[2], 1), "% of variance"))
  
  return(list(
    pca_result = pca_result,
    scores = pca_scores,
    loadings = loadings,
    var_explained = var_explained,
    interactive_biplot = p1_interactive,
    interactive_loadings = p2_interactive,
    interactive_scree = p3_interactive
  ))
}

# Helper function to convert image to base64 data URI
image_to_base64 <- function(img_path) {
  if (file.exists(img_path)) {
    img_data <- base64enc::base64encode(img_path)
    return(paste0("data:image/png;base64,", img_data))
  } else {
    return("")
  }
}

# Function to generate interactive HTML report
generate_html_report <- function(results_list, file_names, output_dir, padj_cutoff, pathway_database, pca_results = NULL) {
  
  message("\nGenerating interactive HTML report...")
  
  tryCatch({
    # Prepare data
    comparison_names <- c()
    for (i in seq_along(results_list)) {
      if (!is.null(results_list[[i]])) {
        comparison_names <- c(comparison_names, tools::file_path_sans_ext(basename(file_names[i])))
      }
    }
    
    if (length(comparison_names) == 0) {
      message("No results to create HTML report")
      return(NULL)
    }
    
    # Create HTML content
    html_content <- HTML(paste0('
<!DOCTYPE html>
<html>
<head>
    <title>GSEA Analysis Report</title>
    <meta charset="UTF-8">
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 40px;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 8px;
        }
        h3 {
            color: #7f8c8d;
        }
        .summary-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .summary-table th {
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
        }
        .summary-table td {
            padding: 10px 12px;
            border-bottom: 1px solid #ecf0f1;
        }
        .summary-table tr:hover {
            background-color: #f8f9fa;
        }
        .tabs {
            overflow: hidden;
            background-color: #fff;
            border: 1px solid #ddd;
            margin-top: 20px;
        }
        .tab {
            background-color: inherit;
            float: left;
            border: none;
            outline: none;
            cursor: pointer;
            padding: 14px 16px;
            transition: 0.3s;
            font-size: 14px;
        }
        .tab:hover {
            background-color: #e8e8e8;
        }
        .tab.active {
            background-color: #3498db;
            color: white;
        }
        .tabcontent {
            display: none;
            padding: 20px;
            border: 1px solid #ddd;
            border-top: none;
            background-color: white;
        }
        .section {
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .plot-container {
            margin: 20px 0;
            text-align: center;
        }
        .plot-container img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 5px;
        }
        .info-box {
            background-color: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 15px 0;
        }
    </style>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <h1>', pathway_database, ' GSEA Analysis Report</h1>
    <div class="info-box">
        <strong>Analysis Parameters:</strong><br>
        Pathway Database: ', pathway_database, '<br>
        Adjusted p-value cutoff for plots: ', padj_cutoff, '<br>
        Number of comparisons: ', length(comparison_names), '
    </div>
'))
    
    # Section 1: Summary Table
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>Summary: Significant Pathways per Comparison</h2>
        <table class="summary-table">
            <tr>
                <th>Comparison</th>
                <th>Total Pathways (padj < ', padj_cutoff, ')</th>
                <th>Upregulated (NES > 0)</th>
                <th>Downregulated (NES < 0)</th>
                <th>Significant at FDR < 0.1</th>
            </tr>
'))
    
    # Fill summary table
    for (i in seq_along(results_list)) {
      if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
        comp_name <- tools::file_path_sans_ext(basename(file_names[i]))
        df <- results_list[[i]]@result
        
        sig_padj <- df[df$p.adjust < padj_cutoff, ]
        sig_fdr <- df[df$p.adjust < 0.1, ]
        up_count <- sum(sig_padj$NES > 0)
        down_count <- sum(sig_padj$NES < 0)
        
        html_content <- tagList(html_content, HTML(paste0('
            <tr>
                <td><strong>', comp_name, '</strong></td>
                <td>', nrow(sig_padj), '</td>
                <td style="color: #e74c3c;">', up_count, '</td>
                <td style="color: #3498db;">', down_count, '</td>
                <td>', nrow(sig_fdr), '</td>
            </tr>
')))
      }
    }
    
    html_content <- tagList(html_content, HTML('
        </table>
    </div>
'))
    
    # Section 2: Cross-Comparison Plots
    # Convert images to base64
    upset_up_base64 <- image_to_base64(file.path(output_dir, "upset_plot_upregulated_pathways.png"))
    upset_down_base64 <- image_to_base64(file.path(output_dir, "upset_plot_downregulated_pathways.png"))
    comp_up_base64 <- image_to_base64(file.path(output_dir, "comparative_dotplot_upregulated_top10.png"))
    comp_down_base64 <- image_to_base64(file.path(output_dir, "comparative_dotplot_downregulated_top10.png"))
    pca_biplot_base64 <- image_to_base64(file.path(output_dir, "pca_comparisons_biplot.png"))
    pca_loadings_base64 <- image_to_base64(file.path(output_dir, "pca_pathway_loadings.png"))
    pca_scree_base64 <- image_to_base64(file.path(output_dir, "pca_scree_plot.png"))
    
    html_content <- tagList(html_content, HTML(paste0('
    <div class="section">
        <h2>Cross-Comparison Analysis</h2>
        <h3>UpSet Plots - Pathway Overlap</h3>
        <div class="plot-container">
            <h4>Upregulated Pathways</h4>
            <img src="', upset_up_base64, '" alt="Upregulated Upset Plot">
        </div>
        <div class="plot-container">
            <h4>Downregulated Pathways</h4>
            <img src="', upset_down_base64, '" alt="Downregulated Upset Plot">
        </div>
        
        <h3>Comparative Dotplots</h3>
        <div class="plot-container">
            <h4>Top 10 Upregulated Pathways Across Comparisons</h4>
            <img src="', comp_up_base64, '" alt="Comparative Upregulated Dotplot">
        </div>
        <div class="plot-container">
            <h4>Top 10 Downregulated Pathways Across Comparisons</h4>
            <img src="', comp_down_base64, '" alt="Comparative Downregulated Dotplot">
        </div>
        
        <h3>Principal Component Analysis (PCA) - NES-based Clustering</h3>
        <p>PCA reveals which comparisons cluster together based on their pathway enrichment patterns (NES values). <strong>Interactive plots:</strong> Hover for details, zoom, and pan.</p>
')))
    
    # Add interactive PCA biplot if available
    if (!is.null(pca_results) && !is.null(pca_results$interactive_biplot)) {
      biplot_json <- plotly::plotly_json(pca_results$interactive_biplot, jsonedit = FALSE)
      
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Comparison Clustering (Interactive PCA Biplot)</h4>
            <div id="pca-biplot" style="width:100%;height:600px;"></div>
            <script>
              var biplotData = ', biplot_json, ';
              Plotly.newPlot("pca-biplot", biplotData.data, biplotData.layout);
            </script>
            <p style="color: #666; font-size: 0.9em;">Each point represents a comparison. Proximity indicates similar enrichment patterns. Hover for details.</p>
        </div>
')))
    } else {
      # Fallback to static image if PCA not available
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Comparison Clustering (PCA Biplot)</h4>
            <img src="', pca_biplot_base64, '" alt="PCA Biplot of Comparisons">
            <p style="color: #666; font-size: 0.9em;">Each point represents a comparison. Proximity indicates similar enrichment patterns.</p>
        </div>
')))
    }
    
    # Add interactive pathway loadings plot if available
    if (!is.null(pca_results) && !is.null(pca_results$interactive_loadings)) {
      loadings_json <- plotly::plotly_json(pca_results$interactive_loadings, jsonedit = FALSE)
      
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Top Pathways Driving Variation (Interactive)</h4>
            <div id="pca-loadings" style="width:100%;height:600px;"></div>
            <script>
              var loadingsData = ', loadings_json, ';
              Plotly.newPlot("pca-loadings", loadingsData.data, loadingsData.layout);
            </script>
            <p style="color: #666; font-size: 0.9em;">Each point shows a pathway that drives differences between comparisons. Hover for pathway details and contribution scores.</p>
        </div>
')))
    } else {
      # Fallback to static image if PCA not available
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Top Pathways Driving Variation</h4>
            <img src="', pca_loadings_base64, '" alt="PCA Pathway Loadings">
            <p style="color: #666; font-size: 0.9em;">Arrows show pathways that most strongly influence PC1 and PC2.</p>
        </div>
')))
    }
    
    # Add interactive scree plot if available
    if (!is.null(pca_results) && !is.null(pca_results$interactive_scree)) {
      scree_json <- plotly::plotly_json(pca_results$interactive_scree, jsonedit = FALSE)
      
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Variance Explained (Interactive Scree Plot)</h4>
            <div id="pca-scree" style="width:100%;height:500px;"></div>
            <script>
              var screeData = ', scree_json, ';
              Plotly.newPlot("pca-scree", screeData.data, screeData.layout);
            </script>
            <p style="color: #666; font-size: 0.9em;">Shows the % of variance explained by each principal component. Hover for details.</p>
        </div>
    </div>
')))
    } else {
      # Fallback to static image if PCA not available
      html_content <- tagList(html_content, HTML(paste0('
        <div class="plot-container">
            <h4>Variance Explained (Scree Plot)</h4>
            <img src="', pca_scree_base64, '" alt="PCA Scree Plot">
            <p style="color: #666; font-size: 0.9em;">Shows the % of variance explained by each principal component.</p>
        </div>
    </div>
')))
    }
    
    # Section 3: Interactive Volcano Plots
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>Interactive Volcano Plots</h2>
        <p>Click on points to see pathway details. All tested pathways are shown.</p>
        <div class="tabs" id="volcano-tabs">
'))
    
    for (i in seq_along(comparison_names)) {
      active_class <- if(i == 1) " active" else ""
      html_content <- tagList(html_content, HTML(paste0('
            <button class="tab', active_class, '" onclick="openTab(event, \'volcano-', i, '\', \'volcano-tabs\')">',
                                                         comparison_names[i], '</button>
')))
    }
    
    html_content <- tagList(html_content, HTML('
        </div>
'))
    
    # Generate volcano plots - embed directly as JSON
    for (i in seq_along(results_list)) {
      if (!is.null(results_list[[i]]) && nrow(results_list[[i]]@result) > 0) {
        comp_name <- comparison_names[i]
        df <- results_list[[i]]@result
        
        # Prepare volcano plot data
        df$log10padj <- -log10(df$p.adjust)
        df$significant <- ifelse(df$p.adjust < padj_cutoff, 
                                ifelse(df$NES > 0, "Upregulated", "Downregulated"), 
                                "Not Significant")
        
        # Create plotly volcano plot
        volcano_plot <- plot_ly(data = df, x = ~NES, y = ~log10padj,
                               type = 'scatter', mode = 'markers',
                               color = ~significant,
                               colors = c("Upregulated" = "#e74c3c", 
                                        "Downregulated" = "#3498db", 
                                        "Not Significant" = "gray"),
                               text = ~paste("Pathway:", Description,
                                           "<br>NES:", round(NES, 3),
                                           "<br>P.adjust:", format(p.adjust, scientific = TRUE)),
                               hoverinfo = 'text',
                               marker = list(size = 6, opacity = 0.7)) %>%
          layout(title = paste("Volcano Plot:", comp_name),
                 xaxis = list(title = "Normalized Enrichment Score (NES)"),
                 yaxis = list(title = "-log10(adjusted p-value)"),
                 hovermode = 'closest')
        
        # Convert plotly to JSON
        volcano_json <- plotly::plotly_json(volcano_plot, jsonedit = FALSE)
        
        display_style <- if(i == 1) "block" else "none"
        
        # Embed plot directly using plotly.js
        html_content <- tagList(html_content, 
                               HTML(paste0('<div id="volcano-', i, '" class="tabcontent" style="display:', display_style, ';">
                                          <div id="volcano-plot-', i, '" style="width:100%;height:600px;"></div>
                                          <script>
                                            var plotData = ', volcano_json, ';
                                            Plotly.newPlot("volcano-plot-', i, '", plotData.data, plotData.layout);
                                          </script>
                                      </div>')))
        
        message(paste("Volcano plot", i, "embedded directly in HTML"))
      }
    }
    
    html_content <- tagList(html_content, HTML('
    </div>
'))
    
    # Section 4: Individual GSEA Dotplots with tabs
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>GSEA Dotplots by Comparison</h2>
        <div class="tabs" id="dotplot-tabs">
'))
    
    for (i in seq_along(comparison_names)) {
      active_class <- if(i == 1) " active" else ""
      html_content <- tagList(html_content, HTML(paste0('
            <button class="tab', active_class, '" onclick="openTab(event, \'dotplot-', i, '\', \'dotplot-tabs\')">',
                                                         comparison_names[i], '</button>
')))
    }
    
    html_content <- tagList(html_content, HTML('
        </div>
'))
    
    for (i in seq_along(comparison_names)) {
      comp_name <- comparison_names[i]
      display_style <- if(i == 1) "block" else "none"
      
      # Convert image to base64
      dotplot_base64 <- image_to_base64(file.path(output_dir, comp_name, paste0(comp_name, "_gsea_dotplot.png")))
      
      html_content <- tagList(html_content, HTML(paste0('
        <div id="dotplot-', i, '" class="tabcontent" style="display:', display_style, ';">
            <div class="plot-container">
                <img src="', dotplot_base64, '" alt="', comp_name, ' Dotplot">
            </div>
        </div>
')))
    }
    
    html_content <- tagList(html_content, HTML('
    </div>
'))
    
    # Section 5: Ridge Plots with tabs
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>Ridge Plots by Comparison</h2>
        <div class="tabs" id="ridge-tabs">
'))
    
    for (i in seq_along(comparison_names)) {
      active_class <- if(i == 1) " active" else ""
      html_content <- tagList(html_content, HTML(paste0('
            <button class="tab', active_class, '" onclick="openTab(event, \'ridge-', i, '\', \'ridge-tabs\')">',
                                                         comparison_names[i], '</button>
')))
    }
    
    html_content <- tagList(html_content, HTML('
        </div>
'))
    
    for (i in seq_along(comparison_names)) {
      comp_name <- comparison_names[i]
      display_style <- if(i == 1) "block" else "none"
      
      # Convert image to base64
      ridge_base64 <- image_to_base64(file.path(output_dir, comp_name, paste0(comp_name, "_gsea_ridgeplot.png")))
      
      html_content <- tagList(html_content, HTML(paste0('
        <div id="ridge-', i, '" class="tabcontent" style="display:', display_style, ';">
            <div class="plot-container">
                <img src="', ridge_base64, '" alt="', comp_name, ' Ridge Plot">
            </div>
        </div>
')))
    }
    
    html_content <- tagList(html_content, HTML('
    </div>
'))
    
    # Section 6: Tree Plots with tabs
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>Tree Plots by Comparison</h2>
        <div class="tabs" id="tree-tabs">
'))
    
    for (i in seq_along(comparison_names)) {
      active_class <- if(i == 1) " active" else ""
      html_content <- tagList(html_content, HTML(paste0('
            <button class="tab', active_class, '" onclick="openTab(event, \'tree-', i, '\', \'tree-tabs\')">',
                                                         comparison_names[i], '</button>
')))
    }
    
    html_content <- tagList(html_content, HTML('
        </div>
'))
    
    for (i in seq_along(comparison_names)) {
      comp_name <- comparison_names[i]
      display_style <- if(i == 1) "block" else "none"
      
      # Convert image to base64
      tree_base64 <- image_to_base64(file.path(output_dir, comp_name, paste0(comp_name, "_gsea_treeplot.png")))
      
      html_content <- tagList(html_content, HTML(paste0('
        <div id="tree-', i, '" class="tabcontent" style="display:', display_style, ';">
            <div class="plot-container">
                <img src="', tree_base64, '" alt="', comp_name, ' Tree Plot">
            </div>
        </div>
')))
    }
    
    html_content <- tagList(html_content, HTML('
    </div>
'))
    
    # Section 7: Enrichment Maps with tabs
    html_content <- tagList(html_content, HTML('
    <div class="section">
        <h2>Enrichment Maps by Comparison</h2>
        <div class="tabs" id="emap-tabs">
'))
    
    for (i in seq_along(comparison_names)) {
      active_class <- if(i == 1) " active" else ""
      html_content <- tagList(html_content, HTML(paste0('
            <button class="tab', active_class, '" onclick="openTab(event, \'emap-', i, '\', \'emap-tabs\')">',
                                                         comparison_names[i], '</button>
')))
    }
    
    html_content <- tagList(html_content, HTML('
        </div>
'))
    
    for (i in seq_along(comparison_names)) {
      comp_name <- comparison_names[i]
      display_style <- if(i == 1) "block" else "none"
      
      # Convert image to base64
      emap_base64 <- image_to_base64(file.path(output_dir, comp_name, paste0(comp_name, "_gsea_emapplot.png")))
      
      html_content <- tagList(html_content, HTML(paste0('
        <div id="emap-', i, '" class="tabcontent" style="display:', display_style, ';">
            <div class="plot-container">
                <img src="', emap_base64, '" alt="', comp_name, ' Enrichment Map">
            </div>
        </div>
')))
    }
    
    html_content <- tagList(html_content, HTML('
    </div>
'))
    
    # Add JavaScript for tab functionality
    html_content <- tagList(html_content, HTML('
    <script>
    function openTab(evt, tabName, tabGroupId) {
        var i, tabcontent, tablinks;
        var tabGroup = document.getElementById(tabGroupId);
        
        // Get all elements with class="tabcontent" that are siblings
        tabcontent = tabGroup.parentElement.getElementsByClassName("tabcontent");
        for (i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = "none";
        }
        
        // Get all elements with class="tab" within this tab group
        tablinks = tabGroup.getElementsByClassName("tab");
        for (i = 0; i < tablinks.length; i++) {
            tablinks[i].className = tablinks[i].className.replace(" active", "");
        }
        
        document.getElementById(tabName).style.display = "block";
        evt.currentTarget.className += " active";
    }
    </script>
'))
    
    # Close HTML
    html_content <- tagList(html_content, HTML('
</body>
</html>
'))
    
    # Save HTML file
    html_file <- file.path(output_dir, "gsea_report.html")
    save_html(html_content, file = html_file)
    
    message(paste("Interactive HTML report saved to:", html_file))
    
  }, error = function(e) {
    message(paste("Could not create HTML report:", e$message))
  })
}

# Main execution
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("
Usage: Rscript preranked_gsea_analysis.R <input_directory> [options]

Required:
  input_directory       Path to directory containing differential expression TSV files

Options:
  --organism            'human' or 'mouse' (default: human)
  --rank-by             Column to rank by (default: stat)
                        Options: log2FoldChange, stat
  --pathway-database    Pathway database to use (default: GO)
                        Options: GO (Gene Ontology: Biological Processes), Reactome
  --gmt-file            Path to custom GMT file for gene sets (overrides --pathway-database)
                        GMT format: each line has gene_set_name, description, genes...
  --output-dir          Output directory (default: gsea_results)
  --padj-cutoff         Adjusted p-value cutoff for plotting (default: 0.25)
                        Determines which pathways appear in plots
  --min-gset-size       Minimum gene set size (default: 10)
  --max-gset-size       Maximum gene set size (default: 2000)

Examples:
  # Using GO for human
  Rscript preranked_gsea_analysis.R ./deg_tables --organism human --pathway-database GO
  
  # Using Reactome for mouse
  Rscript preranked_gsea_analysis.R ./deg_tables --organism mouse --pathway-database Reactome
  
  # Using custom GMT file
  Rscript preranked_gsea_analysis.R ./deg_tables --gmt-file my_gene_sets.gmt
  
Expected columns in input files:
  baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene, contrast

Output Structure:
  - Each input file gets its own subdirectory with individual results
  - A merged wide-format table with all results is created in the main output directory
")
    quit(status = 1)
  }
  
  # Set defaults
  input_dir <- args[1]
  organism <- "human"
  rank_by <- "stat"
  pathway_database <- "GO"
  output_dir <- "gsea_results"
  padj_cutoff <- 0.25
  min_gset_size <- 10
  max_gset_size <- 2000
  gmt_file <- NULL
  
  # Parse optional arguments
  if (length(args) > 1) {
    for (i in 2:length(args)) {
      if (args[i] == "--organism" && i < length(args)) {
        organism <- args[i + 1]
      } else if (args[i] == "--rank-by" && i < length(args)) {
        rank_by <- args[i + 1]
      } else if (args[i] == "--pathway-database" && i < length(args)) {
        pathway_database <- args[i + 1]
      } else if (args[i] == "--gmt-file" && i < length(args)) {
        gmt_file <- args[i + 1]
      } else if (args[i] == "--output-dir" && i < length(args)) {
        output_dir <- args[i + 1]
      } else if (args[i] == "--padj-cutoff" && i < length(args)) {
        padj_cutoff <- as.numeric(args[i + 1])
      } else if (args[i] == "--min-gset-size" && i < length(args)) {
        min_gset_size <- as.integer(args[i + 1])
      } else if (args[i] == "--max-gset-size" && i < length(args)) {
        max_gset_size <- as.integer(args[i + 1])
      }
    }
  }
  
  # If GMT file is provided, override pathway_database
  if (!is.null(gmt_file)) {
    if (!file.exists(gmt_file)) {
      stop(paste("GMT file not found:", gmt_file))
    }
    pathway_database <- gmt_file  # Use GMT file path as pathway_database
  }
  
  # Check if input directory exists
  if (!dir.exists(input_dir)) {
    stop(paste("Input directory not found:", input_dir))
  }
  
  # Find all TSV files in the directory
  tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  if (length(tsv_files) == 0) {
    stop(paste("No TSV files found in directory:", input_dir))
  }
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n=== Starting Batch GSEA Analysis ===\n")
  cat(paste("Input directory:", input_dir, "\n"))
  cat(paste("Files found:", length(tsv_files), "\n"))
  cat(paste("Organism:", organism, "\n"))
  if (grepl("\\.gmt$", pathway_database, ignore.case = TRUE)) {
    cat(paste("Pathway database: Custom GMT file -", basename(pathway_database), "\n"))
  } else {
    cat(paste("Pathway database:", pathway_database, "\n"))
  }
  cat(paste("Ranking metric:", rank_by, "\n"))
  cat(paste("Output directory:", output_dir, "\n\n"))
  
  # Store all results for merging
  all_gsea_results <- list()
  
  # Process each file
  for (i in seq_along(tsv_files)) {
    file <- tsv_files[i]
    base_name <- tools::file_path_sans_ext(basename(file))
    
    cat(paste0("\n[", i, "/", length(tsv_files), "] Processing: ", basename(file), "\n"))
    cat(paste(rep("=", 60), collapse = ""), "\n")
    
    # Create subdirectory for this file's results
    file_output_dir <- file.path(output_dir, base_name)
    
    # Run GSEA
    tryCatch({
      results <- perform_gsea(
        de_file = file,
        organism = organism,
        rank_by = rank_by,
        output_dir = file_output_dir,
        padj_cutoff = padj_cutoff,
        min_gset_size = min_gset_size,
        max_gset_size = max_gset_size,
        pathway_database = pathway_database
      )
      
      all_gsea_results[[i]] <- results
      
    }, error = function(e) {
      message(paste("ERROR processing", basename(file), ":", e$message))
      all_gsea_results[[i]] <- NULL
    })
  }
  
  # Merge all results into wide format
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  merged_file <- file.path(output_dir, "merged_gsea_results_wide.tsv")
  merged_results <- merge_gsea_results(all_gsea_results, tsv_files, merged_file)
  
  # Create comparative dotplots
  create_comparative_dotplots(all_gsea_results, tsv_files, output_dir)
  
  # Create upset plot for pathway overlap
  create_upset_plot(all_gsea_results, tsv_files, output_dir, padj_cutoff)
  
  # Perform PCA on NES values across comparisons
  pca_results <- create_nes_pca(all_gsea_results, tsv_files, output_dir, padj_cutoff)
  
  # Generate interactive HTML report
  generate_html_report(all_gsea_results, tsv_files, output_dir, padj_cutoff, pathway_database, pca_results)
  
  cat("\n=== Analysis Complete ===\n")
  cat(paste("Total files processed:", length(tsv_files), "\n"))
  cat(paste("Successful analyses:", sum(!sapply(all_gsea_results, is.null)), "\n"))
  cat(paste("Results location:", output_dir, "\n"))
  cat(paste("HTML report:", file.path(output_dir, "gsea_report.html"), "\n"))
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
