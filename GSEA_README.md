# Pre-ranked GSEA Analysis with ClusterProfiler - Batch Processing

## Installation

Before running the script, install the required R packages:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
    "clusterProfiler",
    "enrichplot",
    "org.Hs.eg.db",  # Human annotation
    "org.Mm.eg.db",  # Mouse annotation
    "DOSE",          # For similarity calculations (treeplot, emapplot)
    "ReactomePA"     # For Reactome pathway analysis
))

# Install CRAN packages
install.packages(c("ggplot2", "ComplexUpset", "plotly", "htmltools", "base64enc", "htmlwidgets"))
```

## Usage

The script now processes **all TSV files** in a directory at once and creates:
1. Individual results for each file in separate subdirectories
2. A merged wide-format table with results from all files

### Basic Usage (Human)
```bash
Rscript preranked_gsea_analysis.R /path/to/deg_tables/
```

### Mouse Analysis
```bash
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --organism mouse
```

### Using Reactome Pathways
```bash
# Human Reactome pathways
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --pathway-database Reactome

# Mouse Reactome pathways
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --organism mouse --pathway-database Reactome
```

### Using Custom GMT Files
You can use any gene set collection in GMT (Gene Matrix Transposed) format:

```bash
# Using custom GMT file
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --gmt-file my_gene_sets.gmt
```

**GMT File Format:**
Each line represents one gene set:
```
GENESET_NAME    DESCRIPTION    GENE1    GENE2    GENE3    ...
```

Example:
```
IMMUNE_RESPONSE    Immune response genes    CD4    CD8A    IL2    IFNG    TNF
CELL_CYCLE    Cell cycle genes    CDK1    CDK2    CCNA2    CCNB1    CCND1
```

**Where to get GMT files:**
- [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/) - Molecular Signatures Database
- [Enrichr](https://maayanlab.cloud/Enrichr/) - Download gene set libraries
- [KEGG](https://www.genome.jp/kegg/) - Download pathway gene sets
- Create your own custom gene sets

**Note:** When using GMT files:
- Gene symbols in the GMT file must match those in your differential expression data
- The `--organism` parameter is still required for certain operations but doesn't affect GMT file reading
- GMT file overrides the `--pathway-database` parameter

### Using log2FoldChange for ranking
```bash
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --organism human --rank-by log2FoldChange
```

### Custom Output Directory
```bash
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ --organism mouse --output-dir my_gsea_results
```

### All Options
```bash
# Using built-in pathway databases
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ \
  --organism human \
  --pathway-database GO \
  --rank-by stat \
  --output-dir gsea_results \
  --padj-cutoff 0.25 \
  --min-gset-size 10 \
  --max-gset-size 2000

# Using custom GMT file
Rscript preranked_gsea_analysis.R /path/to/deg_tables/ \
  --gmt-file my_gene_sets.gmt \
  --rank-by stat \
  --output-dir gsea_results \
  --padj-cutoff 0.25 \
  --min-gset-size 10 \
  --max-gset-size 2000
```

**Available Options:**
- `--organism`: Species database (`human` or `mouse`, default: `human`)
- `--pathway-database`: Built-in database (`GO` or `Reactome`, default: `GO`)
- `--gmt-file`: Path to custom GMT file (overrides `--pathway-database`)
- `--rank-by`: Column for ranking genes (`stat` or `log2FoldChange`, default: `stat`)
- `--output-dir`: Directory for results (default: `gsea_results`)
- `--padj-cutoff`: p-value cutoff for plots (default: `0.25`)
- `--min-gset-size`: Minimum genes per set (default: `10`)
- `--max-gset-size`: Maximum genes per set (default: `2000`)

## Input File Format

Your directory should contain TSV files with these columns:
- `baseMean`: Mean of normalized counts
- `log2FoldChange`: Log2 fold change
- `lfcSE`: Standard error of log2 fold change
- `stat`: Wald statistic
- `pvalue`: P-value
- `padj`: Adjusted p-value
- `gene`: Gene symbol (e.g., "TP53", "BRCA1")
- `contrast`: Comparison name

Example directory structure:
```
deg_tables/
├── treatment1_vs_control.tsv
├── treatment2_vs_control.tsv
└── treatment3_vs_control.tsv
```

## Ranking Metrics

You can rank genes by either:
- `stat`: Wald statistic (default, recommended for GSEA as it incorporates both effect size and significance)
- `log2FoldChange`: Magnitude of change

The `stat` column is preferred because it accounts for both the magnitude of change and the confidence in that change.

## Pathway Databases

You can choose between two pathway databases using the `--pathway-database` flag:

### Gene Ontology: Biological Processes (GO:BP) - Default
- **Pros**: 
  - Comprehensive coverage of biological processes
  - Well-established and widely used
  - Hierarchical structure allows for broad to specific analysis
  - Large number of gene sets (typically 5,000-10,000)
- **Cons**: 
  - Can have redundant terms
  - Some terms are very broad or non-specific
- **Use when**: You want comprehensive biological process annotation

### Reactome
- **Pros**:
  - Curated biological pathways with detailed mechanistic information
  - Less redundancy than GO
  - Focuses on specific molecular interactions and reactions
  - Easier to interpret results (fewer, more specific pathways)
- **Cons**:
  - Smaller number of pathways (typically 1,000-2,000)
  - Less comprehensive coverage than GO
- **Use when**: You want specific, mechanistic pathway information

**Recommendation**: Start with GO:BP for broad discovery, then use Reactome for more specific mechanistic insights.

**Note**: Reactome requires gene symbols to be converted to Entrez IDs internally. The script handles this automatically.

## Output Structure

The script creates the following output structure:

```
gsea_results/
├── treatment1_vs_control/
│   ├── treatment1_vs_control_gsea_all_results.tsv
│   ├── treatment1_vs_control_gsea_significant_FDR0.1.tsv
│   ├── treatment1_vs_control_gsea_dotplot.pdf
│   ├── treatment1_vs_control_gsea_enrichment_plot_upregulated.pdf
│   ├── treatment1_vs_control_gsea_enrichment_plot_downregulated.pdf
│   ├── treatment1_vs_control_gsea_ridgeplot.pdf
│   ├── treatment1_vs_control_gsea_treeplot.pdf
│   └── treatment1_vs_control_gsea_emapplot.pdf
├── treatment2_vs_control/
│   ├── treatment2_vs_control_gsea_all_results.tsv
│   ├── treatment2_vs_control_gsea_significant_FDR0.1.tsv
│   └── ...
├── treatment3_vs_control/
│   └── ...
├── merged_gsea_results_wide.tsv  ← Combined results from all files
├── comparative_dotplot_upregulated_top10.pdf  ← Cross-comparison plot
├── comparative_dotplot_downregulated_top10.pdf  ← Cross-comparison plot
├── upset_plot_upregulated_pathways.pdf  ← Pathway overlap visualization
├── upset_plot_downregulated_pathways.pdf  ← Pathway overlap visualization
├── pathway_overlap_summary.tsv  ← Overlap statistics
└── gsea_report.html  ← Interactive HTML report with all visualizations
```

### Individual File Outputs
For each input file, a subdirectory is created containing:
1. `*_gsea_all_results.tsv`: Complete GSEA results for all gene sets tested
2. `*_gsea_significant_FDR0.1.tsv`: Filtered results with only significant gene sets (FDR < 0.1)
3. `*_gsea_dotplot.pdf`: Dot plot of enriched pathways
4. `*_gsea_enrichment_plot_upregulated.pdf`: Enrichment plots for top 20 most significant upregulated pathways (positive NES) passing padj_cutoff
5. `*_gsea_enrichment_plot_downregulated.pdf`: Enrichment plots for top 20 most significant downregulated pathways (negative NES) passing padj_cutoff
6. `*_gsea_ridgeplot.pdf`: Ridge plot showing expression distributions
7. `*_gsea_treeplot.pdf`: Hierarchical clustering of GO terms based on semantic similarity
8. `*_gsea_emapplot.pdf`: Enrichment map showing network of enriched terms with edges representing gene overlap

### Cross-Comparison Outputs
In the main output directory:
1. `merged_gsea_results_wide.tsv`: Wide-format table with all results across all comparisons
2. `comparative_dotplot_upregulated_top10.pdf`: Dotplot showing top 10 upregulated pathways from each comparison and how they perform across ALL comparisons (NES = color, -log10(padj) = size)
3. `comparative_dotplot_downregulated_top10.pdf`: Dotplot showing top 10 downregulated pathways from each comparison and how they perform across ALL comparisons (NES = color, -log10(padj) = size)
4. `upset_plot_upregulated_pathways.pdf`: UpSet plot showing the overlap of significant upregulated pathways (at padj_cutoff) across all comparisons
5. `upset_plot_downregulated_pathways.pdf`: UpSet plot showing the overlap of significant downregulated pathways (at padj_cutoff) across all comparisons
6. `pathway_overlap_summary.tsv`: Summary statistics showing the number of upregulated and downregulated pathways per comparison
7. **`pca_comparisons_biplot.pdf`**: PCA biplot showing clustering of comparisons based on NES values
8. **`pca_pathway_loadings.pdf`**: Plot showing top pathways that drive variation between comparisons
9. **`pca_scree_plot.pdf`**: Scree plot showing variance explained by each principal component
10. `pca_comparison_scores.tsv`: PC coordinates for each comparison
11. `pca_pathway_loadings.tsv`: Pathway contributions to principal components
12. `pca_variance_explained.tsv`: Variance explained by each PC
13. **`gsea_report.html`**: Interactive HTML report with all visualizations organized in sections with tabs

### Merged Wide-Format Table
The `merged_gsea_results_wide.tsv` file contains:
- **Rows**: All unique gene sets tested across all input files
- **Columns**: 
  - `ID`: GO term ID
  - `Description`: GO term description
  - For each input file: `filename_NES`, `filename_pvalue`, `filename_padj`

Example columns:
```
ID          Description                    treatment1_NES  treatment1_pvalue  treatment1_padj  treatment2_NES  treatment2_pvalue  treatment2_padj
GO:0006915  apoptotic process             2.34            0.001              0.02             1.89            0.005              0.08
GO:0008283  cell proliferation            -1.67           0.010              0.15             -2.01           0.002              0.03
```

This format allows easy comparison of enrichment scores across all your contrasts/conditions.

### Comparative Dotplots
The script creates two cross-comparison dotplots:

**`comparative_dotplot_upregulated_top10.pdf`**: 
- Takes the top 10 upregulated pathways from each comparison (at FDR<0.1)
- Creates a union of all unique pathways across comparisons
- Shows how each pathway performs across ALL comparisons
- **X-axis**: Comparisons
- **Y-axis**: Pathways
- **Dot color**: NES (red = upregulated, blue = downregulated, white = neutral)
- **Dot size**: -log10(adjusted p-value) (larger = more significant)
- Missing dots indicate pathway not significant in that comparison

**`comparative_dotplot_downregulated_top10.pdf`**: Same as above but for downregulated pathways

These plots let you quickly identify:
- Pathways consistently enriched across multiple conditions
- Condition-specific pathway responses
- The strength and direction of enrichment across comparisons

### UpSet Plots
The script creates two separate UpSet plots:

**`upset_plot_upregulated_pathways.pdf`**:
- Shows ALL significant **upregulated** pathways (positive NES) at your `--padj-cutoff` threshold across comparisons
- Visualizes which upregulated pathways are shared between different combinations of comparisons
- Red/dark red color scheme to indicate upregulation

**`upset_plot_downregulated_pathways.pdf`**:
- Shows ALL significant **downregulated** pathways (negative NES) at your `--padj-cutoff` threshold across comparisons
- Visualizes which downregulated pathways are shared between different combinations of comparisons
- Blue color scheme to indicate downregulation

**Both plots include:**
- Bars on the left showing the total number of significant pathways per comparison
- Connected dots in the matrix showing which comparisons share pathways
- Bars on top showing the size of each intersection (e.g., pathways found in comparison A AND B but not C)
- Ordered by frequency to highlight the most common overlap patterns

This is particularly useful for identifying:
- Core pathways that are consistently upregulated or downregulated across all or most comparisons
- Unique pathways specific to individual comparisons
- Pairwise or multi-way overlaps between specific combinations of conditions
- Differences in the overlap patterns between upregulated vs downregulated pathways

### Principal Component Analysis (PCA) 🆕

The script now performs **PCA on NES values** across all comparisons to reveal global patterns:

**`pca_comparisons_biplot.pdf`**:
- Shows clustering of **comparisons** based on their pathway enrichment patterns
- Each point = one comparison
- **Proximity indicates similarity**: Comparisons close together have similar biological responses
- PC1 and PC2 axes show the main sources of variation

**`pca_pathway_loadings.pdf`**:
- Shows the **top 20 pathways** that drive the differences between comparisons
- Arrows indicate pathways with strongest contributions to PC1 and PC2
- Longer arrows = pathways that most differentiate your conditions
- Identifies which biological processes separate your experimental groups

**`pca_scree_plot.pdf`**:
- Shows % variance explained by each principal component
- Helps assess how many dimensions capture the biological variation
- PC1 typically captures the largest source of variation

**What PCA Reveals:**
- Which experimental conditions produce similar pathway enrichment patterns
- Which biological processes most strongly differentiate your comparisons
- Unexpected similarities or differences between conditions
- Main axes of biological variation in your dataset

**Data Files:**
- `pca_comparison_scores.tsv`: PC coordinates for each comparison (for custom plotting)
- `pca_pathway_loadings.tsv`: All pathway contributions to PCs
- `pca_variance_explained.tsv`: Variance by component

The PCA uses Normalized Enrichment Scores (NES) from **all tested pathways**, centered and scaled. This provides an unbiased view of the relationships between your comparisons.

**🎯 Interactive in HTML Report!**
The PCA biplot and scree plot are **fully interactive** in the HTML report with:
- Hover tooltips showing exact coordinates and values
- Zoom and pan capabilities
- Built-in download to PNG
- Dynamic exploration of your data

**See `PCA_FEATURE_DOCUMENTATION.md` for detailed interpretation guide.**
**See `INTERACTIVE_PCA_FEATURE.md` for interactive features guide.**

### Interactive HTML Report

**`gsea_report.html`**:

The script automatically generates a comprehensive, interactive HTML report that consolidates all analysis results in a single, easy-to-navigate document. The report includes:

1. **Summary Table**: Overview of significant pathways for each comparison, split by upregulated/downregulated
2. **Cross-Comparison Section**: 
   - UpSet plots showing pathway overlaps
   - Comparative dotplots
   - **Interactive PCA plots** 🆕 - biplot and scree plot with hover tooltips, zoom, and pan
3. **Interactive Volcano Plots** (tabbed by comparison): 
   - Shows ALL tested pathways (not just significant ones)
   - Interactive - hover to see pathway details
   - Click and zoom functionality
   - Color-coded by significance and direction
4. **Individual GSEA Plots** (each type in separate tabbed sections):
   - Dotplots
   - Ridge plots
   - Tree plots
   - Enrichment maps

**To view**: Simply open `gsea_report.html` in any web browser. The report is **completely self-contained** with all images and interactive plots embedded. No external files needed!

**Benefits**:
- Single self-contained HTML file for easy sharing
- **Interactive PCA, volcano, and other plots** for data exploration
- Hover tooltips show exact values and details
- Zoom and pan capabilities on interactive plots
- Tabbed interface makes it easy to compare across conditions
- Professional appearance suitable for reports and presentations
- No need to open multiple PDF files

## Visualization Types

The script generates multiple complementary visualizations:

**Individual Comparison Visualizations:**

1. **Interactive Volcano Plots** (in HTML report): Shows ALL tested pathways (not just significant ones) with NES on x-axis and -log10(adjusted p-value) on y-axis. Hover over points to see pathway details. Great for exploring borderline pathways and getting an overview of enrichment patterns.

2. **Dotplot**: Shows enriched pathways with dot size representing gene count and color representing adjusted p-value. Good for quick overview of top pathways.

3. **Enrichment plots (upregulated/downregulated)**: Classic GSEA running enrichment score plots showing where genes from each pathway are located in your ranked list. Split by enrichment direction for clarity. Shows the top 20 most significant pathways (or all if fewer than 20) for each direction.

4. **Ridge plot**: Shows the distribution of genes within each pathway across your ranked list. Useful for understanding the expression patterns within pathways.

5. **Tree plot**: Hierarchical clustering of GO terms based on semantic similarity. Helps identify groups of related pathways and reduces redundancy in interpretation.

6. **Enrichment map (emapplot)**: Network visualization where nodes are pathways and edges connect pathways that share genes. The thickness of edges represents the degree of gene overlap. Excellent for understanding relationships between enriched biological processes.

**Cross-Comparison Visualizations:**

7. **Comparative upregulated dotplot**: Shows the union of top 10 upregulated pathways from all comparisons. Each pathway is shown across all comparisons with dot color indicating NES and dot size indicating significance. Perfect for identifying consistently upregulated pathways or condition-specific responses.

8. **Comparative downregulated dotplot**: Same as above but for downregulated pathways. These plots are ideal for manuscript figures showing pathway enrichment patterns across multiple experimental conditions.

9. **UpSet plots (upregulated and downregulated)**: Two separate plots showing the overlap of significant pathways (at your `--padj-cutoff` threshold) across all comparisons. Unlike Venn diagrams, UpSet plots can effectively show overlaps among many sets. Each plot displays:
   - Bar chart showing the number of pathways per comparison
   - Matrix showing which comparisons share pathways
   - Bar chart showing the size of each intersection
   Separating by direction allows you to see if upregulated and downregulated pathways have different overlap patterns across conditions.

## Tips

1. **Gene symbols**: Ensure your gene column contains valid gene symbols (not Ensembl IDs). If you have Ensembl IDs, you'll need to convert them first.

2. **Ranking metric**: The `stat` column is generally better for GSEA than `log2FoldChange` because it incorporates statistical significance.

3. **Two result tables**: The script outputs both all results and filtered significant results. The "all results" table is useful for exploring borderline results or applying custom cutoffs, while the "significant" table (FDR < 0.1) is ready for direct reporting. The `--padj-cutoff` parameter controls which pathways appear in the visualization plots.

4. **Pathway database choice**: GO:BP provides comprehensive coverage and is great for discovery, while Reactome offers more specific mechanistic pathways. Consider running both and comparing results. GO typically yields more significant pathways but with more redundancy.

5. **Merged table**: The merged wide-format table is perfect for creating heatmaps or comparing enrichment patterns across multiple conditions. You can easily filter this table to find pathways enriched in multiple conditions.

6. **Enrichment plots**: The script creates two separate enrichment plots - one for upregulated pathways (positive NES) and one for downregulated pathways (negative NES). Each plot shows the top 20 most significant pathways (sorted by adjusted p-value) that pass the `--padj-cutoff` threshold, or all pathways if there are fewer than 20. This keeps the plots readable while focusing on the most important findings.

7. **Network and clustering plots**: The treeplot shows hierarchical clustering of GO terms based on semantic similarity, helping identify groups of related pathways. The enrichment map (emapplot) displays pathways as a network where nodes are GO terms and edges represent gene overlap - connected pathways share genes and may represent related biological processes.

8. **Comparative dotplots**: After processing all files, the script creates two powerful visualization showing the top pathways across ALL your comparisons. These plots make it easy to see which pathways are consistently enriched (appearing in multiple comparisons) versus condition-specific, and the relative strength of enrichment across conditions.

9. **UpSet plots for pathway overlap**: The upset plots use your `--padj-cutoff` value to determine which pathways are "significant" for overlap analysis. Two separate plots are created - one for upregulated pathways and one for downregulated pathways. This separation allows you to identify whether upregulated and downregulated pathways show different patterns of overlap across your conditions. The plots are ordered by intersection frequency to highlight the most common patterns.

10. **Gene set size**: Default range is 10-2000 genes per set. Adjust with `--min-gset-size` and `--max-gset-size` if needed.

11. **Interactive HTML report**: The HTML report is generated automatically and provides a convenient way to explore all results in one place. The interactive volcano plots are particularly useful for identifying borderline pathways that didn't meet the significance threshold but might be biologically interesting. Simply open the HTML file in any modern web browser.

12. **Progress tracking**: The script shows progress for each file being processed, making it easy to monitor batch analyses.

## Example Analysis Workflow

```bash
# 1. Organize your DEG tables in a directory
mkdir deg_tables
cp treatment*_vs_control.tsv deg_tables/

# 2. Run batch GSEA analysis
Rscript preranked_gsea_analysis.R deg_tables/ --organism human

# 3. Review results
# - Open gsea_report.html in your web browser for interactive exploration
# - Check individual subdirectories for detailed results per comparison
# - Open merged_gsea_results_wide.tsv for cross-comparison analysis
# - Use the merged table to identify common pathways across conditions
# - Share the HTML report with collaborators or include in presentations
```

## Troubleshooting

**Problem**: "Gene symbols not recognized"
- Solution: Make sure you're using the correct organism and that gene symbols match the annotation database (e.g., human uses "TP53", mouse uses "Trp53")

**Problem**: No significant results for some files
- Solution: Check the individual file output directories. Some comparisons may have fewer significant pathways than others, which is normal.

**Problem**: Missing packages
- Solution: Install all required packages as shown in the Installation section above

**Problem**: Script skips some files
- Solution: Check that all files are in TSV format with proper column headers. Error messages will indicate which files failed and why.

**Problem**: Reactome returns fewer pathways than expected
- Solution: Reactome requires Entrez IDs, so gene symbols that can't be converted are dropped. Check the console message showing how many genes were successfully converted. If the conversion rate is low, verify your gene symbols are current and match the organism database.

**Problem**: Images not showing in HTML report
- Solution: Make sure all PNG files have been generated successfully. The HTML report references PNG versions of plots that should be created alongside the PDF versions. Check that PNG files exist in the output directory and subdirectories. If only PDFs exist, there may have been an error during PNG generation.

**Problem**: Upset plots appear blank
- Solution: Check the console output for debug messages showing how many pathways are being plotted. If you see "Need at least 2 comparisons with upregulated/downregulated pathways", there aren't enough overlapping pathways. The upset plot requires at least 2 comparisons with some shared significant pathways. If the debug output shows data but plots are still blank, there may be an issue with the data structure being passed to UpSetR.
