# rna-seq_toolkit
Collection of scripts for downstream RNAseq analysis

DEG Analysis Guide

# Enhanced DESeq2 Script with Difference-of-Differences (Interaction) Contrasts

## What's New?

This enhanced version of your DESeq2 script now supports **difference-of-differences** (also called **interaction contrasts**), allowing you to test questions like:

- "Is the treatment effect different between disease and healthy samples?"
- "Does the response to stimulus A differ between two conditions?"
- "Is the fold change in Group A vs B different from the fold change in Group C vs D?"

## How It Works

### Simple Contrasts (Original Functionality)

Your original simple contrasts still work exactly as before:

```
treat                   ctrl
Healthy_ADP-h_1hr      Healthy_un-stimulate
ROSAH_ADP-h_1hr        ROSAH_un-stimulate
```

This compares two groups directly and asks: "What genes are differentially expressed between treat and ctrl?"

### Difference-of-Differences Contrasts (NEW!)

The new interaction contrasts follow this format:

```
treat                                                                    ctrl
(GroupA-GroupB)-(GroupC-GroupD)                                          .
```

For example:
```
treat                                                                                ctrl
(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)      .
```

This asks: **"Is the effect of ADP-h stimulation different between Healthy and ROSAH samples?"**

In other words, it tests whether the genes that respond to ADP-h treatment in Healthy samples respond *differently* in ROSAH samples.

## Biological Interpretation

### Example 1: Treatment effect modification by disease status

```
(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)
```

**What it tests:** Genes that show a *different* response to ADP-h treatment in Healthy vs ROSAH samples.

**Biological question:** Which genes respond differently to the stimulus depending on disease status?

**Interpretation of results:**
- **Positive log2FC:** Gene is more upregulated (or less downregulated) by ADP-h in Healthy compared to ROSAH
- **Negative log2FC:** Gene is more upregulated (or less downregulated) by ADP-h in ROSAH compared to Healthy

### Example 2: Drug effect comparison

```
(Healthy_ADP-h+diABZi_1hr-Healthy_ADP-h_1hr)-(ROSAH_ADP-h+diABZi_1hr-ROSAH_ADP-h_1hr)
```

**What it tests:** Whether the effect of adding diABZi to ADP-h differs between Healthy and ROSAH.

**Biological question:** Does diABZi modulate the response differently in the two cell types?

## Usage

### 1. Create Your Contrasts File

You can mix simple and interaction contrasts in the same file:

**contrasts.txt:**
```
treat                                                                                ctrl
Healthy_ADP-h_1hr                                                                   Healthy_un-stimulate
ROSAH_ADP-h_1hr                                                                     ROSAH_un-stimulate
(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)      .
(Healthy_ADP-h_4hr-Healthy_un-stimulate)-(ROSAH_ADP-h_4hr-ROSAH_un-stimulate)      .
```

**Note:** For interaction contrasts, put a `.` or leave the `ctrl` column empty.

### 2. Run the Script

Same as before:

```bash
Rscript deseq2_cpm_html3_FIXED_v4_with_interactions.R \
  --counts counts.tsv \
  --design all_design.txt \
  --contrasts contrasts_with_interactions.txt \
  --outdir results \
  --sample-col SampleID \
  --min-cpm 1 \
  --min-samples 2 \
  --sig-alpha 0.05
```

### 3. Outputs

The script produces the same outputs as before, but now includes interaction contrasts:

- **Volcano plots:** One per contrast (including interactions)
- **TSV files:** Per-contrast results in `contrast_tables/`
- **Combined CSV:** All contrasts in wide and long format
- **HTML report:** Interactive visualizations

## Technical Details

### Model Design

For **simple contrasts** (GroupA vs GroupB):
- Model: `~ Patient_ID + Group` (if Patient_ID blocking is feasible)
- Falls back to `~ Group` if patient blocking isn't possible

For **interaction contrasts** `(A-B)-(C-D)`:
- Model: `~ Patient_ID + Group` (includes all 4 groups)
- Uses a custom contrast vector: `coef(A) - coef(B) - coef(C) + coef(D)`
- Reference level is set to GroupD
- The contrast vector represents: `(GroupA - GroupD) - (GroupB - GroupD) - (GroupC - GroupD)`
  Which simplifies to: `GroupA - GroupB - GroupC` (in coefficient space)

### Sample Filtering

For interaction contrasts, CPM filtering is applied across **all 4 groups** together, which is more stringent than individual comparisons. Genes must be expressed in at least `min_samples` across the combined sample set.

### Patient Blocking

The script attempts to use Patient_ID blocking for interaction contrasts if:
1. Patient_ID column exists in the design
2. There are at least 2 patients in the combined 4-group dataset

The blocking doesn't require patients to span all 4 groups (which may not be feasible in your design), but it does provide better control for patient-level effects.

## Contrast Format Rules

**Valid formats:**

✅ `(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)`

**Important:**
- Group names cannot contain parentheses `()` themselves
- Use exactly this pattern: `(A-B)-(C-D)` with hyphens between groups
- Spaces around hyphens are OK: `(A - B) - (C - D)` also works

## Real-World Example Using Your Data

Based on your design file, here are some biologically meaningful interaction contrasts:

```
# Does ADP-h treatment effect differ by disease status?
(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)
(Healthy_ADP-h_4hr-Healthy_un-stimulate)-(ROSAH_ADP-h_4hr-ROSAH_un-stimulate)
(Healthy_ADP-h_24hr-Healthy_un-stimulate)-(ROSAH_ADP-h_24hr-ROSAH_un-stimulate)

# Does diABZi treatment effect differ by disease status?
(Healthy_diABZi_1hr-Healthy_un-stimulate)-(ROSAH_diABZi_1hr-ROSAH_un-stimulate)
(Healthy_diABZi_4hr-Healthy_un-stimulate)-(ROSAH_diABZi_4hr-ROSAH_un-stimulate)
(Healthy_diABZi_24hr-Healthy_un-stimulate)-(ROSAH_diABZi_24hr-ROSAH_un-stimulate)

# Does the effect of adding diABZi to ADP-h differ by disease status?
(Healthy_ADP-h+diABZi_1hr-Healthy_ADP-h_1hr)-(ROSAH_ADP-h+diABZi_1hr-ROSAH_ADP-h_1hr)

# Does the time effect (4hr vs 1hr) of ADP-h differ by disease status?
(Healthy_ADP-h_4hr-Healthy_ADP-h_1hr)-(ROSAH_ADP-h_4hr-ROSAH_ADP-h_1hr)
```

## Interpreting Results

For an interaction contrast `(A-B)-(C-D)`:

**Positive log2FoldChange:**
- The difference A-B is **larger** than the difference C-D
- The effect in the first pair (A-B) is stronger

**Negative log2FoldChange:**
- The difference C-D is **larger** than the difference A-B  
- The effect in the second pair (C-D) is stronger

**Example:**
For `(Healthy_ADP-h_1hr-Healthy_un-stimulate)-(ROSAH_ADP-h_1hr-ROSAH_un-stimulate)`:

- Gene with log2FC = +2, padj < 0.05:
  - This gene is upregulated by ADP-h treatment MORE in Healthy than in ROSAH
  - Or equivalently: the treatment effect is stronger in Healthy samples

- Gene with log2FC = -2, padj < 0.05:
  - This gene is upregulated by ADP-h treatment MORE in ROSAH than in Healthy
  - Or equivalently: the treatment effect is stronger in ROSAH samples

## Backwards Compatibility

All your existing contrasts files will work without modification. The script automatically detects whether a contrast is simple or interaction-based by looking for the `(A-B)-(C-D)` pattern.

## Troubleshooting

**Error: "Invalid interaction contrast format"**
- Check that you're using the exact format: `(GroupA-GroupB)-(GroupC-GroupD)`
- Ensure group names match exactly what's in your design file's `Group` column

**Error: "one or more groups empty"**
- One or more of the 4 groups in your interaction contrast has no samples
- Double-check group names against your design file

**Warning: "cannot find all required coefficients in model"**
- The model couldn't estimate coefficients for all groups
- This can happen if groups are confounded with other variables
- Try a simpler model or check your design

## Questions?

The key difference between simple and interaction contrasts:

- **Simple contrast:** "Is gene X different between A and B?"
- **Interaction contrast:** "Is the difference between A and B *different from* the difference between C and D?"

Interaction contrasts are powerful for finding genes that respond differently to treatments depending on condition, genotype, disease status, etc.


GSEA Guide

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


Timecourse analysis guide

# Batch Processing Guide - Multiple Timecourse Analyses

## Overview

The `batch_timecourse.sh` wrapper script allows you to run **multiple timecourse analyses in a single command**. Perfect for analyzing multiple treatment groups or testing multiple group comparisons!

## What It Does

### Option 1: Multiple Single-Group Analyses
Run separate timecourse analyses for each treatment group:
- Each group gets its own analysis
- Each group gets its own output directory
- All use the same parameters (method, clustering, enrichment, etc.)

### Option 2: Multiple Interaction Analyses  
Run separate interaction tests for multiple group pairs:
- Each pair gets its own interaction analysis
- Each comparison gets its own output directory
- All use the same parameters

## Quick Start

### Analyze Multiple Groups
```bash
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --groups wt,ASC_KO,BST1_KO,GSDMD_KO \
  --output batch_results
```

### Analyze Multiple Interactions
```bash
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --interactions "wt-ASC_KO,wt-BST1_KO,wt-GSDMD_KO" \
  --output batch_results
```

## Full Usage

```bash
./batch_timecourse.sh [OPTIONS]

Required:
  --counts FILE               Path to counts file
  --design FILE               Path to design file
  
Batch Options (choose one):
  --groups GROUP1,GROUP2,...  Run separate timecourse for each group
  --interactions PAIR1,PAIR2  Run interaction analyses for pairs (format: wt-ASC_KO)
  
Optional:
  --script PATH               Path to timecourse_analysis_cli.R
  --time-column NAME          Time column name [default: Time]
  --group-column NAME         Group column name [default: Treatment]
  --output DIR                Base output directory [default: batch_results]
  --method STRING             Analysis method: limma, deseq2, or both [default: both]
  --spline-df INT             Spline degrees of freedom [default: 3]
  --min-count INT             Minimum count threshold [default: 10]
  --min-samples INT           Minimum samples threshold [default: 3]
  --cluster-genes             Enable gene clustering
  --n-clusters INT            Number of clusters [default: 6]
  --run-enrichment            Enable pathway enrichment
  --enrichment-db STRING      Enrichment database [default: reactome]
  --organism STRING           Organism: human or mouse [default: human]
  --verbose                   Verbose output
```

## Output Structure

### For Multiple Groups
```
batch_results/
├── wt/
│   ├── limma/
│   ├── deseq2/
│   └── figures/
├── ASC_KO/
│   ├── limma/
│   ├── deseq2/
│   └── figures/
├── BST1_KO/
│   ├── limma/
│   ├── deseq2/
│   └── figures/
└── GSDMD_KO/
    ├── limma/
    ├── deseq2/
    └── figures/
```

### For Multiple Interactions
```
batch_results/
├── wt_vs_ASC_KO/
│   ├── limma/
│   ├── deseq2/
│   └── figures/
├── wt_vs_BST1_KO/
│   ├── limma/
│   ├── deseq2/
│   └── figures/
└── wt_vs_GSDMD_KO/
    ├── limma/
    ├── deseq2/
    └── figures/
```

## Your Data Examples

### Analyze All KO Lines
```bash
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --groups wt,ASC_KO,BST1_KO,GSDMD_KO \
  --method deseq2 \
  --cluster-genes \
  --n-clusters 6 \
  --run-enrichment \
  --enrichment-db reactome \
  --organism human \
  --output all_KO_lines
```

### Compare Each KO to WT
```bash
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --interactions "wt-ASC_KO,wt-BST1_KO,wt-GSDMD_KO" \
  --method deseq2 \
  --cluster-genes \
  --n-clusters 6 \
  --run-enrichment \
  --enrichment-db reactome \
  --organism human \
  --output KO_vs_WT_comparisons
```

### With and Without Nigericin
```bash
# Individual analyses
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --groups wt,wt_Nigericin,ASC_KO,ASC_KO_Nigericin \
  --output Nigericin_effects

# Interaction analyses
./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --interactions "wt-wt_Nigericin,ASC_KO-ASC_KO_Nigericin" \
  --output Nigericin_interactions
```

## How It Works

### For Single Groups
1. Reads your design file
2. For each group specified:
   - Filters samples to that group only
   - Creates group-specific design file
   - Runs complete timecourse analysis
   - Saves to group-specific directory
3. Reports success/failure for each

### For Interactions
1. For each pair specified:
   - Filters samples to those two groups
   - Runs interaction analysis (time × group)
   - Saves to comparison-specific directory
2. Reports success/failure for each

## Benefits

### 1. Time Savings
Run all analyses with one command instead of dozens:
```bash
# Old way (manual)
Rscript timecourse_analysis_cli.R ... --output wt
Rscript timecourse_analysis_cli.R ... --output ASC_KO
Rscript timecourse_analysis_cli.R ... --output BST1_KO
Rscript timecourse_analysis_cli.R ... --output GSDMD_KO

# New way (batch)
./batch_timecourse.sh --groups wt,ASC_KO,BST1_KO,GSDMD_KO ...
```

### 2. Consistency
All analyses use identical parameters:
- Same filtering thresholds
- Same statistical method
- Same clustering settings
- Same enrichment database
- Eliminates parameter inconsistencies

### 3. Organization
Automatic directory structure:
- Clear naming (group names or comparisons)
- Easy to find results
- Ready for publication

### 4. Error Handling
- Continues if one analysis fails
- Reports which succeeded/failed
- Logs for troubleshooting

## Advanced Usage

### Custom Parameters Per Analysis
If you need different parameters for different groups, run the script multiple times:

```bash
# WT with 8 clusters
./batch_timecourse.sh \
  --counts counts.txt --design design.txt \
  --groups wt \
  --n-clusters 8 \
  --output results

# KOs with 6 clusters  
./batch_timecourse.sh \
  --counts counts.txt --design design.txt \
  --groups ASC_KO,BST1_KO,GSDMD_KO \
  --n-clusters 6 \
  --output results
```

### Parallel Processing
Run analyses in parallel (if you have multiple cores):

```bash
# Run each group in background
for group in wt ASC_KO BST1_KO GSDMD_KO; do
  Rscript timecourse_analysis_cli.R \
    --counts counts.txt \
    --design design.txt \
    --group-column Treatment \
    --output results_${group} \
    ... &
done
wait  # Wait for all to complete
```

### Subset by Other Criteria
If your design has additional columns (e.g., Batch, Replicate):

```bash
# Analyze only specific batches
grep "Batch1" design.txt > design_batch1.txt

./batch_timecourse.sh \
  --counts counts.txt \
  --design design_batch1.txt \
  --groups wt,ASC_KO,BST1_KO \
  --output batch1_results
```

## Troubleshooting

### "Group not found"
Check group names match exactly (case-sensitive):
```bash
# View unique groups in your design
cut -f3 timecourse_design.txt | sort | uniq
```

### "Not enough samples"
Some groups may have too few samples after filtering:
- Check sample counts per group
- May need to relax `--min-samples` threshold

### "Script not found"
Specify full path to R script:
```bash
./batch_timecourse.sh \
  --script /path/to/timecourse_analysis_cli.R \
  ...
```

### One Analysis Fails
The batch script continues even if one fails:
- Check console output for error messages
- Failed analyses are counted separately
- Successful ones are still saved

### Memory Issues
If analyzing many large groups, may run out of memory:
- Run fewer groups at a time
- Use `--method limma` (less memory than DESeq2)
- Increase system memory

## Console Output

Example output when running:

```
================================================================================
BATCH TIMECOURSE ANALYSIS
================================================================================

Configuration:
  Counts file: RawCountFile_rsemgenes_reformat.txt
  Design file: timecourse_design.txt
  Time column: Time
  Group column: Treatment
  Method: deseq2
  Output directory: batch_results

Will analyze 4 groups:
  - wt
  - ASC_KO
  - BST1_KO
  - GSDMD_KO

================================================================================
Analyzing group: wt
================================================================================

  Samples in wt: 12
  Timepoints: 0, 1, 4, 18

[... analysis output ...]

✓ Completed: wt

================================================================================
Analyzing group: ASC_KO
================================================================================

  Samples in ASC_KO: 12
  Timepoints: 0, 1, 4, 18

[... analysis output ...]

✓ Completed: ASC_KO

[... continues for all groups ...]

================================================================================
BATCH ANALYSIS COMPLETE
================================================================================
  Total groups: 4
  Completed: 4
  Failed: 0
  Time elapsed: 1847s

Results saved to: batch_results
```

## Comparison Table

| Feature | Manual | Batch Script |
|---------|--------|--------------|
| **Commands needed** | One per analysis | One total |
| **Parameter consistency** | Manual checking | Automatic |
| **Output organization** | Manual naming | Automatic |
| **Error handling** | Stops on error | Continues |
| **Time tracking** | Manual | Automatic |
| **Parallelization** | Complex | Simple |

## Tips & Best Practices

### 1. Test First
Run one group manually to confirm parameters:
```bash
Rscript timecourse_analysis_cli.R \
  --counts counts.txt \
  --design design.txt \
  --method deseq2 \
  --output test_run
  
# If successful, use same parameters in batch
```

### 2. Use Descriptive Output Names
```bash
--output KO_timecourse_Dec2025
--output WT_vs_KO_interactions_Dec2025
```

### 3. Save Your Command
Create a script file:
```bash
#!/bin/bash
# My batch analysis - Dec 2025

./batch_timecourse.sh \
  --counts RawCountFile_rsemgenes_reformat.txt \
  --design timecourse_design.txt \
  --groups wt,ASC_KO,BST1_KO,GSDMD_KO \
  --method deseq2 \
  --cluster-genes \
  --n-clusters 6 \
  --run-enrichment \
  --enrichment-db reactome \
  --organism human \
  --output all_KO_lines_Dec2025 \
  --verbose
```

### 4. Monitor Progress
Use `--verbose` and redirect output:
```bash
./batch_timecourse.sh ... --verbose 2>&1 | tee batch_analysis.log
```

### 5. Backup Results
After completion:
```bash
tar -czf batch_results_backup_$(date +%Y%m%d).tar.gz batch_results/
```

## Post-Analysis

After batch processing completes:

### 1. Check All Outputs
```bash
# List all result directories
ls -lh batch_results/

# Check each has expected files
for dir in batch_results/*/; do
  echo "$dir:"
  ls "$dir"/limma/*.txt 2>/dev/null | wc -l
  ls "$dir"/deseq2/*.txt 2>/dev/null | wc -l
done
```

### 2. Compare Across Groups
```bash
# Compare number of significant genes
for dir in batch_results/*/; do
  group=$(basename "$dir")
  n_sig=$(tail -n +2 "$dir"/deseq2/deseq2_time_effect_results.txt | \
          awk '$NF < 0.05' | wc -l)
  echo "$group: $n_sig significant genes"
done
```

### 3. Create Summary Report
```R
# R script to summarize batch results
library(data.table)

results_dirs <- list.dirs("batch_results", recursive=FALSE)

summary_list <- lapply(results_dirs, function(dir) {
  group <- basename(dir)
  results_file <- file.path(dir, "deseq2/deseq2_time_effect_results.txt")
  
  if (file.exists(results_file)) {
    res <- fread(results_file)
    n_sig <- sum(res$padj < 0.05, na.rm=TRUE)
    return(data.frame(Group=group, N_significant=n_sig))
  }
})

summary_df <- do.call(rbind, summary_list)
print(summary_df)
write.csv(summary_df, "batch_results/summary.csv", row.names=FALSE)
```

## Summary

✅ **One command** runs multiple analyses
✅ **Consistent parameters** across all groups
✅ **Automatic organization** of results
✅ **Error handling** continues on failure
✅ **Progress tracking** shows status
✅ **Time efficient** for large-scale studies

The batch wrapper makes it easy to analyze complex experimental designs with multiple treatment groups or conduct comprehensive pairwise comparisons!

---

**Version**: 1.0 (Batch Processing)
**Date**: December 2025
