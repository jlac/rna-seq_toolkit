# pca_analysis.R

PCA of RNA-seq raw counts with top gene loadings per principal component.

Takes a raw counts matrix (genes × samples), filters low-count genes, normalizes
with DESeq2 VST or rlog, runs PCA via `prcomp`, and exports the top *N* genes
driving each of the first *K* PCs ranked by absolute rotation loading. Also
produces scatter plots, a scree plot, a heatmap of top loading genes, and all
coordinates and variance-explained tables.

---

## Requirements

R ≥ 4.0 with the following packages:

| Package | Source |
|---|---|
| `optparse` | CRAN |
| `DESeq2` | Bioconductor |
| `ggplot2` | CRAN |
| `ggrepel` | CRAN |
| `RColorBrewer` | CRAN |
| `pheatmap` | CRAN |

Install Bioconductor packages if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

---

## Quick start

```bash
# Minimal — counts only, no metadata coloring
Rscript pca_analysis.R -c counts.tsv

# With metadata for colored/shaped PCA
Rscript pca_analysis.R \
  -c counts.tsv \
  -m metadata.tsv \
  --color_by Condition \
  --shape_by Batch \
  --label_samples \
  -o my_pca_results
```

---

## Input file formats

### Counts matrix (`-c / --counts`, required)

Tab- or comma-delimited. First column is gene IDs; remaining columns are
integer raw counts, one per sample. Column headers are sample IDs.

```
Gene        Sample_1    Sample_2    Sample_3
BRCA1       412         389         510
TP53        8923        9102        8844
MYC         204         198         215
```

Duplicate gene IDs are handled automatically (see `--dup_method`).

### Sample metadata (`-m / --metadata`, optional)

Tab- or comma-delimited. First column is sample IDs that match the counts
column headers. Additional columns are grouping variables (condition, batch,
sex, etc.) used for PCA aesthetics.

```
SampleID    Condition   Batch
Sample_1    Treated     B1
Sample_2    Control     B1
Sample_3    Treated     B2
```

If metadata is provided but `--color_by` is not set, the first metadata
column after the sample ID is used for coloring by default.

---

## Options

### Input / output

| Flag | Description | Default |
|---|---|---|
| `-c, --counts` | Path to raw counts matrix | *required* |
| `-m, --metadata` | Path to sample metadata | none |
| `-o, --output_dir` | Output directory | `pca_results` |
| `--prefix` | Filename prefix for all outputs | `pca` |

### Duplicate gene handling

| Flag | Description | Default |
|---|---|---|
| `--dup_method` | How to resolve duplicate gene IDs: `sum` (combine counts across duplicates) or `make_unique` (append `_dup1`, `_dup2`, … suffixes) | `sum` |

**When to use each:**

- `sum` — appropriate when duplicates represent the same gene quantified
  across multiple features (e.g., multiple Ensembl IDs mapping to the same
  gene symbol). This is the standard approach for gene-level summarization.
- `make_unique` — appropriate when you want to preserve each row as a
  distinct feature, e.g., when duplicates represent genuinely different
  genomic loci that happen to share a name.

### Filtering

| Flag | Description | Default |
|---|---|---|
| `--min_count` | Minimum count threshold a gene must reach | `10` |
| `--min_samples` | Number of samples that must meet `--min_count` | `3` |
| `--min_samples_frac` | Fraction of samples that must meet threshold (overrides `--min_samples`) | none |

A gene is retained if at least `--min_samples` (or the fraction
`--min_samples_frac` of all samples) have a count ≥ `--min_count`.

### Normalization

| Flag | Description | Default |
|---|---|---|
| `--norm_method` | `vst` (variance-stabilizing transform) or `rlog` (regularized log) | `vst` |

Both use DESeq2 with `blind = TRUE` (sample-group–agnostic). VST is faster
and recommended for most datasets; rlog may perform better with very small
sample sizes (< 10).

### PCA parameters

| Flag | Description | Default |
|---|---|---|
| `--n_pcs` | Number of PCs to extract top genes for | `5` |
| `--n_top_genes` | Number of top genes per PC (by absolute loading) | `100` |
| `--n_top_var` | Restrict PCA to the top *N* most variable genes | all filtered genes |

### Plot aesthetics

| Flag | Description | Default |
|---|---|---|
| `--color_by` | Metadata column for point color | first metadata column |
| `--shape_by` | Metadata column for point shape | none |
| `--label_samples` | Add sample-name labels to PCA scatter | off |
| `--point_size` | Point size | `3` |
| `--width` | Plot width (inches) | `8` |
| `--height` | Plot height (inches) | `6` |

---

## Outputs

All files are written to `--output_dir` with the `--prefix` prepended.
Example tree with defaults:

```
pca_results/
├── pca_PC1_vs_PC2.pdf          # PCA scatter, PC1 vs PC2
├── pca_PC1_vs_PC2.png
├── pca_PC1_vs_PC3.pdf          # PCA scatter, PC1 vs PC3
├── pca_PC1_vs_PC3.png
├── pca_PC2_vs_PC3.pdf          # PCA scatter, PC2 vs PC3
├── pca_PC2_vs_PC3.png
├── pca_scree.pdf               # Scree plot (bar + cumulative line)
├── pca_scree.png
├── pca_top_loadings_heatmap.pdf  # Heatmap of union of top 20 genes per PC
├── pca_top100_genes_PC1.tsv    # Top 100 genes for PC1
├── pca_top100_genes_PC2.tsv    # Top 100 genes for PC2
├── pca_top100_genes_PC3.tsv    # ...
├── pca_top100_genes_PC4.tsv
├── pca_top100_genes_PC5.tsv
├── pca_top100_genes_all_PCs.tsv  # Combined wide-format table
├── pca_coordinates.tsv         # Sample PCA coordinates + metadata
├── pca_variance_explained.tsv  # Per-PC variance and cumulative %
└── pca_normalized_matrix.tsv   # Full VST/rlog normalized matrix
```

### Per-PC gene loading files

Each `*_genes_PCX.tsv` contains:

| Column | Description |
|---|---|
| `Gene` | Gene identifier |
| `Loading` | Signed rotation loading on that PC |
| `Abs_Loading` | Absolute value of loading |
| `Rank` | Rank (1 = highest absolute loading) |

The combined file places all PCs side by side with columns
`Rank`, `PC1_Gene`, `PC1_Loading`, `PC2_Gene`, `PC2_Loading`, etc.

### Heatmap

The heatmap shows the union of the top 20 loading genes from each PC
(row-scaled z-scores), clustered with Ward's method. If metadata is provided,
sample annotations are shown as a color bar along the top.

---

## Examples

```bash
# Basic run on Biowulf with lscratch
sinteractive --mem=16g --gres=lscratch:20
module load R/4.3
Rscript pca_analysis.R -c my_counts.tsv -o pca_out

# Full options
Rscript pca_analysis.R \
  -c raw_counts.csv \
  -m sample_info.csv \
  --color_by Treatment \
  --shape_by Sex \
  --dup_method sum \
  --min_count 10 \
  --min_samples_frac 0.25 \
  --norm_method vst \
  --n_top_genes 200 \
  --n_pcs 10 \
  --n_top_var 5000 \
  --label_samples \
  --prefix experiment1 \
  --width 10 \
  --height 8 \
  -o experiment1_pca

# rlog normalization for a small dataset (< 10 samples)
Rscript pca_analysis.R \
  -c small_counts.tsv \
  -m meta.tsv \
  --norm_method rlog \
  -o small_pca
```

---

## Method summary

1. **Load** — Read counts matrix; auto-detect TSV vs CSV.
2. **Deduplicate** — Resolve duplicate gene IDs by summing counts (default) or appending suffixes.
3. **Filter** — Remove genes not meeting the count threshold in enough samples.
4. **Normalize** — Apply DESeq2 VST or rlog (blind, design ~ 1) to stabilize variance.
5. **[Optional] Subset** — Restrict to top *N* most variable genes if `--n_top_var` is set.
6. **PCA** — `prcomp(t(matrix), center = TRUE, scale. = TRUE)` on the normalized, transposed matrix.
7. **Extract loadings** — Rank genes by absolute rotation loading per PC; export top *N*.
8. **Plot** — PCA scatter plots for PC pairs (1v2, 1v3, 2v3), scree plot, and heatmap of top loading genes.
9. **Export** — All coordinates, variance explained, normalized matrix, and gene lists as TSVs.
