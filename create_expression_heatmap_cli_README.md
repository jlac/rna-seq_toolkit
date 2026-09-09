# create_expression_heatmap_cli.R

A command-line R script for making four kinds of figures from RNA-seq results:

| Plot type | Input | What it shows |
|---|---|---|
| `heatmap` | expression matrix + design | Genes × samples (or × group averages) |
| `violin` | expression matrix + design | Distribution per gene, split by group |
| `boxplot` | expression matrix + design | Same as violin, box-and-whisker |
| `volcano` | DE results table | log2FC vs −log10(p), with optional gene highlighting |

Every option is set on the command line. `--help` prints the full list with defaults.

---

## Installation

```r
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages(c("circlize", "RColorBrewer", "viridis", "optparse",
                   "ggplot2", "reshape2"))
```

Two optional packages:

- **ggrepel** — non-overlapping gene labels on volcano plots. Without it the script falls back to plain text labels that may overlap, and prints a note.
- **ggsignif** — significance brackets on violin/boxplots (`--add_stats`). Without it that flag is skipped with a warning.

```r
install.packages(c("ggrepel", "ggsignif"))
```

The script does not pin a ComplexHeatmap version. It uses long-stable API, so any 2.x release should work.

---

## Input files

### Expression matrix (`-e`)

Tab-delimited. First column = gene IDs, remaining columns = samples. Values should already be normalized (TPM, CPM, FPKM, or normalized counts).

```
gene_id    sample1    sample2    sample3
A1BG       12.4       8.1        15.2
A1CF       0.0        0.3        0.1
```

Duplicate gene IDs are handled by `--duplicate_genes` rather than failing. Your file has 149 duplicates, so this matters.

### Design file (`-d`)

Tab-delimited. First column must be named `Sample` and match the expression matrix column names. Any number of additional metadata columns.

```
Sample     Group                       BigGroup         Pop
sample1    Filaria-neg_Media_pop-neg   Media_pop-neg    PopNeg
sample2    Filaria-neg_Media_Pop-POS   Media_Pop-POS    PopPos
```

Your file has `Group`, `BigGroup`, and `Pop` — these are the values you can pass to `-f`, `-a`, `--group_by`, `--average_by`, and `--split_by`.

### DE results table (`--de_file`, volcano only)

Tab- or comma-delimited (auto-sniffed from the header line). Needs a gene ID column, a log2 fold-change column, and a p-value column. DESeq2, edgeR, and limma output all work without specifying column names.

---

## Quick start

```bash
# See all options
Rscript create_expression_heatmap_cli.R --help

# Heatmap, top 1000 variable genes, PopNeg vs PopPos
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt \
  -d evitta_design.txt \
  -o heatmap.pdf

# Violin plot of specific genes
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt -d evitta_design.txt \
  --plot_type violin -g "IFNG,IL6,TNF" --group_by Pop -o violin.pdf

# Volcano plot with highlighted genes (no -e or -d needed)
Rscript create_expression_heatmap_cli.R \
  --plot_type volcano --de_file deseq2_results.txt \
  --highlight_genes "IFNG,IL6,TNF" -o volcano.pdf
```

**Quoting:** comma-separated lists go inside **one** set of quotes.
`-v "A,B"` is correct. `-v "A","B"` is not — the shell splits that into two arguments and only the first reaches the script.

---

## Option reference

### Input and output

| Option | Type | Default | Description |
|---|---|---|---|
| `-e`, `--expression` | file | `Reformat_TPMCountFile_rsemgenes.txt` | Expression matrix. Not read in volcano mode. |
| `-d`, `--design` | file | `evitta_design.txt` | Design/metadata file. Not read in volcano mode. |
| `-o`, `--output` | file | `expression_heatmap.pdf` | Output figure. Extension decides format: `.png` writes PNG, anything else writes PDF. |
| `--save_matrix` | flag | off | Also write the data behind the figure to a text file. |
| `--matrix_file` | file | `<output>_matrix.txt` | Custom path for the saved data. In volcano mode the default suffix is `_data.txt`. |
| `--duplicate_genes` | choice | `make_unique` | How to resolve repeated gene IDs. See below. |
| `--plot_type` | choice | `heatmap` | `heatmap`, `violin`, `boxplot`, or `volcano`. |
| `--verbose` | flag | off | Print extra progress detail. |

**`--duplicate_genes` values:**

- `make_unique` — append `_1`, `_2`, … to repeated IDs. Keeps every row; safest default.
- `sum` — add the values across duplicate rows.
- `mean` — average across duplicate rows.
- `first` — keep the first occurrence, discard the rest.

**What `--save_matrix` writes:** the matrix as it exists at plot time — after gene selection, filtering, transformation, scaling, and group averaging. Not the raw input. In volcano mode it writes gene, log2FC, pvalue, neglog10p, category, and highlight status.

---

### Sample filtering — heatmap, violin, boxplot

| Option | Type | Default | Description |
|---|---|---|---|
| `-f`, `--filter_column` | column | `Pop` | Design column to filter on. Pass `none` to keep all samples. |
| `-v`, `--filter_values` | list | `PopNeg,PopPos` | Values in that column to keep. |

The default drops your `Ignore` samples. To keep everything:

```bash
-f none
```

---

### Gene selection — heatmap, violin, boxplot

| Option | Type | Default | Description |
|---|---|---|---|
| `-n`, `--n_genes` | integer | `1000` | Use the N most variable genes. |
| `-g`, `--gene_list` | list | none | Specific genes, comma-separated. Overrides `-n`. |
| `--gene_file` | file | none | Specific genes, one per line. Overrides `-n` and `-g`. |

Precedence: `--gene_file` > `-g` > `-n`.

Genes in your list that aren't in the matrix produce a warning naming the first 10 and are skipped; the plot still gets made from the rest.

**Violin and boxplot require an explicit gene list.** `-n` is rejected in those modes — a violin panel per gene only makes sense for a handful of named genes.

**How top-variable genes are chosen:** variance is computed on the values *as supplied*, before `-t` is applied. On TPM data that means selection is driven by absolute magnitude, so highly-expressed genes dominate. If you want variance ranked on the log scale instead, pre-log your input matrix and run with `-t none`.

---

### Transformation and scaling — heatmap, violin, boxplot

| Option | Type | Default | Description |
|---|---|---|---|
| `-t`, `--transform` | choice | `log2` | `log2`, `log10`, `zscore`, or `none`. |
| `-p`, `--pseudocount` | number | `1` | Added before log transform, to keep zeros finite. |
| `-s`, `--scale` | choice | `zscore` | Per-gene row scaling: `zscore` or `none`. |

These are two separate steps, applied in order: `-t` first, then `-s`.

The defaults (`-t log2 -s zscore`) are the standard choice for expression heatmaps and what you want most of the time.

Two combinations to avoid:

- **`-t none -s zscore`** — this is what caused your earlier clustering crash. Genes with identical values across your samples have zero standard deviation, and dividing by it yields `Inf`, which `hclust` rejects. The script now detects and drops zero-variance genes before scaling, reporting how many, so this no longer errors — but log-transforming first is still the better fix.
- **`-t zscore -s zscore`** — z-scores the data twice. Harmless but pointless; use `-t log2 -s zscore`.

---

### Annotations — heatmap only

| Option | Type | Default | Description |
|---|---|---|---|
| `-a`, `--annotations` | list | `BigGroup,Pop` | Design columns to draw as colored bars above the heatmap. Pass `none` for no annotations. |
| `--annotation_fontsize` | number | `10` | Annotation label font size. |

Colors are assigned automatically: ColorBrewer `Set2` for up to 8 levels, `Set3` for up to 12, `rainbow()` beyond that.

---

### Group averaging — heatmap only

| Option | Type | Default | Description |
|---|---|---|---|
| `--average_groups` | flag | off | Collapse samples into one column per group. |
| `--average_by` | column | `Group` | Design column defining the groups. |
| `--average_function` | choice | `mean` | `mean` or `median`. |

Averaging runs *after* transformation and scaling, so with the defaults you are averaging z-scores, not raw TPM.

Applies to heatmaps only — it's ignored in violin and boxplot mode, where showing the spread within a group is the entire point of the plot.

`median` is the more robust choice when a group contains an outlier sample.

---

### Clustering — heatmap only

| Option | Type | Default | Description |
|---|---|---|---|
| `--cluster_rows` | TRUE/FALSE | `TRUE` | Cluster genes. |
| `--cluster_columns` | TRUE/FALSE | `TRUE` | Cluster samples. |
| `--no_cluster_rows` | flag | — | Shorthand for `--cluster_rows FALSE`. |
| `--no_cluster_columns` | flag | — | Shorthand for `--cluster_columns FALSE`. |

Turning off column clustering keeps samples in design-file order, which is often what you want when groups are already sorted sensibly.

---

### Heatmap appearance

| Option | Type | Default | Description |
|---|---|---|---|
| `--show_row_names` | flag | off | Print gene names down the side. Only legible for roughly ≤100 genes. |
| `--row_fontsize` | number | `6` | Gene name size. |
| `--col_fontsize` | number | `8` | Sample name size. |
| `--color_scheme` | choice | `blue_white_red` | Also `viridis`, `green_black_red`, `purple_white_orange`. |
| `--color_min` | number | `-2` | Value mapped to the low end of the color scale. |
| `--color_max` | number | `2` | Value mapped to the high end. |
| `--split_by` | column | none | Split the heatmap into panels by a design column. |

`--color_min`/`--color_max` defaults of ±2 assume z-scored data. If you run `-s none`, set them to match your actual value range or the heatmap will be almost entirely saturated.

`--color_scheme viridis` ignores `--color_min`/`--color_max`; it maps across the data range.

---

### Violin and boxplot options

| Option | Type | Default | Description |
|---|---|---|---|
| `--group_by` | column | `Group` | Design column defining the groups compared within each gene. |
| `--violin_groups` | list | all | Restrict to specific groups from that column. |
| `--violin_colors` | list | auto | Hex colors, comma-separated, one per group. |
| `--add_points` | flag | off | Overlay individual samples as jittered points. |
| `--add_stats` | flag | off | Add a significance bracket. Requires **ggsignif**. |

`-f`/`-v` and `--violin_groups` both subset samples and can be combined: `-f`/`-v` runs first on the full design, then `--violin_groups` narrows further within `--group_by`.

`--add_points` is worth using by default with small n — it shows whether a wide violin reflects real spread or just three scattered samples.

**`--add_stats` limitation:** it compares only the first two levels of `--group_by` with a t-test. With more than two groups, or if you need a different test or paired comparisons, run the statistics separately and annotate the figure yourself.

---

### Volcano plot: input columns

| Option | Type | Default | Description |
|---|---|---|---|
| `--de_file` | file | none | **Required.** DE results table. |
| `--de_gene_col` | column | auto | Gene ID column. |
| `--de_lfc_col` | column | auto | log2 fold-change column. |
| `--de_pval_col` | column | auto | P-value column. |

Auto-detection is case-insensitive and searches in this order:

- **Gene:** `gene`, `gene_id`, `gene_name`, `geneid`, `genes`, `symbol`, `gene_symbol`, `id` — falling back to the first column.
- **log2FC:** `log2FoldChange`, `logFC`, `log2FC`, `log2_fold_change`, `lfc`, `fold_change_log2`.
- **P-value:** `padj`, `FDR`, `adj.P.Val`, `qvalue`, `q_value`, `adj_pvalue`, `p_adj`, then `pvalue`, `PValue`, `P.Value`, `p_value`, `pval`.

Adjusted p-values are preferred over raw ones, and the y-axis label changes to match whichever was used. The script prints all three chosen columns — check that line before trusting the figure. If detection fails or picks wrong, name the columns explicitly.

---

### Volcano plot: thresholds

| Option | Type | Default | Description |
|---|---|---|---|
| `--lfc_threshold` | number | `1` | Absolute log2FC cutoff. Sets the vertical dashed lines. |
| `--pval_threshold` | number | `0.05` | P-value cutoff. Sets the horizontal dashed line. |

These control both the guide lines and the Down / Not significant / Up coloring. `--lfc_threshold 0.58` is a 1.5-fold change; `1` is 2-fold.

---

### Volcano plot: highlighting

| Option | Type | Default | Description |
|---|---|---|---|
| `--highlight_genes` | list | none | Genes to highlight, comma-separated. |
| `--highlight_file` | file | none | Genes to highlight, one per line. Takes precedence over `--highlight_genes`. |
| `--highlight_color` | hex | `#D62728` | Color for highlighted points. |
| `--highlight_only` | flag | off | Grey out everything else instead of coloring by category. |
| `--label_top_n` | integer | `0` | Also label the N most significant genes. |
| `--no_labels` | flag | off | Keep highlight coloring, drop the text labels. |

Highlighted genes are drawn on top of the rest, at roughly twice the point size, with a thin black outline, and labeled.

`--highlight_only` is the mode to use when the point is *where a gene set sits* rather than which genes are significant — a pathway or panel against a grey background.

Gene IDs must match the DE table exactly. The script reports how many of your requested genes it found, and names the ones it didn't, distinguishing two cases:

- genes absent from the DE table entirely (usually an ID mismatch — symbols vs Ensembl IDs)
- genes present in the table but dropped for `NA` values

---

### Volcano plot: appearance

| Option | Type | Default | Description |
|---|---|---|---|
| `--volcano_colors` | list | `#3C6FBF,#BFBFBF,#C0392B` | Down, not-significant, up. Must be exactly three. |
| `--point_size` | number | `1.2` | Base point size. |
| `--label_size` | number | `3.5` | Gene label text size. |
| `--plot_title` | text | auto | Custom title. Default states the thresholds used. |
| `--xlim` | `min,max` | auto | X-axis limits. |
| `--ylim` | `min,max` | auto | Y-axis limits. |

`--xlim`/`--ylim` use `coord_cartesian`, so they zoom without dropping data.

---

### Figure size

| Option | Type | Default | Description |
|---|---|---|---|
| `--width` | inches | `10` | Figure width. |
| `--height` | inches | `12` | Figure height. |
| `--dpi` | number | `300` | Resolution, PNG output only. |

In violin/boxplot mode, width auto-expands with the number of genes if `--width` is too small; the value you pass acts as a floor.

The `10 × 12` default is portrait, sized for heatmaps. Volcano plots and violin plots usually want something wider and shorter — try `--width 8 --height 6`.

---

## Behavior worth knowing

**Volcano mode ignores `-e` and `-d`.** It reads only `--de_file` and exits before the expression matrix is touched. Gene selection, filtering, transformation, scaling, averaging, clustering, and annotation options do nothing in volcano mode.

**NA handling in volcano mode.** Rows with `NA` log2FC or p-value are dropped, and the count is reported. In DESeq2 output, `NA` in `padj` usually means the gene was removed by independent filtering or flagged as a count outlier — those genes have no significance call and no y-position, so there's nothing to plot.

**P-values of exactly zero.** `-log10(0)` is infinite. Those points are floored at the smallest nonzero p-value in the table, with a note in the output. Their height is a floor, not a measurement — don't read anything into how tall they are relative to each other.

**Zero-variance genes.** Under `-s zscore`, genes with no variation across the selected samples are dropped before scaling and the count is reported. With a short gene list this can remove genes you specifically asked for, so check that line.

**Order of operations** for heatmap/violin/boxplot:

```
load → filter samples → select genes → transform (-t) → scale (-s)
     → average groups → save matrix → plot
```

---

## Troubleshooting

**`duplicate 'row.names' are not allowed`**
An older version of the script. The current one handles duplicates via `--duplicate_genes`.

**`NA/NaN/Inf in foreign function call (arg 10)` during clustering**
Zero-variance genes reaching the scaling step. Use `-t log2` (recommended), or `-s none`, or disable clustering.

**`None of the specified genes found in expression data!`**
ID namespace mismatch — check whether your list uses symbols and the matrix uses Ensembl IDs, or vice versa. Compare: `head -3 Reformat_TPMCountFile_rsemgenes.txt | cut -f1` against `head -3 your_genes.txt`.

**`Could not auto-detect the ... column`**
Name it explicitly with `--de_gene_col` / `--de_lfc_col` / `--de_pval_col`. The error message lists the columns actually present in your file.

**`Violin/boxplot requires a specific gene list`**
Use `-g` or `--gene_file`; `-n` isn't valid in those modes.

**`Averaging column 'X' not found in design`** / **`Group column 'X' not found`**
Check available columns: `head -1 evitta_design.txt`. Yours are `Sample`, `Group`, `BigGroup`, `Pop`.

**Only one column after averaging**
Every remaining sample falls in the same `--average_by` group. Either your `-f`/`-v` filter was too narrow, or you're averaging by a column that doesn't vary within the filtered set.

**Volcano labels overlap**
Install ggrepel: `install.packages("ggrepel")`.

**Nothing highlighted on the volcano plot**
The script reports how many highlight genes it matched. If it's zero, the IDs don't match the DE table's gene column.

---

## Worked examples

```bash
# Filaria-neg vs Filaria-POS, CMV Pop-POS samples only
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt -d evitta_design.txt \
  -f Group -v "Filaria-neg_CMV_Pop-POS,Filaria-POS_CMV_Pop-POS" \
  -n 1000 -a Group -o cmv_comparison.pdf --save_matrix

# DEG list, gene names shown, keeping design order
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt -d evitta_design.txt \
  --gene_file deg_list.txt --show_row_names --row_fontsize 8 \
  --no_cluster_columns -o deg_heatmap.pdf

# One column per treatment group, medians
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt -d evitta_design.txt \
  --average_groups --average_by BigGroup --average_function median \
  -n 1000 -a BigGroup -o averaged.pdf

# Violin plot, two groups, points shown
Rscript create_expression_heatmap_cli.R \
  -e Reformat_TPMCountFile_rsemgenes.txt -d evitta_design.txt \
  --plot_type violin -g "IFNG,IL6,TNF,IL1B" \
  --group_by Pop -f Pop -v "PopNeg,PopPos" \
  --add_points --width 9 --height 6 -o violin.pdf

# Volcano: pathway genes against a grey background
Rscript create_expression_heatmap_cli.R \
  --plot_type volcano --de_file deseq2_results.txt \
  --highlight_file pathway_genes.txt --highlight_only \
  --lfc_threshold 0.58 --pval_threshold 0.05 \
  --width 8 --height 6 -o pathway_volcano.pdf --save_matrix

# Volcano: standard three-color, top 15 labeled, custom columns
Rscript create_expression_heatmap_cli.R \
  --plot_type volcano --de_file edger_results.csv \
  --de_lfc_col logFC --de_pval_col FDR \
  --label_top_n 15 --plot_title "Filaria-POS vs Filaria-neg" \
  --width 8 --height 6 -o volcano_top15.pdf
```
