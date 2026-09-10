# gsea_pathway_dotplot.py

Compare a hand-picked set of pathways across an arbitrary number of GSEA result
tables in a single dot plot.

- **Dot colour** = NES — red for positive (enriched/up), blue for negative (depleted/down)
- **Dot size** = adjusted p-value on a −log10 scale — bigger dot = more significant

Reads clusterProfiler, fgsea, and GSEA-desktop/GSEApy output without
configuration, in **CSV, TSV, TXT, XLSX, XLS, or XLSM** format.

---

## Requirements

```
python >= 3.8
pandas, numpy, matplotlib
openpyxl     # only for .xlsx/.xlsm input
scipy        # only for --order cluster
```

```bash
pip install pandas numpy matplotlib openpyxl scipy
```

---

## Quick start

```bash
python gsea_pathway_dotplot.py \
  -i "ADP-heptose=gsea_adp_heptose.tsv" \
  -i "ROSAH vs healthy donor=gsea_rosah.tsv" \
  -p pathways.txt \
  -o custom_dotplot_ADP-heptose_ROSAH.pdf --also png
```

---

## Input files

### Formats

Delimiter is sniffed automatically, so **`.tsv`, tab-delimited `.txt`, `.csv`,
and semicolon-delimited files all work with no extra flags**. Excel files are
read with `pandas.read_excel`; specify a sheet with `::`:

```bash
-i "Day 3=all_results.xlsx::D3_GSEA"
-i "Day 7=all_results.xlsx::D7_GSEA"
```

Omitting `::Sheet` reads the first sheet. Stray `Unnamed: 0` index columns are
dropped automatically.

### Labelling comparisons

```bash
-i "ROSAH vs healthy donor=path/to/file.tsv"   # explicit label (quote it if it has spaces)
-i path/to/gsea_rosah.tsv                      # label = filename stem -> "gsea_rosah"
```

Repeat `-i` once per comparison. **Columns appear in the plot in the order you
pass them.**

### Column auto-detection

| Field | Column names recognised |
|---|---|
| Pathway ID | `ID`, `pathway`, `NAME`, `Term`, `ont`, `GeneSet`, `Description` |
| Description | `Description`, `NAME`, `Term`, `pathway`, `ID` |
| NES | `NES`, `normalizedEnrichmentScore`, `normalized_enrichment_score` |
| Adjusted p | `p.adjust`, `padj`, `FDR q-val`, `qvalue`, `FDR`, `adj.P.Val`, `q.value`, … |
| Nominal p (fallback) | `pvalue`, `pval`, `NOM p-val`, `P.Value` |
| Set size | `setSize`, `SIZE`, `size` |

Fallback behaviour, each with a warning on stderr:

- no NES column → uses `ES`/`enrichmentScore`
- no adjusted-p column → uses the nominal p-value for dot size
- `FDR q-val` of exactly `0` (common in GSEA desktop output) → floored at the
  `--p-cap` value instead of producing `inf`

Override any of it explicitly with `--id-col`, `--desc-col`, `--nes-col`,
`--padj-col`.

---

## Pathway list

### File form (`-p pathways.txt`)

One pathway per line. Blank lines and `#` comments are ignored. An optional
**tab-separated** (or ` | ` separated) second field overrides the display label:

```
# hallmark sets of interest
HALLMARK_INFLAMMATORY_RESPONSE
HALLMARK_TNFA_SIGNALING_VIA_NFKB
HALLMARK_IL6_JAK_STAT3_SIGNALING
HALLMARK_INTERFERON_GAMMA_RESPONSE
HALLMARK_MYC_TARGETS_V1	MYC targets (v1)
```

### Inline form

```bash
-p "HALLMARK_INFLAMMATORY_RESPONSE,HALLMARK_TNFA_SIGNALING_VIA_NFKB"
```

### No list at all

Omit `-p` and the script plots the union of the top `--top-n` pathways (by
adjusted p) from each table — useful for a first look before you curate.

### How matching works

Names are normalised before comparison: case, underscores, spaces, and
punctuation are stripped, as are collection prefixes (`HALLMARK_`, `GOBP_`,
`GOCC_`, `GOMF_`, `GO_`, `REACTOME_`, `KEGG_`, `WP_`, `BIOCARTA_`, `PID_`,
`HP_`, `MSIGDB_`, `C2_`, `C5_`, `C7_`). Both the ID and the Description column
are checked.

So `HALLMARK_TNFA_SIGNALING_VIA_NFKB` matches `Tnfa Signaling Via Nfkb` and
`TNFA SIGNALING VIA NFKB` across tables that name things differently. Add
`--partial-match` for substring fallback — handy for long GO/Reactome names.
If several rows match, the one with the smallest adjusted p wins.

Pathways found in no table are reported on stderr and drawn as an empty row
unless you pass `--drop-missing`.

---

## Options

### Input / output

| Option | Default | Description |
|---|---|---|
| `-i`, `--input [LABEL=]FILE[::SHEET]` | *required* | GSEA table; repeat once per comparison. Order sets column order. |
| `-p`, `--pathways FILE\|LIST` | — | Pathway list file, or comma-separated names. Omit to use `--top-n`. |
| `--top-n N` | `10` | Without `-p`: union of the top N pathways per table by adjusted p. |
| `-o`, `--output FILE` | `gsea_dotplot.pdf` | Output path. Extension sets the format (`.pdf`, `.png`, `.svg`, `.tiff`, `.eps`). |
| `--also EXT [EXT …]` | — | Write extra formats alongside, e.g. `--also png svg`. |
| `--save-table TSV` | — | Dump the exact plotted values (pathway, comparison, NES, padj, setSize, source ID). |

### Column overrides

| Option | Description |
|---|---|
| `--id-col NAME` | Force the pathway ID column. |
| `--desc-col NAME` | Force the description column. |
| `--nes-col NAME` | Force the NES column. |
| `--padj-col NAME` | Force the adjusted p-value column (e.g. `--padj-col qvalue`). |

### Matching, filtering, ordering

| Option | Default | Description |
|---|---|---|
| `--partial-match` | off | Substring fallback when exact normalised matching fails. |
| `--use-description` | off | Label rows from the Description column instead of the ID. |
| `--max-padj X` | show all | Hide dots with `padj > X` (e.g. `0.05` to blank non-significant cells). |
| `--min-nes X` | show all | Hide dots with `\|NES\| < X`. |
| `--drop-missing` | off | Drop pathways absent from every table instead of leaving a blank row. |
| `--drop-empty-rows` | off | Drop rows left with no dots after `--max-padj`/`--min-nes`. |
| `--mark-missing` | off | Draw a grey **×** where a pathway is absent from a given table. |
| `--order {input,nes,padj,alpha,cluster}` | `input` | Row order. `input` = as listed, first at top. `nes` = mean NES descending. `padj` = best adjusted p. `alpha` = alphabetical. `cluster` = hierarchical clustering on the NES matrix (needs scipy). |
| `--reverse-rows` | off | Flip the row order after sorting. |

### Colour (NES)

| Option | Default | Description |
|---|---|---|
| `--cmap NAME` | `RdBu_r` | Any diverging matplotlib colormap. `RdBu_r` gives blue = down, red = up. Alternatives: `coolwarm`, `bwr`, `seismic`. |
| `--nes-limit X` | from data | Symmetric colour limit `[-X, +X]`. Set explicitly (e.g. `3`) to keep the scale identical across figures. Default rounds the max \|NES\| up to the nearest 0.5, minimum 1. |
| `--color-label TEXT` | `NES` | Colourbar title. |

### Size (adjusted p)

| Option | Default | Description |
|---|---|---|
| `--p-cap X` | `10` | −log10(padj) at which dot size saturates; `10` means everything ≤ 1e-10 gets the largest dot. Also the floor applied to `FDR q-val = 0`. |
| `--min-size A` | `25` | Marker **area** (pt²) for the least significant dot. |
| `--max-size A` | `420` | Marker area for the most significant dot. |
| `--size-gamma G` | `1.0` | Exponent on the size scale. `G < 1` (e.g. `0.6`) inflates small dots — useful when everything is highly significant and the dots all look alike. |
| `--legend-p P [P …]` | `0.5 0.05 0.001 1e-5 1e-10` | Values shown in the size legend. Entries below the `--p-cap` floor are dropped. |
| `--size-label TEXT` | `p.adjust` | Size legend title (e.g. `FDR`, `q-value`). |

### Layout and typography

| Option | Default | Description |
|---|---|---|
| `--label-style {upper,title,sentence,raw}` | `upper` | Row label rendering after prefix stripping and `_` → space. |
| `--title TEXT` | — | Plot title. |
| `--col-width IN` | `0.95` | Inches per comparison column. |
| `--row-height IN` | `0.72` | Inches per pathway row. |
| `--figsize W H` | auto | Override the whole figure size in inches. |
| `--left-margin IN` | auto | Space for pathway labels; raise it if long names are clipped. |
| `--bottom-margin IN` | auto | Space for rotated comparison labels. |
| `--x-rotation DEG` | `45` | Rotation of the comparison labels; `0` centres them horizontally. |
| `--tick-fontsize N` | `13` | Pathway and comparison label size. |
| `--legend-fontsize N` | `11` | Colourbar and size-legend text. |
| `--title-fontsize N` | `14` | Title size. |
| `--font NAME` | matplotlib default | Font family, e.g. `Arial`, `Helvetica`. |
| `--edge-color C` | `black` | Dot outline colour (`none` for no outline). |
| `--edge-width W` | `0.8` | Dot outline width. |
| `--dpi N` | `300` | Raster resolution. |
| `--transparent` | off | Transparent background. |

PDF/EPS output uses `fonttype 42` and SVG uses `svg.fonttype: none`, so text
stays as editable text in Illustrator/Inkscape rather than being converted to
paths.

---

## Recipes

**Reproduce the two-column hallmark figure**

```bash
python gsea_pathway_dotplot.py \
  -i "ADP-heptose=gsea_adp.tsv" \
  -i "ROSAH vs healthy donor=gsea_rosah.tsv" \
  -p hallmark_pathways.txt \
  --nes-limit 3 \
  -o custom_dotplot_ADP-heptose_ROSAH.pdf --also png
```

**Multi-tab workbook, one column per timepoint, clustered rows**

```bash
python gsea_pathway_dotplot.py \
  -i "0 h=gsea.xlsx::T0" -i "6 h=gsea.xlsx::T6" -i "24 h=gsea.xlsx::T24" \
  -p pathways.txt --order cluster --mark-missing \
  -o timecourse_dotplot.pdf
```

**Only significant dots, non-significant cells left blank**

```bash
python gsea_pathway_dotplot.py -i a.tsv -i b.tsv -p pathways.txt \
  --max-padj 0.05 --drop-empty-rows -o sig_only.pdf
```

**Exploratory pass with no curated list**

```bash
python gsea_pathway_dotplot.py -i a.tsv -i b.tsv -i c.tsv \
  --top-n 15 --order nes --save-table plotted_values.tsv -o explore.png
```

**Everything is p < 1e-10 and the dots all look the same**

```bash
--p-cap 20 --size-gamma 0.6 --legend-p 0.05 1e-5 1e-10 1e-20
```

**Long Reactome names running off the left edge**

```bash
--partial-match --label-style sentence --left-margin 5.0 --tick-fontsize 10
```

---

## Notes and gotchas

- Quote labels containing spaces or `=`: `-i "ROSAH vs healthy donor=file.tsv"`.
  A bare `=` inside a path is only treated as a separator when the full string
  isn't an existing file.
- With `--order input`, the first pathway in your list is drawn at the **top**.
- The colour scale is symmetric about zero by construction, so white always
  means NES ≈ 0. Fix `--nes-limit` when you want several figures directly
  comparable.
- Dot size is a function of adjusted p only. If you'd rather encode set size or
  leading-edge count, `--save-table` gives you the tidy dataframe to work from.
- Warnings (unmatched pathways, fallback columns, filtered dots) go to stderr;
  redirect with `2> log.txt` if you're scripting batches.
