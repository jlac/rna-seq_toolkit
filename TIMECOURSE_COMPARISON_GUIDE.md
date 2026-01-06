# Cross-Group Comparison Guide

## Overview

After running batch analyses with multiple groups, the `compare_batch_results.R` script automatically runs to identify:

1. **Genes that change similarly across groups**
2. **Pathways enriched in multiple groups**
3. **Direction of changes** (up/down) across groups
4. **Group-specific vs. shared responses**

## Automatic Execution

When you run batch analyses with `--groups`, the comparison script runs automatically if 2+ groups complete successfully:

```bash
./batch_timecourse.sh \
  --counts counts.txt \
  --design design.txt \
  --groups wt,ASC_KO,BST1_KO,GSDMD_KO \
  --method deseq2 \
  --cluster-genes \
  --run-enrichment \
  -o batch_results
```

**After all groups complete:**
```
===============================================================================
BATCH ANALYSIS COMPLETE
================================================================================
  Total groups: 4
  Completed: 4
  Failed: 0

===============================================================================
RUNNING CROSS-GROUP COMPARISON
===============================================================================

✓ Cross-group comparison complete!
  Results saved to: batch_results/cross_group_comparison
```

## Output Structure

```
batch_results/
├── wt/                          ← Individual group results
├── ASC_KO/
├── BST1_KO/
├── GSDMD_KO/
└── cross_group_comparison/      ← NEW: Comparison results
    ├── COMPARISON_SUMMARY.txt   ← Main summary report
    ├── tables/
    │   ├── significant_genes_per_group.txt
    │   ├── gene_overlap_summary.txt
    │   ├── genes_shared_in_all_groups.txt
    │   ├── genes_unique_to_wt.txt
    │   ├── genes_unique_to_ASC_KO.txt
    │   ├── direction_of_change_summary.txt
    │   └── pathway_comparison.txt
    └── figures/
        ├── gene_overlap_venn.pdf         ← 2-5 groups
        ├── gene_overlap_upset.pdf        ← 3+ groups
        ├── shared_genes_lfc_heatmap.pdf  ← Log fold changes
        └── pathway_comparison_heatmap.pdf ← Pathway patterns
```

## Key Analyses

### 1. Gene Overlap

**Questions answered:**
- How many genes are significant in each group?
- Which genes are shared across all groups?
- Which genes are unique to each group?
- Do shared genes change in the same direction?

**Files:**
- `significant_genes_per_group.txt` - Counts per group
- `gene_overlap_summary.txt` - Shared vs. unique
- `genes_shared_in_all_groups.txt` - Common responders
- `genes_unique_to_*.txt` - Group-specific genes

**Visualizations:**
- **Venn diagram** (2-5 groups): Visual overlap
- **UpSet plot** (3+ groups): Complex overlaps
- **Heatmap**: Log fold changes for shared genes

### 2. Direction of Change

**Questions answered:**
- Do genes go up or down consistently?
- Which genes change oppositely between groups?
- Are there concordant responses?

**Categories identified:**
- **Concordant** - Same direction in all groups
- **Discordant** - Different directions
- **Up in all** - Universally upregulated
- **Down in all** - Universally downregulated

**Example output:**
```
Direction of change (shared genes):
                         Category Count
  Concordant (all same direction)   234
Discordant (different directions)    56
                   Up in all groups   145
                 Down in all groups    89
```

### 3. Pathway Comparison

**Questions answered:**
- Which pathways are enriched in multiple groups?
- Which pathways are group-specific?
- What biological processes are shared vs. unique?

**Files:**
- `pathway_comparison.txt` - All pathways with group membership
- **Heatmap**: Pathway enrichment patterns

**Identifies:**
- **Universal pathways** - Enriched in all groups
- **Shared pathways** - Enriched in 2+ groups
- **Group-specific pathways** - Unique to one group

## Example Use Cases

### Use Case 1: Identify Common Response

**Goal:** Find genes that respond to treatment regardless of genetic background

**Look at:**
- `genes_shared_in_all_groups.txt`
- `direction_of_change_summary.txt`

**Interpretation:**
- Genes in "Up in all groups" = Core response
- Genes in "Concordant" = Consistent biology
- Use these for follow-up validation

### Use Case 2: Find Compensation Mechanisms

**Goal:** Identify genes that change oppositely in knockout vs. wildtype

**Look at:**
- `shared_genes_lfc_heatmap.pdf`
- Filter for "Discordant" genes
- `direction_of_change_summary.txt`

**Interpretation:**
- Genes up in WT but down in KO = Potential compensation
- Opposite patterns suggest regulatory relationships

### Use Case 3: Characterize KO-Specific Biology

**Goal:** Understand what's unique about each knockout

**Look at:**
- `genes_unique_to_ASC_KO.txt`
- `pathway_comparison.txt` (pathways with N_groups = 1)
- Heatmap clusters

**Interpretation:**
- Unique genes = KO-specific effects
- Unique pathways = Specialized biology
- Compare across KOs to find patterns

### Use Case 4: Universal vs. Selective Pathways

**Goal:** Find biological processes that are universally affected vs. context-dependent

**Look at:**
- `pathway_comparison_heatmap.pdf`
- `pathway_comparison.txt` (sort by N_groups)

**Interpretation:**
- Pathways in all groups = Core mechanisms
- Pathways in few groups = Context-specific
- Use for therapeutic target selection

## Manual Execution

You can also run the comparison manually:

```bash
Rscript compare_batch_results.R \
  --input-dir batch_results \
  --output-dir batch_results/cross_group_comparison \
  --method deseq2 \
  --fdr-cutoff 0.05 \
  --top-n 20 \
  --verbose
```

**Options:**
- `--input-dir` - Directory with batch results (required)
- `--output-dir` - Where to save comparison [default: batch_comparison]
- `--method` - limma or deseq2 [default: deseq2]
- `--fdr-cutoff` - Significance threshold [default: 0.05]
- `--top-n` - Top pathways to show [default: 20]
- `--verbose` - Print detailed output

## Reading the Summary Report

The `COMPARISON_SUMMARY.txt` file contains:

### Section 1: Overview
- Analysis date
- Number of groups
- Groups analyzed
- Parameters used

### Section 2: Gene Analysis
```
Significant genes per group:
     Group N_significant
        wt          1234
    ASC_KO          2456
   BST1_KO          1876
  GSDMD_KO          2134

Gene overlap:
                              Category Count
Total unique genes across all groups  4567
      Genes significant in ALL groups   234
              Genes unique to wt        123
          Genes unique to ASC_KO        456
         Genes unique to BST1_KO        234
        Genes unique to GSDMD_KO        345
```

### Section 3: Direction Analysis
```
Direction of change (shared genes):
                         Category Count
  Concordant (all same direction)   178
Discordant (different directions)    56
                   Up in all groups   134
                 Down in all groups    44
```

### Section 4: Pathway Analysis
```
Top 20 shared pathways:
                                Pathway N_groups                    Groups
                  Immune System Response        4  wt, ASC_KO, BST1_KO, GSDMD_KO
             Interferon Signaling        4  wt, ASC_KO, BST1_KO, GSDMD_KO
                   Cytokine Signaling        3       ASC_KO, BST1_KO, GSDMD_KO
                  ...
```

## Visualizations Explained

### 1. Venn Diagram (2-5 groups)
- **Overlapping circles** = Shared genes
- **Non-overlapping regions** = Unique genes
- **Numbers** = Gene counts
- **Colors** = Different groups

**Interpretation:**
- Large overlap = Similar responses
- Small overlap = Divergent responses
- Look at intersection numbers

### 2. UpSet Plot (3+ groups)
- **Horizontal bars** = Total genes per group
- **Vertical bars** = Intersection sizes
- **Dots below** = Which groups are in each intersection
- **Connected dots** = Multiple group intersections

**Interpretation:**
- Tallest bar = Largest intersection
- Many tall bars = Complex overlap patterns
- Use to identify 2-way, 3-way, etc. overlaps

### 3. Log Fold Change Heatmap
- **Rows** = Shared genes
- **Columns** = Groups
- **Color** = Direction/magnitude
  - Red = Upregulated
  - Blue = Downregulated
  - White = No change

**Interpretation:**
- Same color across row = Concordant
- Opposite colors = Discordant
- Clusters = Co-regulated genes
- Look for patterns in hierarchical clustering

### 4. Pathway Comparison Heatmap
- **Rows** = Pathways
- **Columns** = Groups
- **Color** = Enriched (dark) or not (white)
- **Binary** = Present/absent

**Interpretation:**
- Horizontal lines = Universal pathways
- Vertical patterns = Group similarities
- Scattered dots = Group-specific pathways
- Clusters = Related groups

## Advanced Analysis Tips

### Combining with Cluster Information

If you ran with `--cluster-genes`, pathways are analyzed per cluster:

1. Look at `pathway_comparison.txt`
2. Note the cluster patterns (e.g., "Cluster 1_Transient Up")
3. Compare which clusters are enriched across groups
4. Identify if similar temporal patterns show similar pathways

### Finding Compensatory Mechanisms

1. Load `shared_genes_lfc_heatmap.pdf`
2. Find genes with opposite colors (red in WT, blue in KO)
3. Check if these are in `direction_of_change_summary.txt` as "Discordant"
4. Look up these genes in individual group results for details

### Meta-Analysis Strategy

1. **Start broad:** Look at `COMPARISON_SUMMARY.txt`
2. **Identify patterns:** Use heatmaps for visual overview
3. **Drill down:** Use gene/pathway lists for specifics
4. **Cross-reference:** Link back to individual group results
5. **Validate:** Pick top candidates for experimental follow-up

### Customizing Thresholds

For more or fewer genes/pathways:

```bash
# Stricter (fewer genes)
Rscript compare_batch_results.R \
  --input-dir batch_results \
  --fdr-cutoff 0.01 \
  --top-n 10

# Looser (more genes)
Rscript compare_batch_results.R \
  --input-dir batch_results \
  --fdr-cutoff 0.10 \
  --top-n 50
```

## Troubleshooting

### No comparison generated

**Possible causes:**
1. Fewer than 2 groups completed successfully
2. `compare_batch_results.R` not in same directory as `batch_timecourse.sh`
3. R packages not installed

**Solutions:**
```bash
# Check completion
ls batch_results/  # Should see multiple group directories

# Run manually
Rscript compare_batch_results.R --input-dir batch_results

# Install missing packages
R -e "install.packages(c('VennDiagram', 'UpSetR'))"
```

### Empty or missing plots

**Possible causes:**
1. No overlapping genes between groups
2. No enrichment analysis was run
3. Too few significant genes

**Solutions:**
```bash
# Check significant gene counts
cat batch_results/cross_group_comparison/tables/significant_genes_per_group.txt

# Relax threshold
Rscript compare_batch_results.R \
  --input-dir batch_results \
  --fdr-cutoff 0.10
```

### Heatmap too large/small

The heatmap height auto-adjusts based on number of genes/pathways. For better display:

```bash
# Show fewer pathways
Rscript compare_batch_results.R \
  --input-dir batch_results \
  --top-n 15  # Instead of default 20
```

## Example Workflows

### Workflow 1: Quick Overview

```bash
# 1. Run batch analysis
./batch_timecourse.sh --groups wt,ASC_KO,BST1_KO,GSDMD_KO ...

# 2. Check summary
cat batch_results/cross_group_comparison/COMPARISON_SUMMARY.txt

# 3. View main figures
open batch_results/cross_group_comparison/figures/*.pdf

# 4. Get gene lists
head batch_results/cross_group_comparison/tables/genes_shared_in_all_groups.txt
```

### Workflow 2: Detailed Analysis

```bash
# 1. Run batch with enrichment
./batch_timecourse.sh \
  --groups wt,ASC_KO,BST1_KO,GSDMD_KO \
  --cluster-genes \
  --run-enrichment \
  ...

# 2. Review comparison
cat batch_results/cross_group_comparison/COMPARISON_SUMMARY.txt

# 3. Extract specific findings
grep "Up in all" batch_results/cross_group_comparison/tables/direction_of_change_summary.txt

# 4. Find universal pathways
awk '$2 == 4' batch_results/cross_group_comparison/tables/pathway_comparison.txt

# 5. Identify KO-specific genes
wc -l batch_results/cross_group_comparison/tables/genes_unique_to_*.txt
```

### Workflow 3: Publication Figures

```bash
# 1. Run batch analysis
./batch_timecourse.sh --groups ... --run-enrichment ...

# 2. Generate high-quality figures
# Main figures are in:
batch_results/cross_group_comparison/figures/

# Recommended for publication:
# - gene_overlap_upset.pdf (multi-group overlap)
# - shared_genes_lfc_heatmap.pdf (expression patterns)
# - pathway_comparison_heatmap.pdf (biological themes)

# 3. Supplement with individual group results
# Use figures from each group directory
```

## Summary Statistics

The comparison analysis provides:

### Quantitative Metrics
- **N significant genes per group**
- **N shared genes (all groups)**
- **N unique genes (per group)**
- **% concordant direction**
- **N universal pathways**
- **N group-specific pathways**

### Qualitative Insights
- **Biological themes** (from pathway analysis)
- **Regulatory patterns** (from direction analysis)
- **Group relationships** (from clustering)

### Actionable Outputs
- **Gene lists** for validation
- **Pathway lists** for interpretation
- **Visual summaries** for presentation

## Best Practices

1. **Run enrichment** - Pathway comparison requires `--run-enrichment`
2. **Use consistent methods** - Same method for all groups
3. **Check individual results first** - Validate quality before comparison
4. **Consider biology** - Numbers are just numbers; interpret in context
5. **Follow up** - Use comparison to guide next experiments

## Next Steps

After reviewing cross-group comparisons:

1. **Pick candidates** - Select top shared/unique genes
2. **Validate** - qPCR, Western blots on selected genes
3. **Functional studies** - Test top pathways experimentally
4. **Systems biology** - Network analysis on shared genes
5. **Publish** - Use figures in your manuscript!

---

**The comparison analysis helps you move from individual group results to integrated biological understanding!** 🎯
