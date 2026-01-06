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
