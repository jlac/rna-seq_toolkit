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
