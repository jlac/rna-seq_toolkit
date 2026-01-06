#!/bin/bash

# Batch Timecourse Analysis Wrapper
# Runs multiple timecourse analyses or interaction tests from a single command

# Note: We don't use 'set -e' because we want to continue processing
# remaining groups even if one fails. Error handling is explicit.

# Default values
SCRIPT="timecourse_analysis_cli.R"
COUNTS=""
DESIGN=""
TIME_COLUMN="Time"
GROUP_COLUMN="Treatment"
OUTPUT_BASE="batch_results"
METHOD="both"
TIME_NUMERIC=""
SPLINE_DF="3"
MIN_COUNT="10"
MIN_SAMPLES="3"
CLUSTER=""
ENRICHMENT=""
VERBOSE=""

# Parse command line arguments
show_help() {
    cat << EOF
Batch Timecourse Analysis Wrapper

Usage: $0 [OPTIONS]

Required:
  --counts FILE               Path to counts file
  --design FILE               Path to design file
  
Batch Options (choose one):
  --groups GROUP1,GROUP2,...  Run separate timecourse for each group
  --interactions PAIR1,PAIR2  Run interaction analyses for pairs (format: wt-ASC_KO,wt-BST1_KO)
  
Optional:
  --script PATH               Path to timecourse_analysis_cli.R [default: $SCRIPT]
  --time-column NAME          Time column name [default: $TIME_COLUMN]
  --time-numeric              Treat time as numeric (use splines)
  --group-column NAME         Group column name [default: $GROUP_COLUMN]
  -o, --output DIR            Base output directory [default: $OUTPUT_BASE]
  --method STRING             Analysis method: limma, deseq2, or both [default: $METHOD]
  --spline-df INT             Spline degrees of freedom [default: $SPLINE_DF]
  --min-count INT             Minimum count threshold [default: $MIN_COUNT]
  --min-samples INT           Minimum samples threshold [default: $MIN_SAMPLES]
  --top-genes INT             Number of top genes for plots [default: 50]
  --duplicate-genes STRING    How to handle duplicate gene IDs: unique, sum, or first [default: unique]
  --cluster-genes             Enable gene clustering
  --n-clusters INT            Number of clusters [default: 6]
  --run-enrichment            Enable pathway enrichment
  --enrichment-db STRING      Enrichment database: reactome, gobp, or gmt [default: reactome]
  --organism STRING           Organism: human or mouse [default: human]
  --verbose                   Verbose output
  --help                      Show this help message

Examples:
  # Run separate analyses for multiple groups
  $0 --counts counts.txt --design design.txt \\
     --groups wt,ASC_KO,BST1_KO,GSDMD_KO

  # Run multiple interaction analyses  
  $0 --counts counts.txt --design design.txt \\
     --interactions "wt-ASC_KO,wt-BST1_KO,wt-GSDMD_KO"

  # With clustering and enrichment
  $0 --counts counts.txt --design design.txt \\
     --groups wt,ASC_KO,BST1_KO \\
     --cluster-genes --n-clusters 6 \\
     --run-enrichment --enrichment-db reactome --organism human

EOF
    exit 0
}

# Parse arguments
BATCH_GROUPS=""
BATCH_INTERACTIONS=""
N_CLUSTERS="6"
ENRICHMENT_DB="reactome"
ORGANISM="human"
TOP_GENES="50"
DUPLICATE_GENES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --help) show_help ;;
        --counts) COUNTS="$2"; shift 2 ;;
        --design) DESIGN="$2"; shift 2 ;;
        --script) SCRIPT="$2"; shift 2 ;;
        --groups) BATCH_GROUPS="$2"; shift 2 ;;
        --interactions) BATCH_INTERACTIONS="$2"; shift 2 ;;
        --time-column) TIME_COLUMN="$2"; shift 2 ;;
        --time-numeric) TIME_NUMERIC="--time-numeric"; shift ;;
        --group-column) GROUP_COLUMN="$2"; shift 2 ;;
        --output|-o) OUTPUT_BASE="$2"; shift 2 ;;
        --method) METHOD="$2"; shift 2 ;;
        --spline-df) SPLINE_DF="$2"; shift 2 ;;
        --min-count) MIN_COUNT="$2"; shift 2 ;;
        --min-samples) MIN_SAMPLES="$2"; shift 2 ;;
        --top-genes) TOP_GENES="$2"; shift 2 ;;
        --duplicate-genes) DUPLICATE_GENES="$2"; shift 2 ;;
        --cluster-genes) CLUSTER="--cluster-genes"; shift ;;
        --n-clusters) N_CLUSTERS="$2"; shift 2 ;;
        --run-enrichment) ENRICHMENT="--run-enrichment"; shift ;;
        --enrichment-db) ENRICHMENT_DB="$2"; shift 2 ;;
        --organism) ORGANISM="$2"; shift 2 ;;
        --verbose) VERBOSE="--verbose"; shift ;;
        *) echo "Unknown option: $1"; show_help ;;
    esac
done

# Validate required arguments
if [[ -z "$COUNTS" ]] || [[ -z "$DESIGN" ]]; then
    echo "Error: --counts and --design are required"
    show_help
fi

if [[ -z "$BATCH_GROUPS" ]] && [[ -z "$BATCH_INTERACTIONS" ]]; then
    echo "Error: Either --groups or --interactions must be specified"
    show_help
fi

if [[ ! -f "$SCRIPT" ]]; then
    echo "Error: Script not found: $SCRIPT"
    exit 1
fi

if [[ ! -f "$COUNTS" ]]; then
    echo "Error: Counts file not found: $COUNTS"
    exit 1
fi

if [[ ! -f "$DESIGN" ]]; then
    echo "Error: Design file not found: $DESIGN"
    exit 1
fi

# Create base output directory
mkdir -p "$OUTPUT_BASE"

# Build common options
COMMON_OPTS="--counts $COUNTS --design $DESIGN --time-column $TIME_COLUMN --method $METHOD $TIME_NUMERIC --spline-df $SPLINE_DF --min-count $MIN_COUNT --min-samples $MIN_SAMPLES --top-genes $TOP_GENES $VERBOSE"

# Add duplicate genes handling if specified
if [[ -n "$DUPLICATE_GENES" ]]; then
    COMMON_OPTS="$COMMON_OPTS --duplicate-genes $DUPLICATE_GENES"
fi

if [[ -n "$CLUSTER" ]]; then
    COMMON_OPTS="$COMMON_OPTS $CLUSTER --n-clusters $N_CLUSTERS"
fi

if [[ -n "$ENRICHMENT" ]]; then
    COMMON_OPTS="$COMMON_OPTS $ENRICHMENT --enrichment-db $ENRICHMENT_DB --organism $ORGANISM"
fi

# Function to run a single group analysis
run_group() {
    local group=$1
    local output_dir="${OUTPUT_BASE}/${group}"
    
    echo ""
    echo "================================================================================"
    echo "Analyzing group: $group"
    echo "================================================================================"
    echo ""
    
    # Create group-specific design file
    local temp_design="${output_dir}_design.txt"
    mkdir -p "$output_dir"
    
    # Filter design file to this group
    head -1 "$DESIGN" > "$temp_design"
    awk -v group="$group" -v col="$GROUP_COLUMN" '
        NR==1 {
            for(i=1; i<=NF; i++) {
                if($i==col) {
                    col_idx=i
                    break
                }
            }
            next
        }
        col_idx > 0 && $col_idx==group { print }
    ' "$DESIGN" >> "$temp_design"
    
    # Check if any samples were found
    local n_samples=$(tail -n +2 "$temp_design" | wc -l)
    if [[ $n_samples -eq 0 ]]; then
        echo "  Warning: No samples found for group '$group'. Skipping."
        rm "$temp_design"
        return 1
    fi
    
    echo "  Found $n_samples samples in group '$group'"
    
    # Build command with filtered design file
    local CMD="$SCRIPT --counts $COUNTS --design $temp_design --time-column $TIME_COLUMN --method $METHOD $TIME_NUMERIC --spline-df $SPLINE_DF --min-count $MIN_COUNT --min-samples $MIN_SAMPLES --top-genes $TOP_GENES $VERBOSE"
    
    # Add duplicate genes handling if specified
    if [[ -n "$DUPLICATE_GENES" ]]; then
        CMD="$CMD --duplicate-genes $DUPLICATE_GENES"
    fi
    
    if [[ -n "$CLUSTER" ]]; then
        CMD="$CMD $CLUSTER --n-clusters $N_CLUSTERS"
    fi
    
    if [[ -n "$ENRICHMENT" ]]; then
        CMD="$CMD $ENRICHMENT --enrichment-db $ENRICHMENT_DB --organism $ORGANISM"
    fi
    
    CMD="$CMD --output-dir $output_dir"
    
    # Run analysis
    Rscript $CMD || {
        echo "Error analyzing $group"
        rm "$temp_design"
        return 1
    }
    
    # Clean up temp file
    rm "$temp_design"
    
    echo ""
    echo "✓ Completed: $group"
    echo ""
}

# Function to run an interaction analysis
run_interaction() {
    local group1=$1
    local group2=$2
    local comparison="${group1}_vs_${group2}"
    local output_dir="${OUTPUT_BASE}/${comparison}"
    
    echo ""
    echo "================================================================================"
    echo "Analyzing interaction: $group1 vs $group2"
    echo "================================================================================"
    echo ""
    
    mkdir -p "$output_dir"
    
    # Create filtered design file with only these two groups
    local temp_design="${output_dir}_design.txt"
    head -1 "$DESIGN" > "$temp_design"
    awk -v g1="$group1" -v g2="$group2" -v col="$GROUP_COLUMN" '
        NR==1 {
            for(i=1; i<=NF; i++) {
                if($i==col) {
                    col_idx=i
                    break
                }
            }
            next
        }
        col_idx > 0 && ($col_idx==g1 || $col_idx==g2) { print }
    ' "$DESIGN" >> "$temp_design"
    
    # Check if any samples were found
    local n_samples=$(tail -n +2 "$temp_design" | wc -l)
    if [[ $n_samples -eq 0 ]]; then
        echo "  Warning: No samples found for groups '$group1' and '$group2'. Skipping."
        rm "$temp_design"
        return 1
    fi
    
    echo "  Found $n_samples total samples for comparison"
    
    # Build command with filtered design file
    local CMD="$SCRIPT --counts $COUNTS --design $temp_design --time-column $TIME_COLUMN --method $METHOD $TIME_NUMERIC --spline-df $SPLINE_DF --min-count $MIN_COUNT --min-samples $MIN_SAMPLES --top-genes $TOP_GENES $VERBOSE"
    
    # Add duplicate genes handling if specified
    if [[ -n "$DUPLICATE_GENES" ]]; then
        CMD="$CMD --duplicate-genes $DUPLICATE_GENES"
    fi
    
    # Add interaction-specific options
    CMD="$CMD --group-column $GROUP_COLUMN --groups ${group1},${group2} --test-interaction"
    
    if [[ -n "$CLUSTER" ]]; then
        CMD="$CMD $CLUSTER --n-clusters $N_CLUSTERS"
    fi
    
    if [[ -n "$ENRICHMENT" ]]; then
        CMD="$CMD $ENRICHMENT --enrichment-db $ENRICHMENT_DB --organism $ORGANISM"
    fi
    
    CMD="$CMD --output-dir $output_dir"
    
    # Run analysis
    Rscript $CMD || {
        echo "Error analyzing $group1 vs $group2"
        rm "$temp_design"
        return 1
    }
    
    # Clean up temp file
    rm "$temp_design"
    
    echo ""
    echo "✓ Completed: $group1 vs $group2"
    echo ""
}

# Main execution
START_TIME=$(date +%s)

echo ""
echo "================================================================================"
echo "BATCH TIMECOURSE ANALYSIS"
echo "================================================================================"
echo ""
echo "Configuration:"
echo "  Counts file: $COUNTS"
echo "  Design file: $DESIGN"
echo "  Time column: $TIME_COLUMN"
echo "  Group column: $GROUP_COLUMN"
echo "  Method: $METHOD"
echo "  Output directory: $OUTPUT_BASE"
echo ""

if [[ -n "$BATCH_GROUPS" ]]; then
    # Run multiple single-group analyses
    IFS=',' read -ra GROUP_ARRAY <<< "$BATCH_GROUPS"
    
    # Trim whitespace from each group name
    for i in "${!GROUP_ARRAY[@]}"; do
        GROUP_ARRAY[$i]=$(echo "${GROUP_ARRAY[$i]}" | tr -d ' ')
    done
    
    echo "Will analyze ${#GROUP_ARRAY[@]} groups:"
    for group in "${GROUP_ARRAY[@]}"; do
        echo "  - $group"
    done
    echo ""
    
    FAILED=0
    COMPLETED=0
    
    for group in "${GROUP_ARRAY[@]}"; do
        if run_group "$group"; then
            ((COMPLETED++))
        else
            ((FAILED++))
        fi
    done
    
    echo ""
    echo "================================================================================"
    echo "BATCH ANALYSIS COMPLETE"
    echo "================================================================================"
    echo "  Total groups: ${#GROUP_ARRAY[@]}"
    echo "  Completed: $COMPLETED"
    echo "  Failed: $FAILED"
    
    # Run cross-group comparison if multiple groups completed successfully
    if [[ $COMPLETED -ge 2 ]]; then
        echo ""
        echo "================================================================================"
        echo "RUNNING CROSS-GROUP COMPARISON"
        echo "================================================================================"
        echo ""
        
        COMPARE_SCRIPT="$(dirname "$SCRIPT")/compare_batch_results.R"
        if [[ ! -f "$COMPARE_SCRIPT" ]]; then
            COMPARE_SCRIPT="compare_batch_results.R"
        fi
        
        if [[ -f "$COMPARE_SCRIPT" ]]; then
            Rscript "$COMPARE_SCRIPT" \
                --input-dir "$OUTPUT_BASE" \
                --output-dir "${OUTPUT_BASE}/cross_group_comparison" \
                --method "$METHOD" \
                --organism "$ORGANISM" \
                $VERBOSE || {
                echo "Warning: Cross-group comparison failed"
            }
            
            if [[ -f "${OUTPUT_BASE}/cross_group_comparison/COMPARISON_SUMMARY.txt" ]]; then
                echo ""
                echo "✓ Cross-group comparison complete!"
                echo "  Results saved to: ${OUTPUT_BASE}/cross_group_comparison"
            fi
        else
            echo "Note: compare_batch_results.R not found - skipping cross-group comparison"
            echo "      Download it separately to enable this feature"
        fi
    fi
    
elif [[ -n "$BATCH_INTERACTIONS" ]]; then
    # Run multiple interaction analyses
    IFS=',' read -ra PAIR_ARRAY <<< "$BATCH_INTERACTIONS"
    
    # Trim whitespace from each pair
    for i in "${!PAIR_ARRAY[@]}"; do
        PAIR_ARRAY[$i]=$(echo "${PAIR_ARRAY[$i]}" | tr -d ' ')
    done
    
    echo "Will analyze ${#PAIR_ARRAY[@]} interaction pairs:"
    for pair in "${PAIR_ARRAY[@]}"; do
        IFS='-' read -ra GROUPS_IN_PAIR <<< "$pair"
        # Trim group names
        GROUPS_IN_PAIR[0]=$(echo "${GROUPS_IN_PAIR[0]}" | tr -d ' ')
        GROUPS_IN_PAIR[1]=$(echo "${GROUPS_IN_PAIR[1]}" | tr -d ' ')
        echo "  - ${GROUPS_IN_PAIR[0]} vs ${GROUPS_IN_PAIR[1]}"
    done
    echo ""
    
    FAILED=0
    COMPLETED=0
    
    for pair in "${PAIR_ARRAY[@]}"; do
        IFS='-' read -ra GROUPS_IN_PAIR <<< "$pair"
        # Trim group names
        GROUPS_IN_PAIR[0]=$(echo "${GROUPS_IN_PAIR[0]}" | tr -d ' ')
        GROUPS_IN_PAIR[1]=$(echo "${GROUPS_IN_PAIR[1]}" | tr -d ' ')
        if [[ ${#GROUPS_IN_PAIR[@]} -ne 2 ]]; then
            echo "Error: Invalid pair format: $pair (use GROUP1-GROUP2)"
            ((FAILED++))
            continue
        fi
        
        if run_interaction "${GROUPS_IN_PAIR[0]}" "${GROUPS_IN_PAIR[1]}"; then
            ((COMPLETED++))
        else
            ((FAILED++))
        fi
    done
    
    echo ""
    echo "================================================================================"
    echo "BATCH ANALYSIS COMPLETE"
    echo "================================================================================"
    echo "  Total pairs: ${#PAIR_ARRAY[@]}"
    echo "  Completed: $COMPLETED"
    echo "  Failed: $FAILED"
    
    # Run cross-group comparison if multiple interaction analyses completed successfully
    if [[ $COMPLETED -ge 2 ]]; then
        echo ""
        echo "================================================================================"
        echo "RUNNING CROSS-GROUP COMPARISON"
        echo "================================================================================"
        echo ""
        
        COMPARE_SCRIPT="$(dirname "$SCRIPT")/compare_batch_results.R"
        if [[ ! -f "$COMPARE_SCRIPT" ]]; then
            COMPARE_SCRIPT="compare_batch_results.R"
        fi
        
        if [[ -f "$COMPARE_SCRIPT" ]]; then
            Rscript "$COMPARE_SCRIPT" \
                --input-dir "$OUTPUT_BASE" \
                --output-dir "${OUTPUT_BASE}/cross_group_comparison" \
                --method "$METHOD" \
                --organism "$ORGANISM" \
                $VERBOSE || {
                echo "Warning: Cross-group comparison failed"
            }
            
            if [[ -f "${OUTPUT_BASE}/cross_group_comparison/COMPARISON_SUMMARY.txt" ]]; then
                echo ""
                echo "✓ Cross-group comparison complete!"
                echo "  Results saved to: ${OUTPUT_BASE}/cross_group_comparison"
            fi
        else
            echo "Note: compare_batch_results.R not found - skipping cross-group comparison"
            echo "      Download it separately to enable this feature"
        fi
    fi
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "  Time elapsed: ${ELAPSED}s"
echo ""
echo "Results saved to: $OUTPUT_BASE"
echo ""
