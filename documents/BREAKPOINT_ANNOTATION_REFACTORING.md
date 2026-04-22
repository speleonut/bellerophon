# Breakpoint Annotation Refactoring - Summary

## Problem Statement
The original `annotate_breakpoints()` function only processed a single breakpoint per variant, but the input data contains TWO breakpoints per variant (break1 and break2 positions). The input format has columns like:
- break1_chromosome, break1_position, break1_orientation
- break2_chromosome, break2_position, break2_orientation

## Solution Implemented

### 1. New Function: `split_breakpoints_for_annotation()`
**Purpose**: Convert multi-breakpoint input into separate breakpoint sets

**Input**:
- Tibble with columns: variant_id, break1_chromosome, break1_position, break1_orientation, break2_chromosome, break2_position, break2_orientation

**Output**:
- List with two elements:
  - `bp1`: Tibble with variant_id, seqname, pos, strand, breakpoint=1
  - `bp2`: Tibble with variant_id, seqname, pos, strand, breakpoint=2

**Key Details**:
- Renames break1_chromosome → seqname, break1_position → pos, break1_orientation → strand
- Renames break2_chromosome → seqname, break2_position → pos, break2_orientation → strand
- Adds breakpoint column (1 for bp1, 2 for bp2) to track which breakpoint each row represents
- Uses `dplyr::select()` and `dplyr::mutate()` for clean data transformation

### 2. New Function: `annotate_breakpoint_set()`
**Purpose**: Extract and encapsulate the core annotation logic for a single breakpoint set

**Input**:
- breakpoint_variants: Data frame with columns variant_id, seqname, pos, strand, breakpoint
- annotations: List from load_gencode_annotations()
- disease_genes: Optional disease gene tibble

**Output**:
- Tibble with annotated breakpoints:
  - variant_id, breakpoint, seqname, pos
  - gene_id, gene_symbol, transcript_id
  - region_type, region_details, is_disease_gene, distance_to_nearest_gene

**Key Details**:
- Contains all original annotation logic (exon overlaps, intron overlaps, nearest gene)
- Automatically extracts breakpoint number from input data
- Returns results with proper breakpoint numbering preserved

### 3. Modified Function: `annotate_breakpoints()`
**Changes**:
1. Added format detection:
   - Checks if input has "break1_chromosome" and "break2_chromosome" columns
   - Distinguishes multi-breakpoint format from single-breakpoint format

2. Added orchestration logic:
   - For multi-breakpoint input:
     - Calls `split_breakpoints_for_annotation()` to split variants
     - Loops through breakpoints 1 and 2
     - Calls `annotate_breakpoint_set()` for each breakpoint
     - Combines results using `dplyr::bind_rows()`
   - For single-breakpoint format:
     - Uses original logic via `annotate_breakpoint_set()`

3. Improved documentation:
   - Updated @param variants to document both input formats
   - Clarified that function now handles dual breakpoints

## Data Flow

```
Input Variants (multi-breakpoint format)
  |
  V
split_breakpoints_for_annotation()
  |
  +---> variants_bp1 (break1 data, breakpoint=1)
  |       |
  |       V
  |     annotate_breakpoint_set() --> results_bp1
  |
  +---> variants_bp2 (break2 data, breakpoint=2)
          |
          V
        annotate_breakpoint_set() --> results_bp2
  |
  V
dplyr::bind_rows(results_bp1, results_bp2)
  |
  V
Output: Combined annotations for all breakpoints
```

## Example Data Transformation

### Input (Single variant with dual breakpoints)
```
variant_id: 30EFE4
break1_chromosome: 8
break1_position: 38425101
break1_orientation: +
break2_chromosome: 8
break2_position: 84860033
break2_orientation: -
```

### After split_breakpoints_for_annotation()
**variants_bp1:**
```
variant_id: 30EFE4
seqname: 8
pos: 38425101
strand: +
breakpoint: 1
```

**variants_bp2:**
```
variant_id: 30EFE4
seqname: 8
pos: 84860033
strand: -
breakpoint: 2
```

### After annotation (combined output)
```
Row 1: variant_id=30EFE4, breakpoint=1, seqname=8, pos=38425101, gene_id=ENSG..., ...
Row 2: variant_id=30EFE4, breakpoint=2, seqname=8, pos=84860033, gene_id=ENSG..., ...
```

## Backward Compatibility

- Function maintains backward compatibility with single-breakpoint input formats
- Existing code calling `annotate_breakpoints()` will continue to work
- Output structure remains unchanged (same columns in results)

## Column Name Mapping

The refactoring uses consistent column naming:
- break1_chromosome → seqname (compatible with GRanges)
- break1_position → pos (compatible with GRanges)
- break1_orientation → strand (compatible with GRanges)
- Prefixes "break1_" and "break2_" allow function to detect multi-breakpoint format

## Testing

Created test_breakpoint_fix.R with validation tests for:
1. Correct splitting of breakpoints
2. Column naming consistency
3. Breakpoint numbering (1 and 2)
4. Data integrity preservation

## Files Modified

- R/breakpoint_annotation.R: Added split_breakpoints_for_annotation(), annotate_breakpoint_set(), refactored annotate_breakpoints()
- test_breakpoint_fix.R: Created validation test script (NEW)
