# Breakpoint Annotation & Fusion Prediction Restructuring

## Overview

Major restructuring of Phase 3 (Breakpoint Annotation) and Phase 4 (Fusion Prediction) to reduce output complexity and provide better decision-making data for fusion candidates.

---

## Phase 3: Breakpoint Annotation Restructuring

### Changes

#### New Output Format: Gene-Level Summary

**Old format**: Multiple rows per breakpoint (one per transcript)
- Created transcript-level explosion of permutations
- Complex to analyze for fusion prediction

**New format**: One row per breakpoint-gene pair with summary counts

Output columns:
- `variant_id`: Variant identifier
- `seqname`, `pos`, `breakpoint`: Genomic coordinates
- `variant_type`: Type of structural variant (from input)
- `is_gene`: TRUE if breakpoint overlaps any gene
- `gene_id`, `gene_symbol`, `gene_strand`: Gene annotation
- `is_exon`: TRUE if overlaps exon(s)
- `num_exon`: Count of exons overlapping this breakpoint for this gene
- `is_intron`: TRUE if overlaps intron(s)
- `num_intron`: Count of introns overlapping for this gene
- `is_cds`: TRUE if overlaps CDS
- `is_disease_gene`: Disease gene flag

#### Key Improvements

1. **Reduced complexity**: Summarizes multiple transcripts into gene-level counts
2. **Preserved intergenic breakpoints**: Intergenic positions included with is_gene=FALSE
3. **Clear annotation counts**: num_exon and num_intron show transcript overlap counts
4. **Better decision input**: Annotations now ideal for fusion prediction logic

#### Implementation Details

**New function**: `annotate_breakpoint_set()` (revised)
- Processes breakpoints per gene
- Counts exons/introns/CDS overlaps
- Returns one row per breakpoint-gene pair
- Includes variant_type from input or sv_variants parameter

**Helper functions**:
- `split_breakpoints_for_annotation()`: Splits break1/break2 into separate sets
- `convert_df_to_granges()`: No changes - works with renamed columns
- `get_gene_info()`: Used for gene lookups
- `find_nearest_gene()`: For future enhancement (currently not used)

---

## Phase 4: Fusion Prediction Restructuring

### Changes

#### 1. Canonical Fusion Decision Matrix

**New function**: `create_canonical_decision_matrix()`

Defines canonical fusions based on:
- `variant_type`: inversion, duplication, deletion, translocation, translocation_inverted
- `strand_bp1`, `strand_bp2`: Strand orientation at each breakpoint
- `is_canonical`: Boolean classification
- `promoter_breakpoint`: Which breakpoint provides promoter

**Logic**:
```
Inversion:   +/- → TRUE (different strands needed for canonical)
Duplication: +/+ → TRUE, -/- → TRUE (same strand, specific direction)
Deletion:    +/+ → TRUE, -/- → TRUE (same strand, specific direction)
Translocation: multiple valid strand combinations
Translocation_inverted: +/- or -/+ → TRUE
```

#### 2. New Fusion Prediction Logic

**Function**: `predict_fusions()` (complete rewrite)

Input: Gene-level breakpoint annotations

Process:
1. Get unique variants by variant_id
2. Extract bp1 and bp2 annotations
3. Generate ALL permutations of bp1_genes × bp2_genes
4. For each permutation:
   - Create fusion_id with format: `GENE1::GENE2_region1_region2_variant_id`
     - region = "e" (exon), "i" (intron), "0" (intergenic)
   - Look up canonicity in decision matrix using variant_type and strands
   - Return row with canonicity classification

Output columns:
- `fusion_id`: Describes both genes and their disruption regions
- `variant_type`: Structural variant type
- `bp1_gene_id`, `bp1_gene_symbol`: Gene at breakpoint 1
- `bp2_gene_id`, `bp2_gene_symbol`: Gene at breakpoint 2
- `bp1_gene_strand`, `bp2_gene_strand`: Strand information
- `is_canonical`: Classification from decision matrix
- `canonical_fusion_gene_promoter_breakpoint`: Which bp provides promoter (if canonical)

#### 3. Reading Frame Prediction (Revised)

**Function**: `predict_reading_frame()` (major updates)

Improvements:
1. **Strand information**: Includes strand in bp1_range and bp2_range creation
2. **Promoter/3'end logic**: Correctly identifies which breakpoint provides promoter vs 3' end
3. **MANE_Select preference**: New helper function `select_mane_cds()`
   - Prefers MANE_Select transcripts for more robust predictions
   - Falls back to first available CDS if no MANE_Select

New parameters:
- `bp_promoter_annot`: Breakpoint providing promoter (not bp1_data)
- `bp_3end_annot`: Breakpoint providing 3' end (not bp2_data)

Logic flow:
- Find downstream CDS from promoter breakpoint using `GenomicRanges::precede()`
- Find upstream CDS from 3' end breakpoint using `GenomicRanges::follow()`
- Compare CDS phase information for in-frame determination

#### 4. Removed Functions

Old functions no longer needed:
- `pair_breakpoints()`: Now implicit in variant grouping
- `classify_fusion_type()`: Replaced by canonical decision matrix
- `get_fusion_id()`: Incorporated into permutation generation
- `get_strand_info_with_sv_type()`: Not needed - strand info in annotations
- `get_strand_info()`: Replaced by direct strand columns
- `check_disease_genes()`: Can be added back if needed

---

## Updated Pipeline Flow

```
Input: multi-breakpoint variants
  |
  V
Phase 3: annotate_breakpoints()
  |
  +--split_breakpoints_for_annotation()
  |
  +--annotate_breakpoint_set() [for each bp set]
  |   |
  |   +--convert_df_to_granges()
  |   +--Query exons/introns/CDS for overlaps
  |   +--Count overlaps per gene
  |   +--return gene-level summary
  |
  V
Output: breakpoint_annotations (gene-level summary)
  - One row per breakpoint-gene pair
  - Includes num_exon, num_intron, is_cds flags
  - Preserves variant_type and intergenic breakpoints
  |
  V
Phase 4: predict_fusions()
  |
  +--create_canonical_decision_matrix()
  |
  +--For each variant:
  |   |
  |   +--Get all bp1 and bp2 annotations
  |   |
  |   +--Generate permutations (bp1_gene × bp2_gene)
  |   |
  |   +--For each permutation:
  |   |   |
  |   |   +--Create fusion_id
  |   |   +--Classify canonicity via decision matrix
  |   |   +--Return fusion result row
  |
  V
Output: fusions (fusion candidates)
  - All possible gene combinations per variant
  - Clear canonicity classification
  - Promoter breakpoint identified
  - Ready for reading frame analysis or further filtering
```

---

## Example Transformation

### Input (variants_processed)
```
variant_id=30EFE4, break1_chromosome=chr8, break1_position=38425101, 
break1_orientation=+, break2_chromosome=chr8, break2_position=84860033, 
break2_orientation=-, variant_type=inversion
```

### Phase 3 Output (breakpoint_annotations)
```
30EFE4, chr8, 38425101, 1, inversion, TRUE,  ENSG00000255201.1, ENSG00000255201, +, 0, 1
30EFE4, chr8, 38425101, 1, inversion, TRUE,  ENSG00000077782.24, FGFR1, -, 0, 30
30EFE4, chr8, 84860033, 2, inversion, TRUE,  ENSG00000184672.12, RALYL, +, 0, 7
```

### Phase 4 Output (fusions) - permutations
```
ENSG00000255201::RALYL_i_i_30EFE4, inversion, ..., FALSE, NA
FGFR1::RALYL_i_i_30EFE4, inversion, ..., TRUE, 1
RALYL::FGFR1_i_i_30EFE4, inversion, ..., TRUE, 2
... (all 6 permutations)
```

---

## Benefits of Restructuring

1. **Cleaner data**: Gene-level summary avoids transcript explosion
2. **Better decision-making**: Permutations with strand info drive canonicity
3. **Transparent logic**: Decision matrix explicitly shows valid combinations
4. **Extensible**: Can add filters/scoring at fusion prediction stage
5. **Accurate reading frame**: Strand info and MANE selection improve predictions

---

## Files Modified

- `R/breakpoint_annotation.R`: Complete restructure of annotate_breakpoint_set()
- `R/fusion_prediction.R`: New predict_fusions() with decision matrix, revised predict_reading_frame()
- `fusion_pipeline.qmd`: Updated function calls and logging for new output formats

---

## Migration Guide (for existing code)

If you have code expecting old breakpoint_annotations format:

**Old expectations**: transcript_id, region_type, region_details columns
**New format**: is_exon, num_exon, is_intron, num_intron, is_cds flags + counts

**Old fusion output**: fusion_type (canonical/non-canonical), transcript1_id, transcript2_id, in_frame
**New format**: is_canonical (boolean), canonical_fusion_gene_promoter_breakpoint, permutations of genes

Would need to be updated to work with new output structure.
