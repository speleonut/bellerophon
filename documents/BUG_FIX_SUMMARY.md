# Bug Fix Summary: Fusion Prediction and Breakpoint Annotation

## Bugs Fixed

### Bug 1: Missed Reciprocal Fusions for Inversions and Translocations

**Issue**: For inversions and translocation_inverted variants, only one fusion gene permutation (bp1::bp2) was being generated, missing the reciprocal fusion (bp2::bp1).

**Example**: For variant 30EFE4 (inversion):
- FGFR1 (ENSG00000077782.24) at breakpoint 1
- RALYL (ENSG00000184672.12) at breakpoint 2

**Missing output** (should be generated):
- RALYL::FGFR1_i_i_30EFE4 (reciprocal)
- RALYL::ENSG00000255201_i_i_30EFE4 (reciprocal)

**Why it matters**: Inversions and translocations can form TWO different fusion genes:
1. Gene A (promoter) → Gene B (coding sequence)
2. Gene B (reciprocal promoter) → Gene A (reciprocal coding sequence)

These have:
- Different promoter elements (different upstream regulation)
- Different 3' coding regions
- Different expression patterns
- Different protein products

**Solution**: Modified `predict_fusions()` in [R/fusion_prediction.R](R/fusion_prediction.R) (lines 155-175):

```r
# Determine if this variant type supports reciprocal fusions
# Inversions and translocation_inverted can form reciprocal fusions (bp2::bp1 in addition to bp1::bp2)
supports_reciprocal <- variant_type %in% c("inversion", "translocation_inverted")

# Process both permutations (bp1::bp2 and bp2::bp1) if supports_reciprocal
# Otherwise only process bp1::bp2
permutations_to_process <- list(list(bp1 = bp1_annot[i, ], bp2 = bp2_annot[j, ]))

if (supports_reciprocal) {
  # Add reciprocal permutation: bp2::bp1
  permutations_to_process[[2]] <- list(bp1 = bp2_annot[j, ], bp2 = bp1_annot[i, ])
}

for (perm_idx in seq_along(permutations_to_process)) {
  # Generate fusion for each permutation
}
```

**Result**: Now generates BOTH permutations for inversions/translocation_inverted, while single-breakpoint variants (duplication, deletion, translocation) continue to generate only bp1::bp2.

---

### Bug 2: In-Frame Fusions Not Detected Due to CDS Overlap Check

**Issue**: The breakpoint annotation module required direct overlap with CDS regions to set `is_cds=TRUE`. However, many clinically important fusion breakpoints occur in **introns that are part of the CDS boundaries** of transcripts.

**Problem code** in [R/breakpoint_annotation.R](R/breakpoint_annotation.R) (old line 223):
```r
cds_overlaps <- IRanges::findOverlaps(variant, annotations$cds, ignore.strand = TRUE)
```

This only finds exact overlaps with CDS exon boundaries. A breakpoint in an intron between two exons won't overlap the CDS GRanges directly.

**Why it matters**: For in-frame reading frame calculation:
- `predict_reading_frame()` needs to find the CDS exons adjacent to the breakpoint
- Previously, introns with no direct CDS overlap were marked `is_cds=FALSE`
- This prevented reading frame calculation, losing critical fusion validation data

**Solution**: Two complementary changes:

1. **Remove direct CDS overlap check** (lines 220-226):
   ```r
   # Query for overlapping exons
   exon_overlaps <- IRanges::findOverlaps(variant, annotations$exons, ignore.strand = TRUE)
   
   # Query for overlapping introns
   intron_overlaps <- IRanges::findOverlaps(variant, annotations$introns, ignore.strand = TRUE)
   
   # Note: CDS overlap is NOT checked here because introns within CDS boundaries won't overlap directly.
   # Instead, predict_reading_frame() will independently retrieve all CDS for candidate genes via tx_lookup.
   ```

2. **Check gene existence via tx_lookup** (lines 279-290):
   ```r
   # Check if gene has CDS (will be determined by predict_reading_frame if needed)
   # Note: We check if the gene exists in any transcripts with CDS data
   has_cds <- FALSE
   if (!is.null(annotations$tx_lookup)) {
     gene_txs <- annotations$tx_lookup %>%
       dplyr::filter(gene_id == !!gene_id) %>%
       dplyr::pull(tx_id)
     # Assume CDS present if we have transcripts for this gene
     # predict_reading_frame will verify actual CDS overlap
     has_cds <- length(gene_txs) > 0
   }
   ```

**Result**: 
- Breakpoint annotation now focuses on exon/intron disruption (structural)
- Reading frame validation moved entirely to `predict_reading_frame()` which:
  - Independently retrieves all transcripts for each gene via tx_lookup
  - Uses `GenomicRanges::precede()` and `GenomicRanges::follow()` to find adjacent CDS
  - Handles introns correctly because it queries full transcript CDS boundaries
  - Determines phase status without requiring direct CDS overlap

---

## Code Changes Summary

### File 1: R/breakpoint_annotation.R

**Lines 220-226**: Remove `cds_overlaps` query, add explanatory comment

**Lines 279-290**: Replace hardcoded CDS overlap check with gene existence check via tx_lookup

### File 2: R/fusion_prediction.R

**Lines 155-175**: Add reciprocal fusion support for inversions/translocation_inverted

**Lines 166-280**: Wrap main fusion generation logic in permutation loop, ensures both forward and reciprocal fusions are processed with same canonicity/reading frame logic

---

## Testing Recommendations

### Test Case 1: Reciprocal Fusions
**Input**: Inversion variant 30EFE4
- Breakpoint 1: FGFR1 (-), ENSG00000255201 (+), both at chr8:38425101
- Breakpoint 2: RALYL (+) at chr8:84860032

**Expected Output** (all 4 fusions):
```
FGFR1::RALYL_i_i_30EFE4        (bp1::bp2)
ENSG00000255201::RALYL_i_i_30EFE4  (bp1::bp2)
RALYL::FGFR1_i_i_30EFE4        (bp2::bp1 reciprocal) ← NEW
RALYL::ENSG00000255201_i_i_30EFE4  (bp2::bp1 reciprocal) ← NEW
```

### Test Case 2: In-Frame Detection
**Scenario**: Duplication variant with intron breakpoints in both genes

**Before fix**:
- `is_cds=FALSE` for intron breakpoints → no reading frame calculation

**After fix**:
- `is_cds=TRUE` (gene has transcripts in tx_lookup) → `predict_reading_frame()` calculates frame → `in_frame=TRUE/FALSE` with valid `frame_offset`

### Test Case 3: Non-reciprocal Variants
**Verify**: Duplication/deletion variants still only generate bp1::bp2 fusions (not reciprocal)

Input: Duplication 4984E3
- SCN2A at bp1, SCN1A at bp2

Expected: Only `SCN2A::SCN1A` (NOT `SCN1A::SCN2A`)

---

## Edge Cases Handled

1. **Intron breakpoints in coding genes**: Now correctly identified for reading frame analysis
2. **Introns with both exon and intron features**: Handled by gene-level summary (one row per gene, all counts included)
3. **Intergenic breakpoints**: bp2 gene_id=NA → fusion won't have reciprocal (no strand to compare)
4. **Missing transcript data**: Returns NA from predict_reading_frame instead of crashing

---

## Backward Compatibility

- ✅ Output format unchanged (same columns)
- ✅ Existing non-inversion variants unaffected (single permutation only)
- ⚠️ Inversion and translocation_inverted fusions now have MORE entries (will expand result sets)
- ✅ `is_cds` flag now reflects "gene has CDS available" rather than "breakpoint overlaps CDS"

