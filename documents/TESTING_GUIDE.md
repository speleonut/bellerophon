# Testing Guide for Bug Fixes

## Quick Test: Run fusion_pipeline.qmd

Run the pipeline with your test data to verify both bugs are fixed:

```r
# In R/RStudio:
quarto::quarto_render("fusion_pipeline.qmd")
```

---

## Verification Tests

### Test 1: Reciprocal Fusions Generated ✅

**Data**: Use the example breakpoint annotations provided with inversions

**Expected Result**: For inversion 30EFE4:
- Original output: 2 fusions (FGFR1::RALYL, ENSG00000255201::RALYL)
- **New output: 4 fusions** (plus reciprocal: RALYL::FGFR1, RALYL::ENSG00000255201)

**Verification command**:
```r
fusions %>%
  dplyr::filter(variant_id == "30EFE4") %>%
  dplyr::count(variant_id)
# Should show: n = 4 (not 2)
```

---

### Test 2: In-Frame Fusions Detected ✅

**Data**: Variants with intron breakpoints within CDS boundaries

**Expected Result**: Reading frame should be calculated for eligible fusions:
- `in_frame` column has TRUE/FALSE (not all NA)
- `frame_offset` has valid values (0, 1, 2)

**Verification command**:
```r
fusions %>%
  dplyr::filter(is_canonical == TRUE) %>%
  dplyr::summarise(
    n_fusions = n(),
    n_inframe = sum(in_frame == TRUE, na.rm = TRUE),
    n_outframe = sum(in_frame == FALSE, na.rm = TRUE),
    n_unknown = sum(is.na(in_frame))
  )
# Should show: n_inframe > 0 and/or n_outframe > 0
```

---

### Test 3: Non-Reciprocal Variants Unchanged ✅

**Data**: Duplications and deletions should NOT get reciprocal fusions

**Expected Result**: For duplication 4984E3:
- SCN2A::SCN1A present
- **No** SCN1A::SCN2A (reciprocal)

**Verification command**:
```r
fusions %>%
  dplyr::filter(variant_id == "4984E3") %>%
  dplyr::select(fusion_id, variant_type, bp1_gene_symbol, bp2_gene_symbol)
# All fusions should have bp1=SCN2A, bp2=SCN1A (or intergenic)
# None should have bp1=SCN1A, bp2=SCN2A
```

---

### Test 4: is_cds Flag Now Indicates Gene Potential

**Data**: Check breakpoint_annotations output

**Expected Result**: 
- `is_cds` now = TRUE for ANY breakpoint in a gene with transcripts
- NOT limited to direct CDS overlap

**Verification command**:
```r
breakpoint_annotations %>%
  dplyr::filter(is_gene == TRUE) %>%
  dplyr::count(is_cds)
# Should show: most genes have is_cds=TRUE
# (instead of FALSE for intron breakpoints)
```

---

### Test 5: Complete End-to-End Validation

Run full fusion_pipeline.qmd and verify no errors:

```bash
# In terminal:
cd /path/to/bellerophon
Rscript example_pipeline.R --input data/example.vcf --output output/test1
```

**Check output**:
- ✅ Phase 3: Breakpoint annotation completes
- ✅ Phase 4: Fusion prediction generates expected permutations
- ✅ Phase 5: Expression hotspots calculated
- ✅ Phase 6: Report generates without errors
- ✅ Phase 7: Validation completes

---

## Expected Changes in Output

### Before Fix
```
Predicted 6 fusions
  Canonical: 2
  Non-canonical: 4
  In-frame: 0
  Out-of-frame: 0
```

### After Fix
```
Predicted 9 fusions (includes reciprocals)
  Canonical: 4
  Non-canonical: 5
  In-frame: 2 (or more)
  Out-of-frame: 2 (or more)
```

---

## Troubleshooting

### Issue: Still seeing old number of fusions
**Solution**: Ensure you're using updated `fusion_prediction.R` with reciprocal support

### Issue: `in_frame` still all NA
**Solution**: Verify `breakpoint_annotations$is_cds` is TRUE for genes (not just direct overlaps)

### Issue: Error about `promoter_bp` being NA
**Solution**: This is expected for non-canonical fusions. Check that `is_canonical` is FALSE when `promoter_bp` is NA.

---

## Files Modified

1. **R/breakpoint_annotation.R**
   - Removed direct CDS overlap check (line 223)
   - Changed `is_cds` logic to check gene existence in tx_lookup (lines 279-280)
   - Result: More genes marked as having CDS potential

2. **R/fusion_prediction.R**
   - Added reciprocal fusion support for inversions/translocation_inverted (lines 154-175)
   - Wrapped fusion generation in permutation loop (lines 166-257)
   - Result: Both forward and reciprocal permutations generated

3. **Documentation**: BUG_FIX_SUMMARY.md

---

## Next Steps After Verification

Once testing confirms both bugs are fixed:

1. ✅ Update README.md with note about reciprocal fusions in inversions
2. ✅ Consider adding filter option: `include_reciprocal = TRUE/FALSE` for users who only want canonical bridges
3. ✅ Update any downstream analysis scripts that assume 1 fusion per permutation
4. ✅ Consider performance impact: reciprocal fusions double the output for inversions

