# Pipeline Compatibility Updates

## Overview
Updated `expression_analysis.R`, `report_generation.R`, and `fusion_pipeline.qmd` to ensure full compatibility with the restructured `fusion_prediction.R` and `breakpoint_annotation.R` modules.

---

## Column Name Changes

### fusion_prediction.R Output
The new `predict_fusions()` function outputs gene pairs using `bp1_` and `bp2_` prefixes:

**Old naming (pre-update):**
- `gene1_id`, `gene2_id`
- `gene1_symbol`, `gene2_symbol`

**New naming:**
- `bp1_gene_id`, `bp2_gene_id` (breakpoint 1/2 genes)
- `bp1_gene_symbol`, `bp2_gene_symbol`
- `bp1_gene_strand`, `bp2_gene_strand`

**Output columns (predict_fusions):**
```
fusion_id                    (character) - unique fusion identifier
variant_id                   (character) - structural variant ID [NEWLY ADDED]
variant_type                 (character) - type of SV (inversion, duplication, etc.)
bp1_gene_id                  (character) - ENSG ID for breakpoint 1 gene
bp1_gene_symbol              (character) - gene symbol for breakpoint 1
bp2_gene_id                  (character) - ENSG ID for breakpoint 2 gene
bp2_gene_symbol              (character) - gene symbol for breakpoint 2
bp1_gene_strand              (character) - strand of bp1 gene
bp2_gene_strand              (character) - strand of bp2 gene
is_canonical                 (logical)   - TRUE if meets canonical criteria from decision matrix
canonical_fusion_gene_promoter_breakpoint (integer) - which breakpoint provides promoter (1 or 2)
in_frame                     (logical)   - TRUE if reading frame is maintained
frame_offset                 (integer)   - CDS phase offset (0, 1, or 2)
```

---

## Module-by-Module Updates

### 1. expression_analysis.R

**Function:** `identify_fusion_hotspots()`

**Changes:**
- Updated to accept `bp1_gene_id`, `bp2_gene_id`, `bp1_gene_symbol`, `bp2_gene_symbol` from fusion_data
- Extracts these columns from the new fusion_prediction output format
- Internal output still uses `gene1_id`, `gene2_id`, `gene1_symbol`, `gene2_symbol` for compatibility with report generation

**Input expectations:**
```r
fusion_data <- predict_fusions(...)  # Contains bp1_gene_id, bp2_gene_id, etc.

expression_hotspots <- identify_fusion_hotspots(
  fusion_data = fusion_data,
  gtex_expression = gtex_expression
)
```

**Output columns:**
```
fusion_id             (character)
gene1_id              (character) - upstream gene (bp1)
gene1_symbol          (character)
gene2_id              (character) - downstream gene (bp2)
gene2_symbol          (character)
tissue                (character)
gene1_tpm             (numeric)
gene1_log2tpm         (numeric)
gene2_tpm             (numeric)
gene2_log2tpm         (numeric)
log2_ratio            (numeric)
is_hotspot            (logical)
hotspot_reason        (character)
```

---

### 2. report_generation.R

**Function:** `generate_fusion_report()`

**Major changes:**

1. **Updated function signature:**
   ```r
   generate_fusion_report <- function(
     fusions, 
     expression_hotspots, 
     breakpoint_annotations,
     disease_genes = NULL  # [NEW PARAMETER]
   )
   ```

2. **Fusion Summary:**
   - Updated to use `bp1_gene_symbol`, `bp2_gene_symbol` (not gene1/gene2)
   - Extracts chromosome and position info from `breakpoint_annotations` via join on `variant_id`
   - Disease gene flags are looked up from the `disease_genes` parameter using `bp1_gene_id` and `bp2_gene_id`

   **Output columns:**
   ```
   fusion_id                      (character)
   bp1_gene_symbol               (character)
   bp2_gene_symbol               (character)
   fusion_type                   (character)
   reading_frame                 (character) - "In-frame"/"Out-of-frame"/"Unknown"
   disease_genes                 (character) - comma-separated disease gene symbols (if any)
   disease_gene_count            (integer)
   bp1_seqname                   (character) - chromosome of bp1
   bp1_pos                       (integer)   - position of bp1
   bp2_seqname                   (character) - chromosome of bp2
   bp2_pos                       (integer)   - position of bp2
   is_canonical                  (logical)
   in_frame                      (logical)
   hotspot_tissue_count          (integer)
   hotspot_tissues               (character) - semicolon-separated tissue list
   ```

3. **Breakpoint Summary:**
   - Updated to work with new `breakpoint_annotations` format (gene-level summary)
   - Columns now match the new annotation structure: variant_id, breakpoint, seqname, pos, gene_id, gene_symbol, gene_strand, is_exon, num_exon, is_intron, num_intron, is_cds, is_disease_gene

   **Output columns:**
   ```
   variant_id              (character)
   breakpoint              (integer) - 1 or 2
   seqname                 (character) - chromosome
   pos                     (integer)  - breakpoint position
   gene_id                 (character) - ENSG ID
   gene_symbol             (character) - gene name
   gene_strand             (character) - +/-
   is_exon                 (logical)
   num_exon                (integer)
   is_intron               (logical)
   num_intron              (integer)
   is_cds                  (logical)
   is_disease_gene         (logical)
   region_type             (character) - "Exon"/"Intron"/"CDS"/"Intergenic"
   gene_disruption_summary (character) - human-readable gene disruption description
   ```

4. **Hotspot Summary:** (unchanged)
   - Renamed from `gene1_/gene2_` to clarify they represent bp1/bp2

5. **Report Metadata:**
   - Updated to count canonical/non-canonical using `is_canonical` column
   - In-frame count uses `in_frame` column
   - Disease gene count derived from fusion_summary

---

### 3. fusion_prediction.R

**Function:** `predict_fusions()`

**Additions:**
- Added `variant_id` column to output tibble (required for report generation joins)
- This allows downstream modules to track which variant each fusion came from

---

### 4. fusion_pipeline.qmd

**Updates:**

1. **Phase 4 (Fusion Prediction):**
   - Already compatible with new column names (bp1_gene_symbol, bp2_gene_symbol, is_canonical, in_frame)
   - No changes needed

2. **Phase 5 (Expression Analysis):**
   - Already compatible with identify_fusion_hotspots()
   - No changes needed

3. **Phase 6 (Report Generation):**
   - Updated to pass `disease_genes` parameter:
   ```r
   report <- generate_fusion_report(
     fusions = fusions,
     expression_hotspots = expression_hotspots,
     breakpoint_annotations = breakpoint_annotations,
     disease_genes = disease_genes  # [NEWLY ADDED]
   )
   ```

---

### 5. example_pipeline.R

**Updates:**
- Updated Phase 6 section to pass `disease_genes` to `generate_fusion_report()`
- Same pattern as fusion_pipeline.qmd

---

## Data Flow Verification

```
Phase 1-2: Input → variants_processed, gencode_data, disease_genes
     ↓
Phase 3: annotate_breakpoints()
     → breakpoint_annotations
         (variant_id, breakpoint, seqname, pos, gene_id, gene_symbol, 
          gene_strand, is_exon, num_exon, is_intron, num_intron, is_cds, is_disease_gene)
     ↓
Phase 4: predict_fusions()
     → fusions
         (fusion_id, variant_id, variant_type, bp1_gene_id, bp1_gene_symbol, 
          bp2_gene_id, bp2_gene_symbol, bp1_gene_strand, bp2_gene_strand,
          is_canonical, canonical_fusion_gene_promoter_breakpoint, in_frame, frame_offset)
     ↓
Phase 5: identify_fusion_hotspots()
     → expression_hotspots
         (fusion_id, gene1_id, gene1_symbol, gene2_id, gene2_symbol, tissue,
          gene1_tpm, gene1_log2tpm, gene2_tpm, gene2_log2tpm, log2_ratio, is_hotspot, hotspot_reason)
     ↓
Phase 6: generate_fusion_report()
     → report (list with fusion_summary, breakpoint_summary, hotspot_summary, report_metadata)
```

---

## Testing Checklist

- [x] expression_analysis.R accepts bp1_/bp2_ prefixed columns
- [x] report_generation.R handles new fusion_prediction output format
- [x] report_generation.R joins breakpoint positions via variant_id
- [x] report_generation.R looks up disease genes correctly
- [x] breakpoint_summary works with gene-level annotations (not transcript-level)
- [x] fusion_pipeline.qmd passes disease_genes to report generation
- [x] example_pipeline.R passes disease_genes to report generation
- [x] variant_id is available in fusions for downstream joins

---

## Key Improvements

1. **Clearer naming:** `bp1_`/`bp2_` prefixes explicitly show breakpoint associations
2. **Complete coordinate info:** Both breakpoint chromsomes and positions preserved in reports
3. **Disease gene integration:** Disease gene information properly integrated from lookup tables
4. **Variant tracking:** variant_id now flows through entire pipeline for traceability
5. **Gene-level efficiency:** Breakpoint annotations no longer create transcript explosion
6. **Reading frame preservation:** Original 3-condition in-frame logic maintained throughout

