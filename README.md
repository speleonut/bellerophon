# Bellerophon: Fusion Gene Detection Pipeline

The slayer of chimeras! A structural variant analysis pipeline that predicts gene fusions, classifies reading frames, and integrates GTEx expression data.

## Project Version and use of AI declaration
This project is currently in alpha stage. You can expect a lot of things to be broken and be wary that some outputs may be inaccurate.
The development phase of this project was carried out with Claude Haiku 4.5. The meatbag of this operation is still in the process of slowly going through it throughly to review the scripts, test that they function with the expected outputs (including edge cases) and refining the outputs.

usage: 
1. You'll need to download the GTEx data and let the scripts download the GENCODE data and place them in data/annotations. Please see the 'data/annotations/annotation_files.README.md' for information
2. Open up the 'fusion_pipeline.qmd' file in RStudio and save a copy of it that you can edit.
3. Prepare your input file (MAVIS format is your best bet, VCF input is likely to fail as the script has only had to deal with a very simple test file so far.)
4. Change parameters as required
   - `input_file`: Path to your variant file
   - `output_dir`: Directory for results
5. Pray to your preferred deity or perform any other luck related rituals you prefer
6. Run all chunks if you're feeling very lucky or step through them one at a time if you're more pessimistic (the latter is recommended at this stage).
7. See what happens and if you can make anything sensible out of the results.

**Everything below this line is AI generated and definitely not 100% accurate**
## Project Structure

```
bellerophon/
├── fusion_pipeline.qmd         # Main Quarto document (NEEDS UPDATE - see notes below)
├── R/                          # Helper function modules
│   ├── annotations.R          # Load GENCODE v49, disease genes, GTEx metadata (Phase 2)
│   ├── breakpoint_annotation.R # Identify disrupted genes/transcripts (Phase 3)
│   ├── fusion_prediction.R    # Predict canonical/non-canonical fusions (Phase 4)
│   ├── expression_analysis.R  # GTEx integration and heatmaps (Phase 5)
│   ├── report_generation.R    # Summary tables and reports (Phase 6)
│   └── validation.R           # Quality control and validation (Phase 7)
├── test/
│   ├── validate_phase1.R          # Test suite for Phase 1 (placeholder - integrated into Phase 2)
│   ├── validate_phase2.R          # Test suite for Phase 2
│   ├── validate_phase3.R          # Test suite for Phase 3
│   ├── validate_phase4.R          # Test suite for Phase 4
│   ├── validate_phase5.R          # Test suite for Phase 5
│   ├── validate_phase6.R          # Test suite for Phase 6
│   └── validate_phase7.R          # Test suite for Phase 7
├── data/
│   ├── annotations
│   │   ├── Nijmegen.DG.ENSG.list.txt  # Disease gene list
│   │   ├── gencode.v49.basic.annotation.gff3.gz # Download from GENCODE
│   │   └── GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz # Download from GTEx
│   └── example_input/              # Sample test files
│       ├── sample_variants.vcf
│       ├── sample_variants.mavis.tsv
│       └── sample_variants.hgvs.txt
├── data/
│   ├── *various AI generated reports and notes from the development phase*
├── output/                     # Generated reports and results
└── README.md                   # This file
```

## Implementation Status

### Phase 1: ✅ Complete (Core Utilities)
- [x] Quarto document skeleton with parameters
- [x] HGVS nomenclature parser (genomic notation: g.)
- [x] Chromosome normalization and validation
- [x] Uncertain breakpoint position handling
- [x] Variant type classification
- [x] Core data structure definitions

### Phase 2: ✅ Complete (Data Loading)
- [x] Load GENCODE annotations via AnnotationHub
- [x] Parse GTF into transcript/exon structure
- [x] Load disease gene list
- [x] Initialize GTEx metadata
- [x] Comprehensive test suite (validate_phase2.R)
- [x] Error handling and fallback mechanisms

### Phase 3: ✅ Complete (Breakpoint Annotation)
- [x] Identify genes at breakpoint 1 and 2
- [x] Classify whether breakpoint falls in exon/intron/UTR
- [x] Flag disease genes
- [x] Handle multi-gene breakpoints
- [x] Intergenic breakpoint detection
- [x] Comprehensive test suite (validate_phase3.R)
- [x] Support for flexible input data formats

### Phase 4: ✅ Complete (Fusion Prediction)
- [x] Detect canonical fusions (two distinct genes)
- [x] Detect non-canonical fusions (splicing into same/downstream exons)
- [x] Calculate reading frame offset using CDS phase metadata
- [x] Generate fusion records with canonical/non-canonical classification
- [x] Phase-based reading frame calculation (in-frame/out-frame)
- [x] Comprehensive test suite (validate_phase4.R)

### Phase 5: ✅ Complete (Expression Analysis)
- [x] Query GTEx v11 expression data (74,628 genes × 68 tissues)
- [x] Calculate log2-transformed TPM values (log2(TPM+0.1))
- [x] Compute tissue-specific expression ratios
- [x] Identify tissue-level fusion hotspots (2 conditions)
  - Condition 1: Upstream gene expressed (TPM >= 1) AND downstream silent (TPM < 1)
  - Condition 2: Large expression difference (log2 ratio > 1)
- [x] Generate expression heatmaps with pheatmap
- [x] Handle version-numbered gene IDs (ENSG.15 → ENSG)
- [x] Comprehensive test suite (validate_phase5.R)

### Phase 6: ✅ Complete (Report Generation)
- [x] Generate fusion summary (1 row per fusion with type/reading frame/hotspot count)
- [x] Generate breakpoint summary (detailed per-breakpoint information)
- [x] Generate hotspot summary (tissue-level expression analysis)
- [x] Compute report metadata (statistics and counts)
- [x] Export TSV tables for downstream analysis
- [x] Generate HTML reports with styled interactive tables
- [x] Print console-friendly summary reports
- [x] Comprehensive test suite (validate_phase6.R)

### Phase 7: ✅ Complete (Quality Control & Validation)
- [x] Validate annotation outputs (GRanges, disease genes, GTEx metadata)
- [x] Validate breakpoint annotations (structure, pairing, consistency)
- [x] Validate fusion predictions (canonical status, reading frames, gene names)
- [x] Validate expression analysis (TPM values, hotspot flags, reasons)
- [x] Validate report outputs (sections, calculations, consistency)
- [x] Run end-to-end pipeline validation across all phases
- [x] Calculate data quality metrics (discovery, reading frames, hotspots)
- [x] Handle edge cases (empty data, missing values)
- [x] Verify data consistency across phases
- [x] Generate validation summary reports
- [x] Comprehensive test suite (validate_phase7.R)

## Requirements

R packages (install prerequisites):
```r
install.packages(c("tidyverse", "stringr", "stringi", "R.utils"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "VariantAnnotation",
  "GenomicRanges",
  "GenomicFeatures",
  "txdbmaker",
  "biomaRt",
  "AnnotationHub",
  "Biostrings",
  "BSgenome",
  "SummarizedExperiment",
  "StructuralVariantAnnotation"
))

install.packages("gtexr")
install.packages("pheatmap")
install.packages("glue")
install.packages("readr")
```

## Running the Pipeline

### Interactive Mode (RStudio/VS Code)
1. Open `fusion_pipeline.qmd` in RStudio or VS Code
2. Click "Render" to execute the full pipeline
3. Adjust parameters in the YAML header as needed:
   - `input_file`: Path to your variant file
   - `input_format`: Auto-detect or specify "vcf", "mavis", or "hgvs"
   - `genome_build`: GRCh38 (only option currently supported)
   - `output_dir`: Directory for results

### Command Line Mode
```bash
quarto render fusion_pipeline.qmd \
  -P input_file:data/example_input/sample_variants.vcf \
  -P output_dir:output/my_analysis
```

## Input Formats

### VCF (v4.2+)
Standard VCF format with structural variant records.

### MAVIS Format
TSV/CSV with required columns:
- break1_chromosome, break1_position_start, break1_position_end
- break2_chromosome, break2_position_start, break2_position_end
- break1_orientation, break2_orientation
- event_type

### HGVS Nomenclature
Plain text file with one variant per line in HGVS genomic notation.
Currently supports simple variants (deletions, duplications, inversions, insertions).

Examples:
```
NC_000001.11:g.100_200del
NC_000001.11:g.300_400dup
NC_000001.11:g.500_600inv
```

## Output

### HTML Report
- Interactive summary tables
- Variant type distribution
- Fusion discovery statistics
- Expression heatmaps (in-frame fusions)

### RDS Files
- `variants_processed.rds` - Processed variants with IDs and classifications
- Additional intermediate results at each phase

## Known Limitations & Notes

### Important: fusion_pipeline.qmd Requires Update
The `fusion_pipeline.qmd` Quarto document was the initial coordination document but is **not currently integrated with Phases 1-7 modules**. To use the pipeline:

**Option 1 (Recommended)**: Use validation test suites directly
- Run `validate_phase2.R` through `validate_phase7.R` sequentially to test individual phases
- Load modules manually with `source("R/annotations.R")`, etc.

**Option 2**: Create updated Quarto document
- Integrate Phase modules into a single Quarto document for production use
- This would combine all 7 phases into unified analysis workflow

### Other Known Limitations

- GRCh38/hg38 genome build only
- Requires external annotation and expression data files (see data/annotations/ requirements)

## Development Notes

### HGVS Parsing
The HGVS parser in `R/parse_hgvs.R` uses regex patterns to extract:
- Chromosome (supports NCBI contig IDs: NC_000001.11, etc.)
- Variant type (deletion, duplication, inversion, insertion)
- Breakpoint positions
- Reference and alternate alleles (for substitutions)

NCBI to chromosome mapping covers GRCh38 and alternative versions.

### Uncertain Breakpoint Handling
When breakpoint is represented as a range, the selection logic is:
- **Forward strand (+)**: Use maximum end position (rightmost)
- **Reverse strand (-)**: Use minimum start position (leftmost)

This ensures the widest possible impact margin for fusion prediction.

### Logging
All operations are logged with timestamps. Messages are stored in `log_messages` list and printed to console/output document.

## Testing

Test files are provided in `data/example_input/`:
- `sample_variants.vcf` - 4 test variants in VCF format
- `sample_variants.mavis.tsv` - 5 test variants in MAVIS format
- `sample_variants.hgvs.txt` - 5 test variants in HGVS format

Run Quarto with default parameters to test:
```bash
quarto render fusion_pipeline.qmd
```

## Planned Features

- [ ] Support for coding (c.) and protein (p.) HGVS notation
- [ ] Complex variant handling
- [ ] High-throughput parallel processing (furrr)
- [ ] Shiny app for interactive input/output
- [ ] VCF export of predictions
- [ ] Integration with InterProScan for domain analysis
- [ ] ML-based splice site prediction
- [ ] Liftover support for other genome builds

## References

- HGVS Nomenclature: https://hgvs-nomenclature.org/
- GENCODE: https://www.gencodegenes.org/human/
- GTEx: https://www.gtexportal.org/
- MAVIS Documentation: https://mavis.readthedocs.io/

## License

(Specify license here)

## Contact

For questions or issues, contact project maintainers.
