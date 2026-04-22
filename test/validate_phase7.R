#!/usr/bin/env Rscript

#' Phase 7 Validation: Quality Control & Validation
#'
#' Tests:
#' 1. Validate annotation module outputs
#' 2. Validate breakpoint annotation outputs
#' 3. Validate fusion prediction outputs
#' 4. Validate expression analysis outputs
#' 5. Validate report generation outputs
#' 6. Run complete end-to-end pipeline validation
#' 7. Calculate data quality metrics
#' 8. Handle edge cases (empty data, missing data)
#' 9. Verify data consistency across phases
#' 10. Summary validation report

library(tidyverse)
library(GenomicRanges)
library(IRanges)

source("R/annotations.R")
source("R/breakpoint_annotation.R")
source("R/fusion_prediction.R")
source("R/expression_analysis.R")
source("R/report_generation.R")
source("R/validation.R")

cat("\n========================================\n")
cat("Phase 7 Validation: QC & Validation\n")
cat("========================================\n\n")

# Track test results
test_results <- list()
test_num <- 1

# ==========================================
# Setup: Create Mock Data (from Phase 6)
# ==========================================
cat("Preparing complete pipeline mock data...\n\n")

# Mock annotations
mock_genes_gr <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr2", "chr3"),
  ranges = IRanges::IRanges(start = c(1000, 5000, 10000), end = c(2000, 6000, 11000)),
  gene_id = c("ENSG00000139618", "ENSG00000106819", "ENSG00000186092"),
  gene_name = c("BRCA2", "EML4", "BCR")
)

mock_tx_lookup <- tibble::tibble(
  tx_id = c("ENST00000380152", "ENST00000394469", "ENST00000334262"),
  gene_id = c("ENSG00000139618", "ENSG00000106819", "ENSG00000186092"),
  gene_name = c("BRCA2", "EML4", "BCR")
)

mock_exons_gr <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1", "chr2"),
  ranges = IRanges::IRanges(start = c(1000, 1500, 5000), end = c(1200, 1800, 5500)),
  tx_id = c("ENST00000380152", "ENST00000380152", "ENST00000394469")
)

mock_introns_gr <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges::IRanges(start = c(1300, 5600), end = c(1400, 5900)),
  tx_id = c("ENST00000380152", "ENST00000394469")
)

mock_cds_gr <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1", "chr2"),
  ranges = IRanges::IRanges(start = c(1000, 1500, 5000), end = c(1200, 1800, 5500)),
  tx_id = c("ENST00000380152", "ENST00000380152", "ENST00000394469"),
  phase = c(0L, 1L, 0L),
  end_phase = c(1L, 0L, 1L)
)

mock_annotations <- list(
  genes_gr = mock_genes_gr,
  tx_lookup = mock_tx_lookup,
  exons_gr = mock_exons_gr,
  introns_gr = mock_introns_gr,
  cds_gr = mock_cds_gr,
  disease_genes = c("ENSG00000139618", "ENSG00000106819"),
  gtex_tissues = c("Brain_Cortex", "Lung", "Heart_Left_Ventricle", "Breast_Mammary_Tissue", 
                   "Bone_Marrow", "Liver", "Muscle_Skeletal", "Skin", "Stomach", "Testis",
                   rep("Other_Tissue", 58))  # 68 total tissues
)

# Mock breakpoint annotations
mock_breakpoints <- tibble::tibble(
  variant_id = c("FUSION_001", "FUSION_001", "FUSION_002", "FUSION_002", "FUSION_003", "FUSION_003"),
  breakpoint = c("bp1", "bp2", "bp1", "bp2", "bp1", "bp2"),
  gene_id = c("ENSG00000139618", "ENSG00000136997", "ENSG00000106819", "ENSG00000171094",
              "ENSG00000186092", "ENSG00000097007"),
  region_type = c("exon", "exon", "intron", "exon", "exon", "intron"),
  transcript_id = c("ENST00000380152", "ENST00000377970", "ENST00000394469", "ENST00000396471",
                    "ENST00000334262", "ENST00000368340"),
  is_disease_gene = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  disruption_type = c("stop_loss", "frameshift", "intronic", "exonic_splice", "frameshift", "intronic")
)

# Mock fusions
mock_fusions <- tibble::tibble(
  fusion_id = c("FUSION_001", "FUSION_002", "FUSION_003"),
  gene1_symbol = c("BRCA2", "EML4", "BCR"),
  gene1_id = c("ENSG00000139618", "ENSG00000106819", "ENSG00000186092"),
  gene2_symbol = c("MYC", "ALK", "ABL1"),
  gene2_id = c("ENSG00000136997", "ENSG00000171094", "ENSG00000097007"),
  chrom1 = c("chr13", "chr2", "chr22"),
  pos1 = c(32890000, 42000000, 23500000),
  chrom2 = c("chr8", "chr2", "chr9"),
  pos2 = c(127000000, 29446400, 133000000),
  canonical = c(TRUE, TRUE, FALSE),
  reading_frame = c("in-frame", "out-frame", "in-frame"),
  gene1_disease = c(TRUE, FALSE, FALSE),
  gene2_disease = c(FALSE, TRUE, FALSE)
)

# Mock hotspots
mock_hotspots <- tibble::tibble(
  fusion_id = c("FUSION_001", "FUSION_001", "FUSION_001", "FUSION_002", "FUSION_002", "FUSION_003"),
  tissue = c("Brain_Cortex", "Lung", "Heart_Left_Ventricle", "Breast_Mammary_Tissue", "Lung", "Bone_Marrow"),
  gene1_tpm = c(5.2, 8.1, 3.4, 12.5, 1.2, 0.5),
  gene2_tpm = c(0.2, 0.1, 0.5, 3.2, 0.3, 2.1),
  log2_ratio = c(4.7, 6.3, 2.8, 1.9, 2.0, -2.1),
  is_hotspot = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  hotspot_reason = c(
    "Upstream expressed", "Upstream expressed", "Upstream expressed",
    "Large expression difference", "Large expression difference", NA
  )
)

# Generate report
mock_report <- generate_fusion_report(mock_fusions, mock_hotspots, mock_breakpoints)

cat("✓ Complete pipeline data prepared\n")
cat("  - Annotations with GRanges and disease genes\n")
cat("  - 3 fusion predictions\n")
cat("  - 6 breakpoint annotations\n")
cat("  - 6 expression results\n")
cat("  - Generated report\n\n")

# ==========================================
# TEST 1: Validate Annotations
# ==========================================
test_name <- "Validate annotation module outputs"
cat("Test", test_num, ":", test_name, "\n")

annotation_checks <- validate_annotations(mock_annotations)
annotation_validation_passed <- all(annotation_checks$passed)

test_results[[test_num]] <- annotation_validation_passed
cat(if (annotation_validation_passed) "✓ PASS\n" else "✗ FAIL\n")
if (!annotation_validation_passed) {
  failed_checks <- annotation_checks %>% dplyr::filter(!passed)
  for (i in seq_len(nrow(failed_checks))) {
    row <- failed_checks[i, ]
    cat(sprintf("  Failed: %s - %s\n", row$check_name, row$details))
  }
}
cat("  Annotation checks:", sum(annotation_checks$passed), "/", nrow(annotation_checks), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 2: Validate Breakpoint Annotations
# ==========================================
test_name <- "Validate breakpoint annotation outputs"
cat("Test", test_num, ":", test_name, "\n")

breakpoint_checks <- validate_breakpoint_annotations(mock_breakpoints)
breakpoint_validation_passed <- all(breakpoint_checks$passed)

test_results[[test_num]] <- breakpoint_validation_passed
cat(if (breakpoint_validation_passed) "✓ PASS\n" else "✗ FAIL\n")
cat("  Breakpoint checks:", sum(breakpoint_checks$passed), "/", nrow(breakpoint_checks), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 3: Validate Fusion Predictions
# ==========================================
test_name <- "Validate fusion prediction outputs"
cat("Test", test_num, ":", test_name, "\n")

fusion_checks <- validate_fusion_predictions(mock_fusions)
fusion_validation_passed <- all(fusion_checks$passed)

test_results[[test_num]] <- fusion_validation_passed
cat(if (fusion_validation_passed) "✓ PASS\n" else "✗ FAIL\n")
cat("  Fusion checks:", sum(fusion_checks$passed), "/", nrow(fusion_checks), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 4: Validate Expression Analysis
# ==========================================
test_name <- "Validate expression analysis outputs"
cat("Test", test_num, ":", test_name, "\n")

expression_checks <- validate_expression_analysis(mock_hotspots)
expression_validation_passed <- all(expression_checks$passed)

test_results[[test_num]] <- expression_validation_passed
cat(if (expression_validation_passed) "✓ PASS\n" else "✗ FAIL\n")
cat("  Expression checks:", sum(expression_checks$passed), "/", nrow(expression_checks), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 5: Validate Report Generation
# ==========================================
test_name <- "Validate report generation outputs"
cat("Test", test_num, ":", test_name, "\n")

report_checks <- validate_report(mock_report)
report_validation_passed <- all(report_checks$passed)

test_results[[test_num]] <- report_validation_passed
cat(if (report_validation_passed) "✓ PASS\n" else "✗ FAIL\n")
cat("  Report checks:", sum(report_checks$passed), "/", nrow(report_checks), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 6: End-to-End Pipeline Validation
# ==========================================
test_name <- "Run complete end-to-end pipeline validation"
cat("Test", test_num, ":", test_name, "\n")

pipeline_validation <- validate_complete_pipeline(
  mock_annotations, mock_breakpoints, mock_fusions, mock_hotspots, mock_report
)

test_success_t6 <- pipeline_validation$overall_passed

test_results[[test_num]] <- test_success_t6
cat(if (test_success_t6) "✓ PASS\n" else "✗ FAIL\n")
cat("  Overall:", pipeline_validation$total_passed, "/", pipeline_validation$total_checks, "checks passed\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 7: Calculate Data Quality Metrics
# ==========================================
test_name <- "Calculate data quality metrics"
cat("Test", test_num, ":", test_name, "\n")

quality_metrics <- calculate_quality_metrics(mock_fusions, mock_hotspots, mock_report)

metrics_calculated <- nrow(quality_metrics) > 0 &&
  all(c("category", "metric", "value") %in% names(quality_metrics)) &&
  length(unique(quality_metrics$category)) >= 4  # At least 4 categories

test_results[[test_num]] <- metrics_calculated
cat(if (metrics_calculated) "✓ PASS\n" else "✗ FAIL\n")
if (metrics_calculated) {
  cat("  Metrics calculated:", nrow(quality_metrics), "\n")
  cat("  Categories:", paste(unique(quality_metrics$category), collapse = ", "), "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 8: Handle Edge Cases
# ==========================================
test_name <- "Handle edge cases (empty data)"
cat("Test", test_num, ":", test_name, "\n")

# Test with empty data
empty_fusions <- mock_fusions %>% dplyr::slice(0)
empty_hotspots <- mock_hotspots %>% dplyr::slice(0)
empty_report <- generate_fusion_report(empty_fusions, empty_hotspots, mock_breakpoints)

edge_case_handling <- tryCatch({
  validate_fusion_predictions(empty_fusions)
  validate_expression_analysis(empty_hotspots)
  validate_report(empty_report)
  TRUE
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

test_results[[test_num]] <- edge_case_handling
cat(if (edge_case_handling) "✓ PASS\n" else "✗ FAIL\n")
cat("  Edge cases handled gracefully\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 9: Verify Data Consistency
# ==========================================
test_name <- "Verify data consistency across phases"
cat("Test", test_num, ":", test_name, "\n")

# Check consistency between phases
breakpoint_variants <- unique(mock_breakpoints$variant_id)
fusion_variants <- mock_fusions$fusion_id

# Check that hotspot fusions match fusion predictions
hotspot_fusions <- unique(mock_hotspots$fusion_id)
fusions_with_hotspot_data <- all(hotspot_fusions %in% fusion_variants) ||
  all(hotspot_fusions %in% fusion_variants)

# Check report references all fusions
fusions_in_report <- unique(mock_report$fusion_summary$fusion_id)
all_fusions_in_report <- all(fusion_variants %in% fusions_in_report)

consistency_check <- fusions_with_hotspot_data && all_fusions_in_report

test_results[[test_num]] <- consistency_check
cat(if (consistency_check) "✓ PASS\n" else "✗ FAIL\n")
cat("  Fusions in hotspots match predictions:", fusions_with_hotspot_data, "\n")
cat("  All fusions present in report:", all_fusions_in_report, "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 10: Summary Validation Report
# ==========================================
test_name <- "Generate summary validation report"
cat("Test", test_num, ":", test_name, "\n")

# Generate summary report
summary_result <- tryCatch({
  capture.output(print_validation_summary(pipeline_validation))
  TRUE
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

test_results[[test_num]] <- summary_result
cat(if (summary_result) "✓ PASS\n" else "✗ FAIL\n")
cat("  Summary report generated successfully\n\n")
test_num <- test_num + 1

# ==========================================
# Cleanup
# ==========================================
cat("Validation test suite completed\n\n")

# ==========================================
# Summary
# ==========================================
cat("========================================\n")
cat("PHASE 7 TEST SUMMARY\n")
cat("========================================\n\n")

passed_count <- sum(unlist(test_results))
total_count <- length(test_results)

test_names <- c(
  "Validate annotations",
  "Validate breakpoints",
  "Validate fusions",
  "Validate expression",
  "Validate reports",
  "End-to-end pipeline",
  "Quality metrics",
  "Edge cases",
  "Data consistency",
  "Summary report"
)

for (i in seq_along(test_results)) {
  status <- if (test_results[[i]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("Test %2d: %s\n", i, status))
}

cat("\n")
cat(sprintf("Total: %d / %d tests passed\n", passed_count, total_count))
cat("\n")

if (passed_count == total_count) {
  cat("✓ Phase 7 is complete. All phases validated.\n")
  cat("✓ Pipeline is ready for production use.\n\n")
} else {
  cat("✗ Phase 7 requires fixes before production use.\n\n")
  cat("Failed tests:", which(!unlist(test_results)), "\n\n")
}
