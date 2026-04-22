#!/usr/bin/env Rscript

#' Phase 6 Validation: Report Generation
#'
#' Tests:
#' 1. Generate basic fusion report from mock data
#' 2. Validate fusion summary table structure
#' 3. Validate breakpoint summary table
#' 4. Validate hotspot summary aggregation
#' 5. Verify report metadata calculations
#' 6. Export report tables to TSV
#' 7. Generate HTML report without errors
#' 8. Print report summary to console
#' 9. Handle missing hotspot data gracefully
#' 10. Integration test with realistic data

library(tidyverse)
library(GenomicRanges)
library(IRanges)

source("R/report_generation.R")

cat("\n========================================\n")
cat("Phase 6 Validation: Report Generation\n")
cat("========================================\n\n")

# Track test results
test_results <- list()
test_num <- 1

# ==========================================
# Setup: Create Mock Data
# ==========================================
cat("Preparing test data...\n\n")

# Mock fusion predictions (Phase 4 output)
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
  gene2_disease = c(FALSE, TRUE, FALSE),
  breakpoint_annotations_1 = c("exon", "intron", "exon"),
  breakpoint_annotations_2 = c("exon", "exon", "intron")
)

# Mock expression hotspots (Phase 5 output)
mock_hotspots <- tibble::tibble(
  fusion_id = c(
    "FUSION_001", "FUSION_001", "FUSION_001",
    "FUSION_002", "FUSION_002",
    "FUSION_003"
  ),
  tissue = c(
    "Brain_Cortex", "Lung", "Heart_Left_Ventricle",
    "Breast_Mammary_Tissue", "Lung",
    "Bone_Marrow"
  ),
  gene1_tpm = c(5.2, 8.1, 3.4, 12.5, 1.2, 0.5),
  gene2_tpm = c(0.2, 0.1, 0.5, 3.2, 0.3, 2.1),
  log2_ratio = c(4.7, 6.3, 2.8, 1.9, 2.0, -2.1),
  is_hotspot = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  hotspot_reason = c(
    "Upstream expressed (gene1_tpm >= 1) and downstream silent (gene2_tpm < 1)",
    "Upstream expressed (gene1_tpm >= 1) and downstream silent (gene2_tpm < 1)",
    "Upstream expressed (gene1_tpm >= 1) and downstream silent (gene2_tpm < 1)",
    "Large expression difference (log2 ratio > 1)",
    "Large expression difference (log2 ratio > 1)",
    NA
  )
)

# Mock breakpoint annotations (Phase 3 output)
mock_breakpoints <- tibble::tibble(
  variant_id = c(
    "FUSION_001", "FUSION_001",
    "FUSION_002", "FUSION_002",
    "FUSION_003", "FUSION_003"
  ),
  breakpoint = c(
    "bp1", "bp2",
    "bp1", "bp2",
    "bp1", "bp2"
  ),
  gene_id = c(
    "ENSG00000139618", "ENSG00000136997",
    "ENSG00000106819", "ENSG00000171094",
    "ENSG00000186092", "ENSG00000097007"
  ),
  region_type = c(
    "exon", "exon",
    "intron", "exon",
    "exon", "intron"
  ),
  transcript_id = c(
    "ENST00000380152", "ENST00000377970",
    "ENST00000394469", "ENST00000396471",
    "ENST00000334262", "ENST00000368340"
  ),
  is_disease_gene = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  disruption_type = c(
    "stop_loss", "frameshift",
    "intronic", "exonic_splice",
    "frameshift", "intronic"
  )
)

cat("✓ Test data prepared\n")
cat("  - 3 fusions (2 canonical, 1 non-canonical)\n")
cat("  - 6 expression tissue results (5 hotspots)\n")
cat("  - 6 breakpoint annotations\n\n")

# ==========================================
# TEST 1: Generate Basic Fusion Report
# ==========================================
test_name <- "Generate basic fusion report from mock data"
cat("Test", test_num, ":", test_name, "\n")

report <- generate_fusion_report(mock_fusions, mock_hotspots, mock_breakpoints)

test_success_t1 <- !is.null(report) &&
  is.list(report) &&
  all(c("fusion_summary", "breakpoint_summary", "hotspot_summary", "report_metadata") %in% names(report)) &&
  nrow(report$fusion_summary) == 3

test_results[[test_num]] <- test_success_t1
cat(if (test_success_t1) "✓ PASS\n" else "✗ FAIL\n")
cat("  Report contains", length(names(report)), "sections\n")
cat("  Fusions in summary:", nrow(report$fusion_summary), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 2: Validate Fusion Summary Structure
# ==========================================
test_name <- "Validate fusion summary table structure"
cat("Test", test_num, ":", test_name, "\n")

required_cols_fusion <- c(
  "fusion_id", "gene1_symbol", "gene2_symbol", "fusion_type",
  "reading_frame", "disease_genes", "disease_gene_count",
  "chrom1", "pos1", "chrom2", "pos2",
  "hotspot_tissue_count", "hotspot_tissues"
)

test_success_t2 <- all(required_cols_fusion %in% names(report$fusion_summary)) &&
  nrow(report$fusion_summary) == 3 &&
  all(report$fusion_summary$disease_gene_count %in% c(0, 1, 2)) &&
  !any(is.na(report$fusion_summary$fusion_id))

test_results[[test_num]] <- test_success_t2
cat(if (test_success_t2) "✓ PASS\n" else "✗ FAIL\n")
cat("  Columns present:", sum(required_cols_fusion %in% names(report$fusion_summary)), "/", length(required_cols_fusion), "\n")
if (nrow(report$fusion_summary) > 0) {
  cat("  Sample (first fusion):\n")
  print(report$fusion_summary %>% dplyr::slice(1))
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 3: Validate Breakpoint Summary
# ==========================================
test_name <- "Validate breakpoint summary table"
cat("Test", test_num, ":", test_name, "\n")

required_cols_bp <- c(
  "variant_id", "breakpoint", "gene_id", "region_type",
  "transcript_id", "is_disease_gene", "disruption_type", "region_display"
)

test_success_t3 <- all(required_cols_bp %in% names(report$breakpoint_summary)) &&
  nrow(report$breakpoint_summary) == 6 &&
  all(report$breakpoint_summary$breakpoint %in% c("bp1", "bp2")) &&
  all(is.logical(report$breakpoint_summary$is_disease_gene))

test_results[[test_num]] <- test_success_t3
cat(if (test_success_t3) "✓ PASS\n" else "✗ FAIL\n")
cat("  Breakpoint rows:", nrow(report$breakpoint_summary), "\n")
cat("  Columns present:", sum(required_cols_bp %in% names(report$breakpoint_summary)), "/", length(required_cols_bp), "\n")
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 4: Validate Hotspot Summary
# ==========================================
test_name <- "Validate hotspot summary aggregation"
cat("Test", test_num, ":", test_name, "\n")

required_cols_hs <- c(
  "fusion_id", "tissue", "gene1_tpm", "gene2_tpm",
  "log2_ratio", "hotspot_status", "hotspot_reason"
)

test_success_t4 <- all(required_cols_hs %in% names(report$hotspot_summary)) &&
  nrow(report$hotspot_summary) == 6 &&
  sum(report$hotspot_summary$hotspot_status == "Hotspot") == 5 &&
  sum(report$hotspot_summary$hotspot_status == "Not hotspot") == 1

test_results[[test_num]] <- test_success_t4
cat(if (test_success_t4) "✓ PASS\n" else "✗ FAIL\n")
cat("  Total tissue results:", nrow(report$hotspot_summary), "\n")
cat("  Hotspots found:", sum(report$hotspot_summary$hotspot_status == "Hotspot"), "\n")
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 5: Verify Report Metadata Calculations
# ==========================================
test_name <- "Verify report metadata calculations"
cat("Test", test_num, ":", test_name, "\n")

metadata <- report$report_metadata[1, ]

test_success_t5 <- metadata$total_fusions[[1]] == 3 &&
  metadata$canonical_fusions[[1]] == 2 &&
  metadata$non_canonical_fusions[[1]] == 1 &&
  metadata$no_fusion_cases[[1]] == 0 &&
  metadata$predicted_in_frame[[1]] == 2 &&
  metadata$total_disease_genes[[1]] == 2 &&
  metadata$total_tissues_analyzed[[1]] == 5 &&
  metadata$total_hotspots_identified[[1]] == 5 &&
  !is.na(metadata$generated_at[[1]])

test_results[[test_num]] <- test_success_t5
cat(if (test_success_t5) "✓ PASS\n" else "✗ FAIL\n")
cat("  Total fusions:", metadata$total_fusions[[1]], "(expected 3)\n")
cat("  Canonical:", metadata$canonical_fusions[[1]], "(expected 2)\n")
cat("  In-frame:", metadata$predicted_in_frame[[1]], "(expected 2)\n")
cat("  Hotspots:", metadata$total_hotspots_identified[[1]], "(expected 5)\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 6: Export Report Tables to TSV
# ==========================================
test_name <- "Export report tables to TSV files"
cat("Test", test_num, ":", test_name, "\n")

output_dir <- "test_output/"
files_written <- tryCatch({
  export_report_tables(report, output_dir = output_dir)
  TRUE
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

# Check if files were created
files_exist <- all(
  file.exists(file.path(output_dir, dir(output_dir)))
)

test_success_t6 <- files_written && files_exist && length(dir(output_dir)) >= 4

test_results[[test_num]] <- test_success_t6
cat(if (test_success_t6) "✓ PASS\n" else "✗ FAIL\n")
if (files_exist) {
  cat("  Files created:", length(dir(output_dir)), "\n")
  cat("  Files:", paste(dir(output_dir), collapse = ", "), "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 7: Generate HTML Report
# ==========================================
test_name <- "Generate HTML report without errors"
cat("Test", test_num, ":", test_name, "\n")

html_file <- "test_output/test_report.html"

html_result <- tryCatch({
  generate_html_report(report, output_file = html_file, title = "Test Fusion Report")
  file.exists(html_file)
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

test_success_t7 <- html_result

test_results[[test_num]] <- test_success_t7
cat(if (test_success_t7) "✓ PASS\n" else "✗ FAIL\n")
if (html_result) {
  file_size <- file.info(html_file)$size
  cat("  HTML report created:", html_file, "\n")
  cat("  File size:", file_size, "bytes\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 8: Print Report Summary
# ==========================================
test_name <- "Print report summary to console"
cat("Test", test_num, ":", test_name, "\n")

print_result <- tryCatch({
  capture.output(print_report_summary(report))
  TRUE
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

test_success_t8 <- print_result

test_results[[test_num]] <- test_success_t8
cat(if (test_success_t8) "✓ PASS\n" else "✗ FAIL\n")
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 9: Handle Empty Hotspot Data
# ==========================================
test_name <- "Handle missing hotspot data gracefully"
cat("Test", test_num, ":", test_name, "\n")

# Create mock data without hotspots
mock_fusions_no_exp <- mock_fusions
mock_hotspots_empty <- tibble::tibble(
  fusion_id = character(),
  tissue = character(),
  gene1_tpm = numeric(),
  gene2_tpm = numeric(),
  log2_ratio = numeric(),
  is_hotspot = logical(),
  hotspot_reason = character()
)

report_no_exp <- tryCatch({
  generate_fusion_report(mock_fusions_no_exp, mock_hotspots_empty, mock_breakpoints)
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  NULL
})

test_success_t9 <- !is.null(report_no_exp) &&
  nrow(report_no_exp$fusion_summary) == 3 &&
  !any(is.na(report_no_exp$fusion_summary$hotspot_tissue_count)) &&
  all(report_no_exp$fusion_summary$hotspot_tissue_count == 0)

test_results[[test_num]] <- test_success_t9
cat(if (test_success_t9) "✓ PASS\n" else "✗ FAIL\n")
cat("  Report generated with empty hotspots\n")
cat("  All hotspot counts are 0:", all(report_no_exp$fusion_summary$hotspot_tissue_count == 0), "\n\n")
test_num <- test_num + 1

# ==========================================
# TEST 10: Comprehensive Report Validation
# ==========================================
test_name <- "Comprehensive report with all sections"
cat("Test", test_num, ":", test_name, "\n")

# Comprehensive validation
is_valid_report <- 
  # Structure check
  all(c("fusion_summary", "breakpoint_summary", "hotspot_summary", "report_metadata") %in% names(report)) &&
  
  # Content checks
  nrow(report$fusion_summary) > 0 &&
  nrow(report$hotspot_summary) > 0 &&
  nrow(report$report_metadata) == 1 &&
  
  # Data integrity
  all(!is.na(report$fusion_summary$fusion_id)) &&
  all(!is.na(report$breakpoint_summary$variant_id)) &&
  all(!is.na(report$hotspot_summary$fusion_id)) &&
  
  # Logical consistency
  n_distinct(report$fusion_summary$fusion_id) == nrow(report$fusion_summary)

test_success_t10 <- is_valid_report

test_results[[test_num]] <- test_success_t10
cat(if (test_success_t10) "✓ PASS\n" else "✗ FAIL\n")
cat("  Report structure valid:", all(c("fusion_summary", "breakpoint_summary", "hotspot_summary", "report_metadata") %in% names(report)), "\n")
cat("  All required data present:", nrow(report$fusion_summary) > 0 && nrow(report$hotspot_summary) > 0, "\n\n")
test_num <- test_num + 1

# ==========================================
# Cleanup
# ==========================================
cat("Cleaning up test files...\n")
unlink(output_dir, recursive = TRUE)
cat("✓ Test output directory removed\n\n")

# ==========================================
# Summary
# ==========================================
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

passed_count <- sum(unlist(test_results))
total_count <- length(test_results)

test_names <- c(
  "Generate basic report",
  "Validate fusion summary",
  "Validate breakpoint summary",
  "Validate hotspot summary",
  "Verify metadata calculations",
  "Export TSV tables",
  "Generate HTML report",
  "Print report summary",
  "Handle empty hotspots",
  "Comprehensive validation"
)

for (i in seq_along(test_results)) {
  status <- if (test_results[[i]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("Test %2d: %s\n", i, status))
}

cat("\n")
cat(sprintf("Total: %d / %d tests passed\n", passed_count, total_count))
cat("\n")

if (passed_count == total_count) {
  cat("✓ Phase 6 is complete and ready for Phase 7.\n\n")
} else {
  cat("✗ Phase 6 requires fixes before proceeding.\n\n")
  cat("Failed tests:", which(!unlist(test_results)), "\n\n")
}
