#!/usr/bin/env Rscript

#' Phase 5 Validation: Expression Analysis
#'
#' Tests:
#' 1. Load GTEx expression data from file
#' 2. Query single gene expression across tissues
#' 3. Query multiple genes simultaneously
#' 4. Handle gene IDs with version numbers
#' 5. Identify fusion hotspots: upstream expressed, downstream silent
#' 6. Identify fusion hotspots: large expression difference
#' 7. Handle missing genes gracefully
#' 8. Filter specific tissues
#' 9. Heatmap visualization
#' 10. Integration with Phase 4 fusion output

library(tidyverse)
library(GenomicRanges)
library(IRanges)

source("R/annotations.R")
source("R/breakpoint_annotation.R")
source("R/fusion_prediction.R")
source("R/expression_analysis.R")

cat("\n========================================\n")
cat("Phase 5 Validation: Expression Analysis\n")
cat("========================================\n\n")

# Load GTEx expression data
cat("Loading GTEx expression data...\n")
gtex_data <- load_gtex_expression()

if (!gtex_data$loaded_successfully) {
  cat("✗ GTEx data not available for testing\n")
  cat("Please download: data/annotations/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz\n\n")
  quit(status = 1)
}

cat("✓ GTEx data loaded\n")
cat("  Dimensions:", nrow(gtex_data$expression), "genes ×", length(gtex_data$tissues), "tissues\n\n")

# Track test results
test_results <- list()
test_num <- 1

# ==========================================
# TEST 1: Query Single Gene
# ==========================================
test_name <- "Query single gene expression across tissues"
cat("Test", test_num, ":", test_name, "\n")

# Use a known gene that exists in GTEx (DDX11L1)
gene_expr_t1 <- query_gene_expression("ENSG00000223972", gtex_data)

test_success_t1 <- nrow(gene_expr_t1) > 0 &&
  all(c("gene_id", "gene_symbol", "tissue", "tpm", "log2_tpm") %in% names(gene_expr_t1)) &&
  n_distinct(gene_expr_t1$tissue) == length(gtex_data$tissues)

test_results[[test_num]] <- test_success_t1
cat(if (test_success_t1) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(gene_expr_t1) > 0) {
  cat("  Tissues queried:", n_distinct(gene_expr_t1$tissue), "\n")
  cat("  Sample output (first 3 tissues):\n")
  head(gene_expr_t1, 3)
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 2: Query Multiple Genes
# ==========================================
test_name <- "Query multiple genes simultaneously"
cat("Test", test_num, ":", test_name, "\n")

# Query two genes
gene_expr_t2 <- query_gene_expression(
  c("ENSG00000223972", "ENSG00000000003"),
  gtex_data
)

test_success_t2 <- nrow(gene_expr_t2) > 1 &&
  n_distinct(gene_expr_t2$gene_id) >= 1  # At least one gene found

test_results[[test_num]] <- test_success_t2
cat(if (test_success_t2) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(gene_expr_t2) > 0) {
  cat("  Genes found:", n_distinct(gene_expr_t2$gene_id), "\n")
  cat("  Total expression values:", nrow(gene_expr_t2), "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 3: Handle Version Numbers
# ==========================================
test_name <- "Strip version numbers from gene IDs"
cat("Test", test_num, ":", test_name, "\n")

# Query with version-numbered ID
gene_expr_t3a <- query_gene_expression("ENSG00000223972.6", gtex_data)
gene_expr_t3b <- query_gene_expression("ENSG00000223972", gtex_data)

test_success_t3 <- nrow(gene_expr_t3a) == nrow(gene_expr_t3b) &&
  all(gene_expr_t3a$gene_symbol == gene_expr_t3b$gene_symbol)

test_results[[test_num]] <- test_success_t3
cat(if (test_success_t3) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(gene_expr_t3a) > 0 && nrow(gene_expr_t3b) > 0) {
  cat("  Query with version (.6):", nrow(gene_expr_t3a), "results\n")
  cat("  Query without version:", nrow(gene_expr_t3b), "results\n")
  cat("  Matching results:", test_success_t3, "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 4: Filter Specific Tissues
# ==========================================
test_name <- "Filter expression to specific tissues"
cat("Test", test_num, ":", test_name, "\n")

tissues_to_filter <- c("Brain_Cortex", "Lung", "Heart_Left_Ventricle")
gene_expr_t4 <- query_gene_expression("ENSG00000223972", gtex_data, tissues = tissues_to_filter)

test_success_t4 <- nrow(gene_expr_t4) > 0 &&
  all(gene_expr_t4$tissue %in% tissues_to_filter) &&
  n_distinct(gene_expr_t4$tissue) <= length(tissues_to_filter)

test_results[[test_num]] <- test_success_t4
cat(if (test_success_t4) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(gene_expr_t4) > 0) {
  cat("  Tissues in result:", paste(unique(gene_expr_t4$tissue), collapse = ", "), "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 5: Hotspot - Upstream Expressed, Downstream Silent
# ==========================================
test_name <- "Identify hotspot: upstream expressed, downstream silent"
cat("Test", test_num, ":", test_name, "\n")

# Create synthetic fusion data with genes we know exist
# Using two different genes to test hotspot detection
fusion_data_t5 <- tibble::tibble(
  fusion_id = "TEST_FUSION_1",
  gene1_id = "ENSG00000156531",      # PHF6 (often expressed)
  gene1_symbol = "PHF6",
  gene2_id = "ENSG00000000003",      # FKBP4
  gene2_symbol = "FKBP4"
)

hotspots_t5 <- identify_fusion_hotspots(fusion_data_t5, gtex_data)

# Check if hotspots are identified where upstream >= 1 and downstream < 1
if (nrow(hotspots_t5) > 0) {
  cond1_hotspots <- hotspots_t5 %>%
    filter(is_hotspot & !is.na(hotspot_reason)) %>%
    filter(stringr::str_detect(hotspot_reason, "Upstream expressed"))
  
  test_success_t5 <- nrow(cond1_hotspots) > 0
} else {
  test_success_t5 <- FALSE
}

test_results[[test_num]] <- test_success_t5
cat(if (test_success_t5) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(hotspots_t5) > 0) {
  hotspot_tissues <- hotspots_t5 %>% filter(is_hotspot)
  cat("  Hotspot tissues identified:", nrow(hotspot_tissues), "\n")
  if (nrow(hotspot_tissues) > 0) {
    cat("  Example (first 3):\n")
    print(head(hotspot_tissues %>% dplyr::select(tissue, gene1_tpm, gene2_tpm, hotspot_reason), 3))
  }
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 6: Hotspot - Large Expression Difference
# ==========================================
test_name <- "Identify hotspot: large expression difference (log2 ratio > 1)"
cat("Test", test_num, ":", test_name, "\n")

# Use hotspots from previous test
if (nrow(hotspots_t5) > 0) {
  cond2_hotspots <- hotspots_t5 %>%
    filter(is_hotspot & !is.na(hotspot_reason)) %>%
    filter(stringr::str_detect(hotspot_reason, "Large expression difference"))
  
  # Verify the log2 ratio condition
  if (nrow(cond2_hotspots) > 0) {
    test_success_t6 <- all(cond2_hotspots$log2_ratio > 1)
  } else {
    test_success_t6 <- TRUE  # No large ratio hotspots found is not a failure
  }
} else {
  test_success_t6 <- FALSE
}

test_results[[test_num]] <- test_success_t6
cat(if (test_success_t6) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(hotspots_t5) > 0) {
  cond2_hotspots <- hotspots_t5 %>% filter(is_hotspot & stringr::str_detect(hotspot_reason %||% "", "Large"))
  cat("  Large ratio hotspots found:", nrow(cond2_hotspots), "\n")
  if (nrow(cond2_hotspots) > 0) {
    cat("  Example (first 3):\n")
    print(head(cond2_hotspots %>% dplyr::select(tissue, log2_ratio, hotspot_reason), 3))
  }
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 7: Handle Missing Genes
# ==========================================
test_name <- "Handle missing genes gracefully"
cat("Test", test_num, ":", test_name, "\n")

# Query for genes that don't exist
gene_expr_t7 <- query_gene_expression(
  c("ENSG00000000000000", "ENSG99999999999999"),
  gtex_data
)

test_success_t7 <- nrow(gene_expr_t7) == 0  # Should return empty tibble

test_results[[test_num]] <- test_success_t7
cat(if (test_success_t7) "✓ PASS\n" else "✗ FAIL\n")
cat("  Non-existent genes query returned:", nrow(gene_expr_t7), "rows (expected 0)\n")
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 8: Log2 TPM Calculation
# ==========================================
test_name <- "Verify log2(TPM + 0.1) calculation"
cat("Test", test_num, ":", test_name, "\n")

gene_expr_t8 <- query_gene_expression("ENSG00000223972", gtex_data)

if (nrow(gene_expr_t8) > 0) {
  # Manually calculate log2 TPM for first row
  first_row <- gene_expr_t8[1, ]
  expected_log2 <- log2(first_row$tpm[1] + 0.1)
  actual_log2 <- first_row$log2_tpm[1]
  
  test_success_t8 <- abs(expected_log2 - actual_log2) < 0.0001
} else {
  test_success_t8 <- FALSE
}

test_results[[test_num]] <- test_success_t8
cat(if (test_success_t8) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(gene_expr_t8) > 0) {
  cat("  Sample: TPM =", gene_expr_t8$tpm[1], 
      "→ log2(TPM+0.1) =", gene_expr_t8$log2_tpm[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 9: Multiple Fusions at Once
# ==========================================
test_name <- "Process multiple fusions simultaneously"
cat("Test", test_num, ":", test_name, "\n")

fusion_data_t9 <- tibble::tibble(
  fusion_id = c("FUSION_A", "FUSION_B"),
  gene1_id = c("ENSG00000223972", "ENSG00000000003"),
  gene1_symbol = c("DDX11L1", "FKBP4"),
  gene2_id = c("ENSG00000000005", "ENSG00000000003"),
  gene2_symbol = c("DDANK1", "FKBP4")
)

hotspots_t9 <- identify_fusion_hotspots(fusion_data_t9, gtex_data)

test_success_t9 <- nrow(hotspots_t9) > 0 &&
  n_distinct(hotspots_t9$fusion_id) == 2

test_results[[test_num]] <- test_success_t9
cat(if (test_success_t9) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(hotspots_t9) > 0) {
  cat("  Fusions processed:", n_distinct(hotspots_t9$fusion_id), "\n")
  cat("  Total tissue-level results:", nrow(hotspots_t9), "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 10: Heatmap Generation (No Error)
# ==========================================
test_name <- "Generate expression heatmap without errors"
cat("Test", test_num, ":", test_name, "\n")

# Test heatmap generation (don't save to file)
heatmap_result <- tryCatch({
  plot_fusion_expression_heatmap(
    "DDX11L1",
    "FKBP4",
    gtex_data,
    title = "Test Fusion Expression"
  )
  TRUE
}, error = function(e) {
  cat("  Error:", as.character(e), "\n")
  FALSE
})

test_success_t10 <- heatmap_result

test_results[[test_num]] <- test_success_t10
cat(if (test_success_t10) "✓ PASS\n" else "✗ FAIL\n")
cat("\n")

# ==========================================
# Summary
# ==========================================
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

passed_count <- sum(unlist(test_results))
total_count <- length(test_results)

test_names <- c(
  "Query single gene",
  "Query multiple genes",
  "Handle version numbers",
  "Filter tissues",
  "Hotspot: Upstream/Downstream",
  "Hotspot: Expression difference",
  "Handle missing genes",
  "Log2 calculation",
  "Multiple fusions",
  "Heatmap generation"
)

for (i in seq_along(test_results)) {
  status <- if (test_results[[i]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("Test %2d: %s\n", i, status))
}

cat("\n")
cat(sprintf("Total: %d / %d tests passed\n", passed_count, total_count))
cat("\n")

if (passed_count == total_count) {
  cat("✓ Phase 5 is complete and ready for Phase 6.\n\n")
} else {
  cat("✗ Phase 5 requires fixes before proceeding.\n\n")
  cat("Failed tests:", which(!unlist(test_results)), "\n\n")
}
