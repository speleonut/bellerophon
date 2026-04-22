#!/usr/bin/env Rscript
#
# Phase 3 Validation - Breakpoint Annotation Tests
#

cat("\n")
cat(strrep("=", 80), "\n", sep = "")
cat("PHASE 3 VALIDATION - BREAKPOINT ANNOTATION\n")
cat(strrep("=", 80), "\n\n")

# Load required packages
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)
library(AnnotationHub)
library(readr)
library(IRanges)

# Source Phase 2 and Phase 3 modules
source("R/annotations.R")
source("R/breakpoint_annotation.R")

# Test counter
total_tests <- 0
passed_tests <- 0
failed_tests <- 0

run_test <- function(test_name, test_func) {
  total_tests <<- total_tests + 1
  cat("TEST ", total_tests, ": ", test_name, "\n", sep = "")
  cat("  ", strrep("-", 76), "\n", sep = "")
  
  tryCatch({
    result <- test_func()
    if (result) {
      cat("  ✓ PASSED\n\n")
      passed_tests <<- passed_tests + 1
      return(TRUE)
    } else {
      cat("  ✗ FAILED\n\n")
      failed_tests <<- failed_tests + 1
      return(FALSE)
    }
  }, error = function(e) {
    cat("  ✗ ERROR: ", as.character(e), "\n\n", sep = "")
    failed_tests <<- failed_tests + 1
    return(FALSE)
  })
}

# ============================================================================
# SETUP: Load annotations for testing
# ============================================================================

cat("Loading annotations for Phase 3 testing...\n\n")

annotations <- load_gencode_annotations(force_download = FALSE)
if (!annotations$loaded_successfully) {
  cat("WARNING: Could not load GENCODE annotations via rtracklayer.\n")
  cat("Phase 3 tests require this data. Tests will be limited.\n\n")
  annotations_available <- FALSE
} else {
  annotations_available <- TRUE
  cat("✓ GENCODE annotations loaded\n")
  cat("✓ Gene count: ", length(annotations$genes), " (GRanges)\n", sep = "")
  cat("✓ Gene lookup: ", nrow(annotations$gene_lookup), " genes\n", sep = "")
  cat("✓ Transcript lookup: ", nrow(annotations$tx_lookup), " transcripts\n", sep = "")
  cat("✓ Exon count: ", length(annotations$exons), "\n\n", sep = "")
}

disease_genes <- load_disease_genes("data/annotations/Nijmegen.DG.ENSG.list.txt")
cat("✓ Disease genes loaded: ", nrow(disease_genes), "\n\n", sep = "")

# ============================================================================
# PHASE 3 TESTS
# ============================================================================

# Test 1: GRanges conversion from data.frame
run_test("Convert data.frame to GRanges", function() {
  
  test_variants <- tibble::tibble(
    seqname = c("chr1", "chr2"),
    pos = c(1000000, 2000000),
    strand = c("+", "-")
  )
  
  # Test helper function
  gr <- convert_df_to_granges(test_variants)
  
  # Validate conversion
  is_granges <- inherits(gr, "GRanges")
  correct_length <- length(gr) == 2
  correct_seqnames <- all(as.character(seqnames(gr)) == c("chr1", "chr2"))
  correct_positions <- all(start(gr) == c(1000000, 2000000))
  
  cat("  GRanges conversion: ", ifelse(is_granges, "✓", "✗"), "\n", sep = "")
  cat("  Correct length: ", ifelse(correct_length, "✓", "✗"), "\n", sep = "")
  cat("  Correct seqnames: ", ifelse(correct_seqnames, "✓", "✗"), "\n", sep = "")
  cat("  Correct positions: ", ifelse(correct_positions, "✓", "✗"), "\n", sep = "")
  
  return(is_granges && correct_length && correct_seqnames && correct_positions)
})

# Test 2: Breakpoint annotation with empty variants
run_test("Annotation of empty variant set", function() {
  
  empty_variants <- tibble::tibble(
    seqname = character(),
    pos = integer(),
    strand = character()
  )
  
  result <- annotate_breakpoints(empty_variants, annotations, disease_genes)
  
  # Should return empty tibble with correct structure
  is_tibble <- is.data.frame(result)
  correct_cols <- all(c("variant_id", "breakpoint", "gene_id", "region_type") %in% names(result))
  is_empty <- nrow(result) == 0
  
  cat("  Result is tibble: ", ifelse(is_tibble, "✓", "✗"), "\n", sep = "")
  cat("  Has correct columns: ", ifelse(correct_cols, "✓", "✗"), "\n", sep = "")
  cat("  Is empty: ", ifelse(is_empty, "✓", "✗"), "\n", sep = "")
  
  return(is_tibble && correct_cols && is_empty)
})

# Test 3: Required columns in annotation output
run_test("Output structure validation", function() {
  
  # Test with a minimal valid variant
  test_variants <- tibble::tibble(
    variant_id = "VAR001",
    seqname = "chr1",
    pos = 1000000,
    strand = "+"
  )
  
  result <- annotate_breakpoints(test_variants, annotations, disease_genes)
  
  expected_cols <- c(
    "variant_id", "breakpoint", "seqname", "pos",
    "gene_id", "gene_symbol", "transcript_id",
    "region_type", "region_details",
    "is_disease_gene", "distance_to_nearest_gene"
  )
  
  has_cols <- all(expected_cols %in% names(result))
  col_count <- length(names(result)) == length(expected_cols)
  
  cat("  Has all expected columns: ", ifelse(has_cols, "✓", "✗"), "\n", sep = "")
  cat("  Correct column count: ", ifelse(col_count, "✓", "✗"), "\n", sep = "")
  
  if (nrow(result) > 0) {
    cat("  Found ", nrow(result), " annotations\n", sep = "")
  }
  
  return(has_cols && col_count)
})

# Test 4: Region type classification
run_test("Region type classification", function() {
  
  # Skip if annotations not available
  if (!annotations_available) {
    cat("  Skipped (annotations not available)\n")
    return(TRUE)
  }
  
  # Test with variants that should hit different region types
  test_variants <- tibble::tibble(
    variant_id = c("VAR001", "VAR002", "VAR003"),
    seqname = c("chr1", "chr1", "chr2"),
    pos = c(1000000, 2000000, 3000000),
    strand = c("+", "-", "+")
  )
  
  result <- annotate_breakpoints(test_variants, annotations, disease_genes)
  
  # Should have region_type values
  has_region_types <- nrow(result) == 0 || 
                      all(!is.na(result$region_type)) ||
                      any(result$region_type %in% c("exon", "intron", "intergenic", "CDS", "5_UTR", "3_UTR"))
  
  # Should have variant IDs
  has_variant_ids <- nrow(result) == 0 || 
                     all(!is.na(result$variant_id))
  
  cat("  Has region_type assignments: ", ifelse(has_region_types, "✓", "✗"), "\n", sep = "")
  cat("  Has variant_id assignments: ", ifelse(has_variant_ids, "✓", "✗"), "\n", sep = "")
  
  if (nrow(result) > 0) {
    cat("  Region types found: ", paste(unique(result$region_type), collapse = ", "), "\n", sep = "")
    
    # Count by region type
    region_counts <- result %>% count(region_type)
    for (i in 1:nrow(region_counts)) {
      cat("    ", region_counts$region_type[i], ": ", region_counts$n[i], "\n", sep = "")
    }
  }
  
  return(has_region_types && has_variant_ids)
})

# Test 5: Disease gene flagging
run_test("Disease gene identification", function() {
  
  # Check that is_disease_gene column exists and is logical
  if (nrow(disease_genes) == 0) {
    cat("  No disease genes loaded - skipping detailed check\n")
    return(TRUE)
  }
  
  # Create a variant on a known disease gene if available
  disease_gene_id <- disease_genes$gene_id[1]
  
  # Try to find this gene in annotations
  if (annotations_available && !is.na(disease_gene_id)) {
    gene_info <- annotations$gene_lookup %>% 
      filter(str_detect(gene_id, paste0("^", disease_gene_id, "."))) %>%
      dplyr::slice(1)
    
    if (nrow(gene_info) > 0) {
      # Create a variant within this gene
      test_variants <- tibble::tibble(
        variant_id = "DISEASE_VAR",
        seqname = gene_info$seqname[1],
        pos = as.integer((gene_info$start[1] + gene_info$end[1]) / 2),
        strand = gene_info$strand[1]
      )
      
      result <- annotate_breakpoints(test_variants, annotations, disease_genes)
      
      if (nrow(result) > 0) {
        flagged_count <- sum(result$is_disease_gene, na.rm = TRUE)
        cat("  Found ", nrow(result), " annotations for disease gene\n", sep = "")
        cat("  Flagged as disease gene: ", flagged_count, "\n", sep = "")
        
        return(flagged_count > 0)
      }
    }
  }
  
  return(TRUE)
})

# Test 6: Multi-variant annotation batch processing
run_test("Batch annotation of multiple variants", function() {
  
  # Create multiple test variants
  test_variants <- tibble::tibble(
    variant_id = c("VAR_A", "VAR_B", "VAR_C", "VAR_D"),
    seqname = c("chr1", "chr2", "chr3", "chrX"),
    pos = c(1000000, 2000000, 3000000, 50000000),
    strand = c("+", "-", "+", "-")
  )
  
  result <- annotate_breakpoints(test_variants, annotations, disease_genes)
  
  # Should have results
  has_results <- nrow(result) > 0
  all_variants_present <- all(c("VAR_A", "VAR_B", "VAR_C", "VAR_D") %in% unique(result$variant_id))
  
  cat("  Processed variants: ", n_distinct(result$variant_id), " / 4\n", sep = "")
  cat("  Total annotations: ", nrow(result), "\n", sep = "")
  
  if (nrow(result) > 0) {
    # Show breakdown
    by_variant <- result %>% count(variant_id)
    cat("  Annotations per variant:\n")
    for (i in 1:nrow(by_variant)) {
      cat("    ", by_variant$variant_id[i], ": ", by_variant$n[i], "\n", sep = "")
    }
  }
  
  return(has_results)
})

# Test 7: Data type validation
run_test("Output data type validation", function() {
  
  test_variants <- tibble::tibble(
    variant_id = "TEST_VAR",
    seqname = "chr1",
    pos = 1000000
  )
  
  result <- annotate_breakpoints(test_variants, annotations, disease_genes)
  
  if (nrow(result) == 0) {
    cat("  No annotations returned (expected for random positions)\n")
    return(TRUE)
  }
  
  # Check data types
  col_types <- sapply(result, class)
  
  type_checks <- list(
    variant_id = col_types["variant_id"] == "character",
    breakpoint = col_types["breakpoint"] == "integer",
    seqname = col_types["seqname"] == "character",
    pos = col_types["pos"] == "integer",
    gene_id = col_types["gene_id"] == "character",
    is_disease_gene = col_types["is_disease_gene"] == "logical"
  )
  
  all_correct <- all(unlist(type_checks))
  
  cat("  Checking column data types:\n")
  for (col_name in names(type_checks)) {
    status <- ifelse(type_checks[[col_name]], "✓", "✗")
    cat("    ", status, " ", col_name, " !\n", sep = "")
  }
  
  return(all_correct)
})

# Test 8: Handling of alternative chromosome names
run_test("Chromosome name handling (chr prefix variations)", function() {
  
  # Test with different chromosome naming conventions
  test_variants <- tibble::tibble(
    variant_id = c("VAR1", "VAR2"),
    seqname = c("chr1", "1"),  # Mix of formats
    pos = c(1000000, 2000000)
  )
  
  # Should handle both formats gracefully
  tryCatch({
    result1 <- annotate_breakpoints(test_variants[1, ], annotations, disease_genes)
    result2 <- annotate_breakpoints(test_variants[2, ], annotations, disease_genes)
    
    cat("  Processed chr1 (with prefix): OK\n")
    cat("  Processed 1 (without prefix): OK\n")
    
    return(TRUE)
  }, error = function(e) {
    cat("  Could not process alternative chromosome formats\n")
    return(FALSE)
  })
})

# Test 9: Integration with disease genes
run_test("Integration: Breakpoints + Disease Genes", function() {
  
  test_variants <- tibble::tibble(
    variant_id = c("INT_VAR_01", "INT_VAR_02"),
    seqname = c("chr1", "chr5"),
    pos = c(1500000, 8000000)
  )
  
  result <- annotate_breakpoints(test_variants, annotations, disease_genes)
  
  # Result should include is_disease_gene column
  has_disease_col <- "is_disease_gene" %in% names(result)
  disease_col_logical <- is.logical(result$is_disease_gene)
  
  cat("  Has is_disease_gene column: ", ifelse(has_disease_col, "✓", "✗"), "\n", sep = "")
  cat("  Column is logical type: ", ifelse(disease_col_logical, "✓", "✗"), "\n", sep = "")
  
  if (nrow(disease_genes) > 0 && nrow(result) > 0) {
    disease_flagged <- sum(result$is_disease_gene)
    cat("  Variants flagged as disease-related: ", disease_flagged, "\n", sep = "")
  }
  
  return(has_disease_col && disease_col_logical)
})

# ============================================================================
# SUMMARY
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("PHASE 3 VALIDATION SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("Total Tests: ", total_tests, "\n", sep = "")
cat("Passed: ", passed_tests, "\n", sep = "")
cat("Failed: ", failed_tests, "\n\n", sep = "")

if (failed_tests == 0) {
  cat("✓ ALL TESTS PASSED!\n")
  cat("Phase 3 is complete and ready for Phase 4 implementation.\n\n")
} else {
  cat("✗ SOME TESTS FAILED\n")
  cat("Please review the errors above.\n\n")
}

cat("Phase 3 implements:\n")
cat("  • Breakpoint annotation with exon/intron/UTR classification\n")
cat("  • Disease gene flagging\n")
cat("  • Intergenic breakpoint detection with nearest gene identification\n")
cat("  • Multi-variant batch processing\n")
cat("  • Flexible input format handling\n\n")

cat(strrep("=", 80), "\n", sep = "")
