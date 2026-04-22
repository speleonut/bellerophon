#!/usr/bin/env Rscript
#
# FINAL PHASE 1 VALIDATION - Comprehensive Test
# This runs all critical tests to ensure Phase 1 is working correctly
#

cat("\n")
cat(strrep("=", 80), "\n", sep = "")
cat("FINAL PHASE 1 VALIDATION\n")
cat(strrep("=", 80), "\n\n")

# Test suite counter
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
# LOAD MODULES
# ============================================================================

cat("Loading modules...\n")

source("R/parse_hgvs.R")
source("R/utils.R")
library(tidyverse)
library(stringr)
library(readr)
library(glue)
source("R/parse_inputs.R")

cat("✓ All modules loaded\n\n")

# ============================================================================
# RUN TESTS
# ============================================================================

# Test 1: HGVS Parser with simple deletion
run_test("HGVS Parser - Deletion", function() {
  result <- parse_hgvs_variant("NC_000001.11:g.100_200del")
  return(result$parsed_successfully && 
         result$variant_type == "deletion" &&
         result$chromosome == "1" &&
         result$position_start == 100 &&
         result$position_end == 200)
})

# Test 2: HGVS Parser with duplication
run_test("HGVS Parser - Duplication", function() {
  result <- parse_hgvs_variant("NC_000001.11:g.300_400dup")
  return(result$parsed_successfully && 
         result$variant_type == "duplication" &&
         result$chromosome == "1")
})

# Test 3: HGVS Parser with inversion
run_test("HGVS Parser - Inversion", function() {
  result <- parse_hgvs_variant("NC_000001.11:g.500_600inv")
  return(result$parsed_successfully && 
         result$variant_type == "inversion" &&
         result$chromosome == "1")
})

# Test 4: HGVS Parser - X chromosome
run_test("HGVS Parser - X Chromosome", function() {
  result <- parse_hgvs_variant("NC_000023.11:g.134381666_134381892del")
  return(result$parsed_successfully && 
         result$variant_type == "deletion" &&
         result$chromosome == "X")
})

# Test 5: Chromosome Normalization - Vector Input
run_test("Chromosome Normalization - Vector", function() {
  input_vec <- c("1", "chr1", "X", "chrX", "MT", "chrM")
  expected <- c("1", "1", "X", "X", "MT", "MT")
  result <- normalize_chromosome(input_vec, style = "NCBI")
  return(all(result == expected))
})

# Test 6: Chromosome Normalization - Single Value
run_test("Chromosome Normalization - Single", function() {
  result_chr1 <- normalize_chromosome("chr1", style = "NCBI")
  result_chrx <- normalize_chromosome("chrX", style = "NCBI")
  result_chrm <- normalize_chromosome("chrM", style = "NCBI")
  return(result_chr1 == "1" && result_chrx == "X" && result_chrm == "MT")
})

# Test 7: Uncertain Breakpoint - Forward Strand
run_test("Uncertain Breakpoint - Forward Strand", function() {
  result <- handle_uncertain_breakpoint(100, 200, "+")
  return(result == 200)  # Forward should use maximum
})

# Test 8: Uncertain Breakpoint - Reverse Strand
run_test("Uncertain Breakpoint - Reverse Strand", function() {
  result <- handle_uncertain_breakpoint(100, 200, "-")
  return(result == 100)  # Reverse should use minimum
})

# Test 9: Uncertain Breakpoint - Point Breakpoint
run_test("Uncertain Breakpoint - Point Breakpoint", function() {
  result <- handle_uncertain_breakpoint(100, 100, "+")
  return(result == 100)
})

# Test 10: Variant Type Classification - Deletion
run_test("Variant Type Classification - Deletion", function() {
  result <- classify_variant_type("DEL", "1", "1", "+", "+")
  return(result == "deletion")
})

# Test 11: Variant Type Classification - Translocation
run_test("Variant Type Classification - Translocation", function() {
  result <- classify_variant_type("BND", "1", "4", "+", "+")
  return(result == "translocation")
})

# Test 12: Variant Type Classification - Inversion (by orientation)
run_test("Variant Type Classification - Inversion", function() {
  result <- classify_variant_type("INV", "1", "1", "+", "-")
  return(result == "inversion")
})

# Test 13: Input Format Detection - VCF
run_test("Input Format Detection - VCF", function() {
  if (!file.exists("data/example_input/sample_variants.vcf")) return(FALSE)
  detection <- detect_input_format("data/example_input/sample_variants.vcf")
  return(detection$format == "vcf" && detection$confidence > 0.8)
})

# Test 14: Input Format Detection - MAVIS
run_test("Input Format Detection - MAVIS", function() {
  if (!file.exists("data/example_input/sample_variants.mavis.tsv")) return(FALSE)
  detection <- detect_input_format("data/example_input/sample_variants.mavis.tsv")
  return(detection$format == "mavis" && detection$confidence > 0.8)
})

# Test 15: Input Format Detection - HGVS
run_test("Input Format Detection - HGVS", function() {
  if (!file.exists("data/example_input/sample_variants.hgvs.txt")) return(FALSE)
  detection <- detect_input_format("data/example_input/sample_variants.hgvs.txt")
  return(detection$format == "hgvs" && detection$confidence > 0.8)
})

# Test 16: HGVS Input Parsing
run_test("HGVS Input Parsing", function() {
  if (!file.exists("data/example_input/sample_variants.hgvs.txt")) return(FALSE)
  variants <- parse_hgvs_input("data/example_input/sample_variants.hgvs.txt")
  return(nrow(variants) > 0 && 
         all(c("break1_chromosome", "break1_position_start", "event_type") %in% names(variants)))
})

# Test 17: MAVIS Input Parsing
run_test("MAVIS Input Parsing", function() {
  if (!file.exists("data/example_input/sample_variants.mavis.tsv")) return(FALSE)
  variants <- parse_mavis_input("data/example_input/sample_variants.mavis.tsv")
  return(nrow(variants) > 0 && nrow(variants) <= 10)
})

# Test 18: NCBI Contig Mapping
run_test("NCBI Contig Mapping", function() {
  chr1 <- ncbi_contig_to_chromosome("NC_000001.11")
  chrx <- ncbi_contig_to_chromosome("NC_000023.11")
  chrm <- ncbi_contig_to_chromosome("NC_012920.1")
  return(chr1 == "1" && chrx == "X" && chrm == "MT")
})

# Test 19: Chromosome Validation
run_test("Chromosome Validation", function() {
  valid_check <- check_valid_chromosomes(c("1", "X", "22", "Y", "MT"))
  return(all(valid_check$is_valid) && nrow(valid_check) == 5)
})

# Test 20: Variant ID Generation
run_test("Variant ID Generation", function() {
  ids <- generate_variant_ids(5, seed = 42)
  return(length(ids) == 5 && 
         all(nchar(ids) == 6) &&
         length(unique(ids)) == 5)  # All unique
})

# ============================================================================
# SUMMARY
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("TEST SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("Total Tests: ", total_tests, "\n", sep = "")
cat("Passed: ", passed_tests, "\n", sep = "")
cat("Failed: ", failed_tests, "\n\n", sep = "")

if (failed_tests == 0) {
  cat("✓ ALL TESTS PASSED!\n")
  cat("Phase 1 is ready for Phase 2 implementation.\n\n")
} else {
  cat("✗ SOME TESTS FAILED\n")
  cat("Please review the errors above.\n\n")
}

cat(strrep("=", 80), "\n", sep = "")
