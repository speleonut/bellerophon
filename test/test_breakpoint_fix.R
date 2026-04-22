#!/usr/bin/env Rscript

# Test script to validate the breakpoint_annotation.R refactoring
# Tests that the new functions work correctly with dual breakpoints

cat("\n=== Testing breakpoint_annotation.R Refactoring ===\n\n")

# Load required libraries
library(dplyr)
library(tibble)
library(GenomicRanges)
library(IRanges)

# Source the modified module
source("R/breakpoint_annotation.R")

# Test 1: split_breakpoints_for_annotation function
cat("Test 1: split_breakpoints_for_annotation\n")
test_variants <- tibble::tibble(
  variant_id = c("30EFE4", "40D498"),
  variant_type = c("inversion", "duplication"),
  break1_chromosome = c("chr2", "chr14"),
  break1_position = c(166007750, 95239053),
  break1_orientation = c("+", "+"),
  break2_chromosome = c("chr2", "chr14"),
  break2_position = c(166077750, 95570195),
  break2_orientation = c("-", "+")
)

split_result <- split_breakpoints_for_annotation(test_variants)

cat("  BP1 result:\n")
print(split_result$bp1)

cat("\n  BP2 result:\n")
print(split_result$bp2)

cat("\n✓ Test 1 passed: Breakpoints split correctly\n\n")

# Test 2: Verify column names
cat("Test 2: Column name mapping\n")
expected_cols <- c("variant_id", "seqname", "pos", "strand", "breakpoint")
bp1_cols <- names(split_result$bp1)
bp2_cols <- names(split_result$bp2)

if (all(expected_cols %in% bp1_cols) && all(expected_cols %in% bp2_cols)) {
  cat("✓ Test 2 passed: Column names are correct\n\n")
} else {
  cat("✗ Test 2 FAILED: Column names don't match\n")
  cat("  Expected: ", paste(expected_cols, collapse=", "), "\n")
  cat("  BP1 got: ", paste(bp1_cols, collapse=", "), "\n")
  cat("  BP2 got: ", paste(bp2_cols, collapse=", "), "\n")
}

# Test 3: Verify breakpoint numbering
cat("Test 3: Breakpoint numbering\n")
if (all(split_result$bp1$breakpoint == 1) && all(split_result$bp2$breakpoint == 2)) {
  cat("✓ Test 3 passed: Breakpoints numbered correctly (1 and 2)\n\n")
} else {
  cat("✗ Test 3 FAILED: Breakpoint numbering incorrect\n")
}

# Test 4: Verify data integrity
cat("Test 4: Data integrity\n")
bp1_match <- (split_result$bp1$seqname == test_variants$break1_chromosome) %>% all()
bp1_pos_match <- (split_result$bp1$pos == test_variants$break1_position) %>% all()
bp2_match <- (split_result$bp2$seqname == test_variants$break2_chromosome) %>% all()
bp2_pos_match <- (split_result$bp2$pos == test_variants$break2_position) %>% all()

if (bp1_match && bp1_pos_match && bp2_match && bp2_pos_match) {
  cat("✓ Test 4 passed: All data mapped correctly\n\n")
} else {
  cat("✗ Test 4 FAILED: Data mapping issues\n")
  cat("  BP1 seqname match: ", bp1_match, "\n")
  cat("  BP1 pos match: ", bp1_pos_match, "\n")
  cat("  BP2 seqname match: ", bp2_match, "\n")
  cat("  BP2 pos match: ", bp2_pos_match, "\n")
}

cat("=== All tests completed ===\n\n")
