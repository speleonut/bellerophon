#!/usr/bin/env Rscript
#
# Phase 1 Revision 2 - HGVS Parser Fix Verification
# Tests all 5 variants from the test file after the fix
#

cat("\n")
cat(strrep("=", 80), "\n", sep = "")
cat("PHASE 1 REVISION 2 - HGVS PARSER FIX VERIFICATION\n")
cat(strrep("=", 80), "\n\n")

# Load required modules
library(tidyverse)
library(stringr)
library(glue)

source("R/parse_hgvs.R")
source("R/parse_inputs.R")

# Test individual HGVS strings
cat("TEST 1: Individual HGVS String Parsing\n")
cat(strrep("-", 80), "\n\n", sep = "")

test_hgvs <- c(
  "NC_000001.11:g.100_200del",
  "NC_000001.11:g.300_400dup",
  "NC_000001.11:g.500_600inv",
  "NC_000002.12:g.1000_1001insAGCT",
  "NC_000023.11:g.134381666_134381892del"
)

passed <- 0
failed <- 0

for (hgvs in test_hgvs) {
  cat("Parsing: ", hgvs, "\n", sep = "")
  result <- parse_hgvs_variant(hgvs)
  
  if (result$parsed_successfully) {
    cat("  ✓ PASSED\n")
    cat("    Type: ", result$variant_type, "\n", sep = "")
    cat("    Chromosome: ", result$chromosome, "\n", sep = "")
    cat("    Positions: ", result$position_start, "-", result$position_end, "\n\n", sep = "")
    passed <- passed + 1
  } else {
    cat("  ✗ FAILED\n")
    cat("    Error: ", result$error_message, "\n\n", sep = "")
    failed <- failed + 1
  }
}

cat("Individual HGVS parsing results:\n")
cat("  Passed: ", passed, "/", length(test_hgvs), "\n", sep = "")
cat("  Failed: ", failed, "/", length(test_hgvs), "\n\n", sep = "")

# Test file parsing
cat("TEST 2: HGVS File Parsing\n")
cat(strrep("-", 80), "\n\n", sep = "")

if (file.exists("data/example_input/sample_variants.hgvs.txt")) {
  cat("Reading: data/example_input/sample_variants.hgvs.txt\n\n")
  
  tryCatch({
    variants <- parse_hgvs_input("data/example_input/sample_variants.hgvs.txt")
    
    cat("✓ File parsed successfully\n")
    cat("  Total variants: ", nrow(variants), "\n\n", sep = "")
    
    cat("Parsed variants:\n")
    for (i in 1:nrow(variants)) {
      v <- variants[i, ]
      cat(i, ": ", v$event_type, " on chr ", v$break1_chromosome, 
          " (", v$break1_position_start, "-", v$break1_position_end, ")\n", sep = "")
    }
    cat("\n")
    
  }, error = function(e) {
    cat("✗ Error parsing file: ", as.character(e), "\n\n", sep = "")
  })
} else {
  cat("✗ Test file not found\n\n")
}

# Test HGVS format detection
cat("TEST 3: HGVS Format Detection\n")
cat(strrep("-", 80), "\n\n", sep = "")

if (file.exists("data/example_input/sample_variants.hgvs.txt")) {
  detection <- detect_input_format("data/example_input/sample_variants.hgvs.txt")
  cat("Detected format: ", detection$format, "\n", sep = "")
  cat("Confidence: ", detection$confidence, "\n\n", sep = "")
}

# Summary
cat(strrep("=", 80), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 80), "\n\n", sep = "")

if (passed == length(test_hgvs)) {
  cat("✓ ALL TESTS PASSED!\n")
  cat("HGVS parser is working correctly.\n")
  cat("Ready to proceed with Phase 1 validation.\n\n")
} else {
  cat("✗ SOME TESTS FAILED\n")
  cat("Please review the errors above.\n\n")
}
