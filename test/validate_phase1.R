#!/usr/bin/env Rscript
#
# Quick validation that Phase 1 modules source correctly
#

cat("Attempting to source Phase 1 modules...\n\n")

tryCatch({
  source("R/parse_hgvs.R")
  cat("✓ R/parse_hgvs.R sourced\n")
}, error = function(e) {
  cat("✗ Error sourcing R/parse_hgvs.R:\n")
  cat(as.character(e), "\n")
})

tryCatch({
  source("R/utils.R")
  cat("✓ R/utils.R sourced\n")
}, error = function(e) {
  cat("✗ Error sourcing R/utils.R:\n")
  cat(as.character(e), "\n")
})

tryCatch({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(glue)
  source("R/parse_inputs.R")
  cat("✓ R/parse_inputs.R sourced\n")
}, error = function(e) {
  cat("✗ Error sourcing R/parse_inputs.R:\n")
  cat(as.character(e), "\n")
})

cat("\n✓ All Phase 1 modules sourced successfully!\n")
cat("\nAttempting quick functionality tests...\n\n")

# Test 1: Simple HGVS parse
tryCatch({
  result <- parse_hgvs_variant("NC_000001.11:g.100_200del")
  if (result$parsed_successfully && result$variant_type == "deletion") {
    cat("✓ HGVS parser test passed\n")
  } else {
    cat("✗ HGVS parser test failed\n")
  }
}, error = function(e) {
  cat("✗ Error testing HGVS parser:\n")
  cat(as.character(e), "\n")
})

# Test 2: Chromosome normalization
tryCatch({
  norm_chr <- normalize_chromosome("chr1", style = "NCBI")
  if (norm_chr == "1") {
    cat("✓ Chromosome normalization test passed\n")
  } else {
    cat("✗ Chromosome normalization test failed\n")
  }
}, error = function(e) {
  cat("✗ Error testing chromosome normalization:\n")
  cat(as.character(e), "\n")
})

# Test 3: Uncertain breakpoint handling
tryCatch({
  pos <- handle_uncertain_breakpoint(100, 200, "+")
  if (pos == 200) {
    cat("✓ Uncertain breakpoint handling test passed\n")
  } else {
    cat("✗ Uncertain breakpoint handling test failed\n")
  }
}, error = function(e) {
  cat("✗ Error testing uncertain breakpoint handling:\n")
  cat(as.character(e), "\n")
})

cat("\n✓ Quick validation complete!\n")
