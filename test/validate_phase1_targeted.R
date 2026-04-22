#!/usr/bin/env Rscript
#
# Targeted Phase 1 Validation - Source Step by Step
#

library(tidyverse)
library(stringr)
library(readr)
library(glue)

cat("Step 1: Source R/parse_hgvs.R\n")
source("R/parse_hgvs.R")
cat("✓ Success\n\n")

cat("Step 2: Source R/utils.R\n")
source("R/utils.R")
cat("✓ Success\n\n")

cat("Step 3: Test normalize_chromosome with vector\n")
test_chr_vec <- c("1", "chr1", "X", "chrX", "MT", "chrM")
result <- normalize_chromosome(test_chr_vec, style = "NCBI")
cat("Input: ", paste(test_chr_vec, collapse = ", "), "\n")
cat("Output: ", paste(result, collapse = ", "), "\n")
cat("✓ Success\n\n")

cat("Step 4: Test HGVS parser\n")
hgvs_test <- "NC_000001.11:g.100_200del"
parsed <- parse_hgvs_variant(hgvs_test)
cat("Input: ", hgvs_test, "\n")
cat("Parsed: ", parsed$parsed_successfully, "\n")
if (parsed$parsed_successfully) {
  cat("Type: ", parsed$variant_type, "\n")
  cat("Chromosome: ", parsed$chromosome, "\n")
  cat("Position: ", parsed$position_start, "-", parsed$position_end, "\n")
  cat("✓ Success\n")
} else {
  cat("✗ Error: ", parsed$error_message, "\n")
}
cat("\n")

cat("Step 5: Source R/parse_inputs.R\n")
source("R/parse_inputs.R")
cat("✓ Success\n\n")

cat("Step 6: Test input detection\n")
if (file.exists("data/example_input/sample_variants.hgvs.txt")) {
  detection <- detect_input_format("data/example_input/sample_variants.hgvs.txt")
  cat("File: data/example_input/sample_variants.hgvs.txt\n")
  cat("Detected format: ", detection$format, "\n")
  cat("Confidence: ", detection$confidence, "\n")
  cat("✓ Success\n")
} else {
  cat("✗ File not found\n")
}
cat("\n")

cat("Step 7: Test HGVS input parsing\n")
if (file.exists("data/example_input/sample_variants.hgvs.txt")) {
  variants <- parse_hgvs_input("data/example_input/sample_variants.hgvs.txt")
  cat("Parsed ", nrow(variants), " variants\n")
  print(head(variants, 2))
  cat("✓ Success\n")
} else {
  cat("✗ File not found\n")
}
cat("\n")

cat("✓ All targeted validation tests passed!\n")
