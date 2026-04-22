#!/usr/bin/env Rscript
#
# Phase 1 Testing Script
# Tests: HGVS parser, chromosome normalization, input detection
#

library(tidyverse)
library(stringr)
library(stringi)
library(glue)
library(readr)

# Source Phase 1 modules
source("R/parse_hgvs.R")
source("R/utils.R")
source("R/parse_inputs.R")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("PHASE 1 TEST SUITE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Test 1: HGVS Parser
cat("TEST 1: HGVS Parser\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

test_hgvs_cases <- c(
  "NC_000001.11:g.100_200del",
  "NC_000001.11:g.300_400dup",
  "NC_000001.11:g.500_600inv",
  "NC_000023.11:g.134381666_134381892del"
)

for (hgvs_str in test_hgvs_cases) {
  result <- parse_hgvs_variant(hgvs_str)
  cat(glue("Input: {hgvs_str}\n"))
  cat(glue("  Parsed: {result$parsed_successfully}\n"))
  if (result$parsed_successfully) {
    cat(glue("  Type: {result$variant_type}\n"))
    cat(glue("  Chromosome: {result$chromosome}\n"))
    cat(glue("  Positions: {result$position_start}-{result$position_end}\n"))
  } else {
    cat(glue("  Error: {result$error_message}\n"))
  }
  cat("\n")
}

# Test 2: Chromosome Normalization
cat("TEST 2: Chromosome Normalization\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

test_chr_cases <- c("1", "chr1", "X", "chrX", "MT", "chrM", "Y", "chrY")
for (chr in test_chr_cases) {
  normalized <- normalize_chromosome(chr, style = "NCBI")
  cat(glue("{chr} -> {normalized}\n"))
}
cat("\n")

# Test 3: Uncertain Breakpoint Handling
cat("TEST 3: Uncertain Breakpoint Handling\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

test_breakpoint_cases <- list(
  list(start = 100, end = 100, orient = "+", desc = "Point breakpoint forward"),
  list(start = 100, end = 200, orient = "+", desc = "Range breakpoint forward"),
  list(start = 100, end = 200, orient = "-", desc = "Range breakpoint reverse"),
  list(start = 100, end = 100, orient = "-", desc = "Point breakpoint reverse")
)

for (case in test_breakpoint_cases) {
  selected <- handle_uncertain_breakpoint(case$start, case$end, case$orient)
  cat(glue("{case$desc}: ({case$start}-{case$end}, {case$orient}) -> {selected}\n"))
}
cat("\n")

# Test 4: Input Format Detection
cat("TEST 4: Input Format Detection\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

test_input_files <- c(
  "data/example_input/sample_variants.vcf",
  "data/example_input/sample_variants.mavis.tsv",
  "data/example_input/sample_variants.hgvs.txt"
)

for (input_file in test_input_files) {
  if (file.exists(input_file)) {
    detection <- detect_input_format(input_file)
    cat(glue("File: {input_file}\n"))
    cat(glue("  Detected format: {detection$format}\n"))
    cat(glue("  Confidence: {detection$confidence}\n"))
    cat("\n")
  } else {
    cat(glue("WARNING: File not found: {input_file}\n\n"))
  }
}

# Test 5: Variant Type Classification
cat("TEST 5: Variant Type Classification\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

test_classification_cases <- list(
  list(et = "DEL", b1c = "1", b2c = "1", b1o = "+", b2o = "+", desc = "Deletion"),
  list(et = "DUP", b1c = "1", b2c = "1", b1o = "+", b2o = "+", desc = "Duplication"),
  list(et = "INV", b1c = "1", b2c = "1", b1o = "+", b2o = "-", desc = "Inversion"),
  list(et = "BND", b1c = "1", b2c = "4", b1o = "+", b2o = "+", desc = "Translocation"),
  list(et = NA, b1c = "1", b2c = "1", b1o = "+", b2o = "+", desc = "Unknown type")
)

for (case in test_classification_cases) {
  classified <- classify_variant_type(case$et, case$b1c, case$b2c, case$b1o, case$b2o)
  cat(glue("{case$desc}: {classified}\n"))
}
cat("\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("PHASE 1 TEST COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
