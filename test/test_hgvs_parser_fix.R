#!/usr/bin/env Rscript
#
# Quick HGVS Parser Test After Fix
#

cat("Testing HGVS Parser after fix...\n\n")

source("R/parse_hgvs.R")

test_cases <- c(
  "NC_000001.11:g.100_200del",
  "NC_000001.11:g.300_400dup",
  "NC_000001.11:g.500_600inv",
  "NC_000002.12:g.1000_1100ins",
  "NC_000023.11:g.134381666_134381892del"
)

for (hgvs_str in test_cases) {
  cat("Testing: ", hgvs_str, "\n", sep = "")
  result <- parse_hgvs_variant(hgvs_str)
  
  if (result$parsed_successfully) {
    cat("  ✓ Parsed successfully\n")
    cat("    Type: ", result$variant_type, "\n", sep = "")
    cat("    Chromosome: ", result$chromosome, "\n", sep = "")
    cat("    Position: ", result$position_start, "-", result$position_end, "\n\n", sep = "")
  } else {
    cat("  ✗ Failed to parse\n")
    cat("    Error: ", result$error_message, "\n\n", sep = "")
  }
}

cat("Test complete!\n")
