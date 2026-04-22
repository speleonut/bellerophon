#!/usr/bin/env Rscript
#
# Phase 1 Setup and Dependency Installer
# Run this script first to ensure all required packages are installed
#

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("BELLEROPHON PHASE 1 SETUP\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# List of required packages
required_cran_packages <- c(
  "tidyverse",      # Core data manipulation (includes dplyr, tibble, stringr)
  "stringr",        # String manipulation
  "stringi",        # Advanced string processing
  "glue",           # String interpolation
  "readr",          # File reading
  "R.utils"         # Utilities including gunzip
)

required_bioc_packages <- c(
  "VariantAnnotation",
  "GenomicRanges",
  "GenomicFeatures",
  "biomaRt",
  "AnnotationHub",
  "Biostrings",
  "BSgenome",
  "SummarizedExperiment",
  "StructuralVariantAnnotation"
)

cat("Checking CRAN packages...\n\n")

for (pkg in required_cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("Installing ", pkg, "...\n", sep = "")
    install.packages(pkg, quiet = TRUE)
    if (require(pkg, character.only = TRUE)) {
      cat("✓ ", pkg, " installed\n\n", sep = "")
    } else {
      cat("✗ Failed to install ", pkg, "\n\n", sep = "")
    }
  } else {
    cat("✓ ", pkg, " already installed\n\n", sep = "")
  }
}

cat("Checking Bioconductor packages...\n\n")
cat("(Note: Bioconductor packages require BiocManager)\n\n")

# Check if BiocManager is available
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", quiet = TRUE)
}

for (pkg in required_bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("Installing ", pkg, " via BiocManager...\n", sep = "")
    BiocManager::install(pkg, quiet = TRUE)
    if (require(pkg, character.only = TRUE)) {
      cat("✓ ", pkg, " installed\n\n", sep = "")
    } else {
      cat("✗ Failed to install ", pkg, "\n\n", sep = "")
    }
  } else {
    cat("✓ ", pkg, " already installed\n\n", sep = "")
  }
}

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SETUP COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Next steps:\n")
cat("1. Run validate_phase1_targeted.R to verify Phase 1 implementation\n")
cat("2. Or open fusion_pipeline.qmd in RStudio and click 'Render'\n")
cat("3. Or run: source('test_phase1.R')\n\n")
