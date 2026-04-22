#!/usr/bin/env Rscript
#
# Phase 2 Validation - Annotation Loading Tests
#

cat("\n")
cat(strrep("=", 80), "\n", sep = "")
cat("PHASE 2 VALIDATION - ANNOTATION LOADING\n")
cat(strrep("=", 80), "\n\n")

# Load required packages
library(tidyverse)
library(GenomicFeatures)
library(txdbmaker)
library(AnnotationHub)
library(biomaRt)

# Source Phase 1 and Phase 2 modules
source("R/annotations.R")

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
# PHASE 2 TESTS
# ============================================================================

# Test 1: GENCODE Loading
run_test("Load GENCODE Annotations", function() {
  cat("  Attempting to load GENCODE from rtracklayer...\n")
  gencode <- load_gencode_annotations(force_download = FALSE)
  
  if (gencode$loaded_successfully) {
    cat("  Source: rtracklayer (direct download from EBI)\n")
    
    # Check required objects
    has_txdb <- !is.null(gencode$txdb)
    has_genes <- !is.null(gencode$genes) && length(gencode$genes) > 0
    has_gene_lookup <- !is.null(gencode$gene_lookup) && nrow(gencode$gene_lookup) > 0
    has_exons <- !is.null(gencode$exons) && length(gencode$exons) > 0
    
    cat("  TxDb object: ", ifelse(has_txdb, "✓", "✗"), "\n", sep = "")
    cat("  Genes GRanges: ", ifelse(has_genes, "✓", "✗"), " (", length(gencode$genes), " genes)\n", sep = "")
    cat("  Gene lookup table: ", ifelse(has_gene_lookup, "✓", "✗"), " (", nrow(gencode$gene_lookup), " genes)\n", sep = "")
    cat("  Exons loaded: ", ifelse(has_exons, "✓", "✗"), " (", length(gencode$exons), " exons)\n", sep = "")
    
    return(has_genes && has_gene_lookup && has_exons)
  } else {
    cat("  Error: ", gencode$error, "\n", sep = "")
    return(FALSE)
  }
})

# Test 2: Disease Gene Loading
run_test("Load Disease Gene List", function() {
  disease_genes <- load_disease_genes("data/annotations/Nijmegen.DG.ENSG.list.txt")
  
  # Should have loaded at least some genes
  has_genes <- nrow(disease_genes) > 0
  
  if (has_genes) {
    cat("  Loaded ", nrow(disease_genes), " disease genes\n", sep = "")
    cat("  Sample genes: ", paste(head(disease_genes$gene_id, 3), collapse = ", "), "\n", sep = "")
  } else {
    cat("  No disease genes loaded\n")
  }
  
  return(has_genes)
})

# Test 3: GTEx Metadata Initialization
run_test("Initialize GTEx Metadata", function() {
  gtex <- initialize_gtex_metadata()
  
  if (gtex$loaded_successfully) {
    cat("  Initialized GTEx metadata\n")
    cat("  Available tissues: ", nrow(gtex$tissues), "\n", sep = "")
    
    # Check tissue structure
    has_tissue_id <- "tissue_id" %in% names(gtex$tissues)
    has_tissue_name <- "tissue_name" %in% names(gtex$tissues)
    has_category <- "tissue_category" %in% names(gtex$tissues)
    
    cat("  Tissue ID column: ", ifelse(has_tissue_id, "✓", "✗"), "\n", sep = "")
    cat("  Tissue name column: ", ifelse(has_tissue_name, "✓", "✗"), "\n", sep = "")
    cat("  Tissue category column: ", ifelse(has_category, "✓", "✗"), "\n", sep = "")
    
    # Show sample tissues
    cat("  Sample tissues:\n")
    for (i in 1:min(3, nrow(gtex$tissues))) {
      cat("    - ", gtex$tissues$tissue_name[i], "\n", sep = "")
    }
    
    return(has_tissue_id && has_tissue_name && has_category)
  } else {
    cat("  Error: ", gtex$error, "\n", sep = "")
    return(FALSE)
  }
})

# Test 4: Disease Gene File Handling
run_test("Disease Gene File - Non-existent File", function() {
  disease_genes <- load_disease_genes("data/nonexistent_file.txt")
  
  # Should return empty tibble gracefully
  return(is_tibble(disease_genes))
})

# Test 5: Tissue Mapping Consistency
run_test("GTEx Tissue ID to Name Mapping", function() {
  gtex <- initialize_gtex_metadata()
  
  if (!gtex$loaded_successfully) return(FALSE)
  
  # Check that tissue IDs and names are consistent
  tissues <- gtex$tissues
  
  # No NA values
  no_na <- !any(is.na(tissues$tissue_id)) && 
           !any(is.na(tissues$tissue_name)) &&
           !any(is.na(tissues$tissue_category))
  
  # Tissue IDs are unique
  unique_ids <- length(unique(tissues$tissue_id)) == nrow(tissues)
  
  cat("  Total tissues: ", nrow(tissues), "\n", sep = "")
  cat("  Unique tissue IDs: ", length(unique(tissues$tissue_id)), "\n", sep = "")
  cat("  No missing values: ", ifelse(no_na, "✓", "✗"), "\n", sep = "")
  cat("  All IDs unique: ", ifelse(unique_ids, "✓", "✗"), "\n", sep = "")
  
  return(no_na && unique_ids)
})

# Test 6: Data Structure Compatibility
run_test("Annotation Objects Data Structure", function() {
  disease_genes <- load_disease_genes("data/annotations/Nijmegen.DG.ENSG.list.txt")
  gtex <- initialize_gtex_metadata()
  
  # Check disease genes structure
  dg_valid <- is_tibble(disease_genes) && 
              "gene_id" %in% names(disease_genes)
  
  # Check GTEx structure
  gtex_valid <- is_tibble(gtex$tissues) &&
                "tissue_id" %in% names(gtex$tissues) &&
                "tissue_name" %in% names(gtex$tissues)
  
  cat("  Disease genes tibble: ", ifelse(dg_valid, "✓", "✗"), "\n", sep = "")
  cat("  GTEx tissues tibble: ", ifelse(gtex_valid, "✓", "✗"), "\n", sep = "")
  
  return(dg_valid && gtex_valid)
})

# ============================================================================
# SUMMARY
# ============================================================================

cat(strrep("=", 80), "\n", sep = "")
cat("PHASE 2 VALIDATION SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("Total Tests: ", total_tests, "\n", sep = "")
cat("Passed: ", passed_tests, "\n", sep = "")
cat("Failed: ", failed_tests, "\n\n", sep = "")

if (failed_tests == 0) {
  cat("✓ ALL TESTS PASSED!\n")
  cat("Phase 2 is complete and ready for Phase 3 implementation.\n\n")
} else {
  cat("✗ SOME TESTS FAILED\n")
  cat("Please review the errors above.\n\n")
}

cat(strrep("=", 80), "\n", sep = "")
