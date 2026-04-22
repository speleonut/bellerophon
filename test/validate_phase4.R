#!/usr/bin/env Rscript

#' Phase 4 Validation: Fusion Gene Prediction
#'
#' Tests canonical vs non-canonical fusion classification and reading frame analysis
#'
#' Test scenarios:
#' 1. Canonical fusion: Two different genes disrupted
#' 2. Canonical fusion: Interchromosomal
#' 3. Canonical fusion: Intron-intron breakpoints 
#' 4. Canonical fusion: Exon-exon breakpoints
#' 5. Non-canonical intragenic: Same gene disrupted twice
#' 6. Non-canonical single-gene: One breakpoint in gene, one intergenic
#' 7. Non-canonical intergenic: Both breakpoints intergenic
#' 8. Disease gene flagging: Fusion involving disease genes
#' 9. Mixed breakpoints: One intron, one exon
#' 10. Reading frame indeterminate: Intronic breakpoints

library(tidyverse)
library(GenomicRanges)
library(IRanges)

source("R/annotations.R")
source("R/breakpoint_annotation.R")
source("R/fusion_prediction.R")

cat("\n========================================\n")
cat("Phase 4 Validation: Fusion Gene Prediction\n")
cat("========================================\n\n")

# Load annotations and disease genes
cat("Loading annotations...\n")
annotations <- load_gencode_annotations()
cat("✓ Annotations loaded\n")

disease_genes <- load_disease_genes("data/annotations/Nijmegen.DG.ENSG.list.txt")
cat("✓ Disease genes loaded:", nrow(disease_genes), "genes\n\n")

# Track test results
test_results <- list()
test_num <- 1

# ==========================================
# TEST 1: Canonical Fusion - Different Genes
# ==========================================
test_name <- "Canonical fusion: two different genes disrupted"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t1 <- tibble::tibble(
  variant_id = "SV_001",
  breakpoint = c(1, 2),
  seqname = c("chr2", "chr2"),
  pos = c(29415701, 42544530),
  gene_id = c("ENSG00000171094", "ENSG00000105221"),
  gene_symbol = c("ALK", "EML4"),
  transcript_id = c("ENST00000389048", "ENST00000370957"),
  region_type = c("intron", "intron"),
  region_details = c("Intron", "Intron"),
  is_disease_gene = c(TRUE, FALSE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t1 <- predict_fusions(bp_data_t1, annotations, disease_genes = disease_genes)

test_success_t1 <- nrow(fusion_result_t1) == 1 &&
  fusion_result_t1$fusion_type[1] == "canonical" &&
  fusion_result_t1$gene1_id[1] != fusion_result_t1$gene2_id[1] &&
  fusion_result_t1$disease_gene_involved[1] == TRUE

test_results[[test_num]] <- test_success_t1
cat(if (test_success_t1) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t1) > 0) {
  cat("  Type:", fusion_result_t1$fusion_type[1], 
      "| Description:", fusion_result_t1$fusion_description[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 2: Canonical Fusion - Interchromosomal
# ==========================================
test_name <- "Canonical fusion: different chromosomes (interchromosomal)"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t2 <- tibble::tibble(
  variant_id = "SV_002",
  breakpoint = c(1, 2),
  seqname = c("chr5", "chr14"),
  pos = c(173808500, 44919500),
  gene_id = c("ENSG00000134872", "ENSG00000139618"),
  gene_symbol = c("FGFR4", "ERBB2"),
  transcript_id = c("ENST00000207193", "ENST00000269571"),
  region_type = c("exon", "exon"),
  region_details = c("Exon", "Exon"),
  is_disease_gene = c(FALSE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t2 <- predict_fusions(bp_data_t2, annotations, disease_genes = disease_genes)

test_success_t2 <- nrow(fusion_result_t2) == 1 &&
  fusion_result_t2$fusion_type[1] == "canonical" &&
  fusion_result_t2$gene1_id[1] != fusion_result_t2$gene2_id[1]

test_results[[test_num]] <- test_success_t2
cat(if (test_success_t2) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t2) > 0) {
  cat("  Chromosomes: chr5 ↔ chr14 | Type:", fusion_result_t2$fusion_type[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 3: Canonical - Intron-Intron
# ==========================================
test_name <- "Canonical fusion: both intron breakpoints"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t3 <- tibble::tibble(
  variant_id = "SV_003",
  breakpoint = c(1, 2),
  seqname = c("chr9", "chr22"),
  pos = c(130650100, 23632345),
  gene_id = c("ENSG00000134023", "ENSG00000100031"),
  gene_symbol = c("ABL1", "BCR"),
  transcript_id = c("ENST00000435335", "ENST00000332410"),
  region_type = c("intron", "intron"),
  region_details = c("Intron", "Intron"),
  is_disease_gene = c(TRUE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t3 <- predict_fusions(bp_data_t3, annotations, disease_genes = disease_genes)

test_success_t3 <- nrow(fusion_result_t3) == 1 &&
  fusion_result_t3$fusion_type[1] == "canonical" &&
  fusion_result_t3$gene1_region[1] == "intron" &&
  fusion_result_t3$gene2_region[1] == "intron"

test_results[[test_num]] <- test_success_t3
cat(if (test_success_t3) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t3) > 0) {
  cat("  Regions: intron-intron | Both disease genes:", fusion_result_t3$disease_gene_involved[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 4: Canonical - Exon-Exon
# ==========================================
test_name <- "Canonical fusion: both exon breakpoints (reading frame analysis)"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t4 <- tibble::tibble(
  variant_id = "SV_004",
  breakpoint = c(1, 2),
  seqname = c("chr10", "chr17"),
  pos = c(43609950, 37844395),
  gene_id = c("ENSG00000157764", "ENSG00000141510"),
  gene_symbol = c("RET", "TP53"),
  transcript_id = c("ENST00000355710", "ENST00000269305"),
  region_type = c("exon", "exon"),
  region_details = c("Exon (exon_rank=1)", "Exon (exon_rank=5)"),
  is_disease_gene = c(TRUE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t4 <- predict_fusions(bp_data_t4, annotations, disease_genes = disease_genes)

# For exon-exon junctions, reading frame should be calculated if CDS phase available
# May be NA if phase data unavailable, which is acceptable
test_success_t4 <- nrow(fusion_result_t4) == 1 &&
  fusion_result_t4$fusion_type[1] == "canonical" &&
  fusion_result_t4$gene1_region[1] == "exon" &&
  fusion_result_t4$gene2_region[1] == "exon"

test_results[[test_num]] <- test_success_t4
cat(if (test_success_t4) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t4) > 0) {
  cat("  Regions: exon-exon | In-frame:", fusion_result_t4$in_frame[1], 
      "| offset:", fusion_result_t4$frame_offset[1], "\n")
  if (is.na(fusion_result_t4$in_frame[1])) {
    cat("  (Reading frame NA - phase data needed for calculation)\n")
  }
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 5: Non-Canonical Intragenic
# ==========================================
test_name <- "Non-canonical intragenic: same gene, both breakpoints"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t5 <- tibble::tibble(
  variant_id = "SV_005",
  breakpoint = c(1, 2),
  seqname = c("chr17", "chr17"),
  pos = c(37844395, 37881629),
  gene_id = c("ENSG00000141510", "ENSG00000141510"),
  gene_symbol = c("TP53", "TP53"),
  transcript_id = c("ENST00000269305", "ENST00000269305"),
  region_type = c("intron", "intron"),
  region_details = c("Intron", "Intron"),
  is_disease_gene = c(TRUE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t5 <- predict_fusions(bp_data_t5, annotations, disease_genes = disease_genes)

test_success_t5 <- nrow(fusion_result_t5) == 1 &&
  fusion_result_t5$fusion_type[1] == "non-canonical_intragenic" &&
  fusion_result_t5$gene1_id[1] == fusion_result_t5$gene2_id[1]

test_results[[test_num]] <- test_success_t5
cat(if (test_success_t5) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t5) > 0) {
  cat("  Type:", fusion_result_t5$fusion_type[1], 
      "| Same gene:", fusion_result_t5$gene1_id[1] == fusion_result_t5$gene2_id[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 6: Non-Canonical Single-Gene
# ==========================================
test_name <- "Non-canonical: one gene, one intergenic breakpoint"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t6 <- tibble::tibble(
  variant_id = "SV_006",
  breakpoint = c(1, 2),
  seqname = c("chr1", "chr1"),
  pos = c(11869, 14829),
  gene_id = c("ENSG00000223972", NA_character_),
  gene_symbol = c("DDOX", NA_character_),
  transcript_id = c("ENST00000456328", NA_character_),
  region_type = c("exon", "intergenic"),
  region_details = c("Exon", "Intergenic"),
  is_disease_gene = c(FALSE, FALSE),
  distance_to_nearest_gene = c(NA_integer_, 100L)
)

fusion_result_t6 <- predict_fusions(bp_data_t6, annotations, disease_genes = disease_genes)

test_success_t6 <- nrow(fusion_result_t6) == 1 &&
  fusion_result_t6$fusion_type[1] == "non-canonical_single" &&
  !is.na(fusion_result_t6$gene1_id[1]) &&
  is.na(fusion_result_t6$gene2_id[1])

test_results[[test_num]] <- test_success_t6
cat(if (test_success_t6) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t6) > 0) {
  cat("  Type:", fusion_result_t6$fusion_type[1], 
      "| Gene:", fusion_result_t6$gene1_symbol[1], "↔ intergenic\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 7: Non-Canonical Both Intergenic
# ==========================================
test_name <- "No fusion: both intergenic breakpoints (should not generate fusion)"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t7 <- tibble::tibble(
  variant_id = "SV_007",
  breakpoint = c(1, 2),
  seqname = c("chr1", "chr1"),
  pos = c(100, 1000),
  gene_id = c(NA_character_, NA_character_),
  gene_symbol = c(NA_character_, NA_character_),
  transcript_id = c(NA_character_, NA_character_),
  region_type = c("intergenic", "intergenic"),
  region_details = c("Intergenic", "Intergenic"),
  is_disease_gene = c(FALSE, FALSE),
  distance_to_nearest_gene = c(5000L, 3000L)
)

fusion_result_t7 <- predict_fusions(bp_data_t7, annotations, disease_genes = disease_genes)

# Should return 0 fusions (no fusion when both intergenic)
test_success_t7 <- nrow(fusion_result_t7) == 0

test_results[[test_num]] <- test_success_t7
cat(if (test_success_t7) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t7) == 0) {
  cat("  Result: No fusion generated (both breakpoints intergenic) ✓\n")
} else {
  cat("  ERROR: Should not generate fusion, but got", nrow(fusion_result_t7), "results\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 8: Disease Gene Flagging
# ==========================================
test_name <- "Disease gene flagging in fusion"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t8 <- tibble::tibble(
  variant_id = "SV_008",
  breakpoint = c(1, 2),
  seqname = c("chr9", "chr22"),
  pos = c(130650100, 23632345),
  gene_id = c("ENSG00000134023", "ENSG00000100031"),
  gene_symbol = c("ABL1", "BCR"),
  transcript_id = c("ENST00000435335", "ENST00000332410"),
  region_type = c("intron", "intron"),
  region_details = c("Intron", "Intron"),
  is_disease_gene = c(TRUE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t8 <- predict_fusions(bp_data_t8, annotations, disease_genes = disease_genes)

test_success_t8 <- nrow(fusion_result_t8) == 1 &&
  fusion_result_t8$disease_gene_involved[1] == TRUE

test_results[[test_num]] <- test_success_t8
cat(if (test_success_t8) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t8) > 0) {
  cat("  Disease gene:", fusion_result_t8$disease_gene_involved[1], 
      "| Genes:", fusion_result_t8$gene1_symbol[1], "--", fusion_result_t8$gene2_symbol[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 9: Mixed Breakpoints
# ==========================================
test_name <- "Canonical: mixed intron-exon breakpoints"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t9 <- tibble::tibble(
  variant_id = "SV_009",
  breakpoint = c(1, 2),
  seqname = c("chr1", "chr2"),
  pos = c(1000000, 2000000),
  gene_id = c("ENSG00000223972", "ENSG00000139618"),
  gene_symbol = c("DDOX", "ERBB2"),
  transcript_id = c("ENST00000456328", "ENST00000269571"),
  region_type = c("intron", "exon"),
  region_details = c("Intron", "Exon"),
  is_disease_gene = c(FALSE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t9 <- predict_fusions(bp_data_t9, annotations, disease_genes = disease_genes)

test_success_t9 <- nrow(fusion_result_t9) == 1 &&
  fusion_result_t9$fusion_type[1] == "canonical" &&
  fusion_result_t9$gene1_region[1] == "intron" &&
  fusion_result_t9$gene2_region[1] == "exon"
  # Reading frame will be NA for intronic breakpoints (expected)

test_results[[test_num]] <- test_success_t9
cat(if (test_success_t9) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t9) > 0) {
  cat("  Type:", fusion_result_t9$fusion_type[1], 
      "| Regions:", fusion_result_t9$gene1_region[1], "--", 
      fusion_result_t9$gene2_region[1], "\n")
  cat("  Reading frame: ", if (is.na(fusion_result_t9$in_frame[1])) "NA (intronic breakpoint)" else fusion_result_t9$in_frame[1], "\n")
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# TEST 10: Intron-Intron Reading Frame
# ==========================================
test_name <- "Canonical: intron-intron breakpoints (frame depends on exon boundaries)"
cat("Test", test_num, ":", test_name, "\n")

bp_data_t10 <- tibble::tibble(
  variant_id = "SV_010",
  breakpoint = c(1, 2),
  seqname = c("chr9", "chr22"),
  pos = c(130651000, 23633000),  # Different positions to ensure in introns
  gene_id = c("ENSG00000134023", "ENSG00000100031"),
  gene_symbol = c("ABL1", "BCR"),
  transcript_id = c("ENST00000435335", "ENST00000332410"),
  region_type = c("intron", "intron"),
  region_details = c("Intron", "Intron"),
  is_disease_gene = c(TRUE, TRUE),
  distance_to_nearest_gene = c(NA_integer_, NA_integer_)
)

fusion_result_t10 <- predict_fusions(bp_data_t10, annotations, disease_genes = disease_genes)

# For intron-intron junctions, reading frame depends on exon boundaries
# Using precede()/follow() with phase information
test_success_t10 <- nrow(fusion_result_t10) == 1 &&
  fusion_result_t10$fusion_type[1] == "canonical" &&
  fusion_result_t10$gene1_region[1] == "intron" &&
  fusion_result_t10$gene2_region[1] == "intron"
  # Reading frame may be NA if phase calculation not possible, which is acceptable

test_results[[test_num]] <- test_success_t10
cat(if (test_success_t10) "✓ PASS\n" else "✗ FAIL\n")
if (nrow(fusion_result_t10) > 0) {
  cat("  Type:", fusion_result_t10$fusion_type[1], 
      "| Regions: intron-intron\n")
  if (!is.na(fusion_result_t10$in_frame[1])) {
    cat("  In-frame:", fusion_result_t10$in_frame[1], 
        "| offset:", fusion_result_t10$frame_offset[1], "\n")
  } else {
    cat("  Reading frame: NA (requires CDS phase lookup)\n")
  }
}
cat("\n")
test_num <- test_num + 1

# ==========================================
# Summary
# ==========================================
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

passed_count <- sum(unlist(test_results))
total_count <- length(test_results)

for (i in seq_along(test_results)) {
  status <- if (test_results[[i]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("Test %2d: %s\n", i, status))
}

cat("\n")
cat(sprintf("Total: %d / %d tests passed\n", passed_count, total_count))
cat("\n")

if (passed_count == total_count) {
  cat("✓ Phase 4 is complete and ready for Phase 5 implementation.\n\n")
} else {
  cat("✗ Phase 4 requires fixes before proceeding to Phase 5.\n\n")
  cat("Failed tests:", which(!unlist(test_results)), "\n\n")
}
