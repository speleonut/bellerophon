#' Phase 7: Quality Control & Validation
#'
#' End-to-end pipeline validation, data integrity checks, and QC metrics.
#'
#' @import tidyverse
#' @import GenomicRanges

# ==========================================
# Validate Annotations
# ==========================================

#' Validate Annotation Module Output
#'
#' Checks the structure and content of loaded GENCODE annotations.
#'
#' @param annotations List returned from load_gencode_annotations()
#'
#' @return Tibble with validation results (check_name, passed, details)
#'
#' @export
validate_annotations <- function(annotations) {
  
  checks <- tibble::tibble()
  
  # Check list structure
  required_elements <- c("genes_gr", "tx_lookup", "exons_gr", "introns_gr", "cds_gr", 
                         "disease_genes", "gtex_tissues")
  missing_elements <- setdiff(required_elements, names(annotations))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "List structure complete",
      passed = length(missing_elements) == 0,
      details = if (length(missing_elements) == 0) 
        "All required elements present" 
      else 
        paste("Missing:", paste(missing_elements, collapse = ", "))
    )
  )
  
  # Validate GRanges objects
  grange_checks <- list(
    genes_gr = annotations$genes_gr,
    exons_gr = annotations$exons_gr,
    introns_gr = annotations$introns_gr,
    cds_gr = annotations$cds_gr
  )
  
  for (name in names(grange_checks)) {
    obj <- grange_checks[[name]]
    passed <- is(obj, "GRanges") && length(obj) > 0
    checks <- checks %>% dplyr::bind_rows(
      tibble::tibble(
        check_name = paste0(name, " is valid GRanges"),
        passed = passed,
        details = if (passed) 
          paste0(length(obj), " ranges") 
        else 
          "Not a GRanges or empty"
      )
    )
  }
  
  # Validate tx_lookup structure
  tx_lookup_valid <- is.data.frame(annotations$tx_lookup) &&
    all(c("tx_id", "gene_id", "gene_name") %in% names(annotations$tx_lookup)) &&
    nrow(annotations$tx_lookup) > 0
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "tx_lookup has required columns",
      passed = tx_lookup_valid,
      details = if (tx_lookup_valid) 
        paste0(nrow(annotations$tx_lookup), " transcripts") 
      else 
        "Missing columns or empty"
    )
  )
  
  # Validate CDS has phase metadata
  cds_has_phase <- "phase" %in% names(GenomicRanges::mcols(annotations$cds_gr)) &&
    "end_phase" %in% names(GenomicRanges::mcols(annotations$cds_gr))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "CDS contains phase metadata",
      passed = cds_has_phase,
      details = if (cds_has_phase) "phase and end_phase columns present" else "Missing phase columns"
    )
  )
  
  # Validate disease genes
  disease_genes_valid <- is.character(annotations$disease_genes) && length(annotations$disease_genes) > 0
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Disease gene list loaded",
      passed = disease_genes_valid,
      details = if (disease_genes_valid) 
        paste0(length(annotations$disease_genes), " disease genes") 
      else 
        "Invalid or empty"
    )
  )
  
  # Validate GTEx tissues
  gtex_tissues_valid <- is.character(annotations$gtex_tissues) && length(annotations$gtex_tissues) > 50
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "GTEx tissues metadata loaded",
      passed = gtex_tissues_valid,
      details = if (gtex_tissues_valid) 
        paste0(length(annotations$gtex_tissues), " tissues") 
      else 
        "Invalid or incomplete"
    )
  )
  
  checks
}


# ==========================================
# Validate Breakpoint Annotations
# ==========================================

#' Validate Breakpoint Annotation Output
#'
#' Checks the structure and logical consistency of breakpoint annotations.
#'
#' @param breakpoint_annotations Tibble from annotate_breakpoints()
#'
#' @return Tibble with validation results
#'
#' @export
validate_breakpoint_annotations <- function(breakpoint_annotations) {
  
  checks <- tibble::tibble()
  
  # Check basic structure
  required_cols <- c("variant_id", "breakpoint", "gene_id", "region_type", "transcript_id")
  missing_cols <- setdiff(required_cols, names(breakpoint_annotations))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Required columns present",
      passed = length(missing_cols) == 0,
      details = if (length(missing_cols) == 0) 
        "All columns present" 
      else 
        paste("Missing:", paste(missing_cols, collapse = ", "))
    )
  )
  
  # Check no missing values in critical columns
  no_missing_critical <- all(!is.na(breakpoint_annotations$variant_id)) &&
    all(!is.na(breakpoint_annotations$breakpoint)) &&
    all(!is.na(breakpoint_annotations$gene_id))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "No missing values in critical columns",
      passed = no_missing_critical,
      details = if (no_missing_critical) 
        "All variant_id, breakpoint, gene_id populated" 
      else 
        "Found missing values"
    )
  )
  
  # Check breakpoint values are valid
  valid_breakpoints <- all(breakpoint_annotations$breakpoint %in% c("bp1", "bp2", NA))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Breakpoint values valid",
      passed = valid_breakpoints,
      details = if (valid_breakpoints) 
        "All breakpoints are 'bp1' or 'bp2'" 
      else 
        "Invalid breakpoint values found"
    )
  )
  
  # Check region types are reasonable
  valid_regions <- all(unique(breakpoint_annotations$region_type) %in% c("exon", "intron", "utr", "intergenic"))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Region types valid",
      passed = valid_regions,
      details = paste("Regions:", paste(unique(breakpoint_annotations$region_type), collapse = ", "))
    )
  )
  
  # Check pairing: variants should have both bp1 and bp2
  variant_breakpoint_pairing <- breakpoint_annotations %>%
    dplyr::group_by(variant_id) %>%
    dplyr::summarise(
      has_bp1 = any(breakpoint == "bp1", na.rm = TRUE),
      has_bp2 = any(breakpoint == "bp2", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!(has_bp1 & has_bp2)) %>%
    nrow()
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Variants have both bp1 and bp2",
      passed = variant_breakpoint_pairing == 0,
      details = if (variant_breakpoint_pairing == 0) 
        "All variants properly paired" 
      else 
        paste(variant_breakpoint_pairing, "variants missing breakpoints")
    )
  )
  
  checks
}


# ==========================================
# Validate Fusion Predictions
# ==========================================

#' Validate Fusion Prediction Output
#'
#' Checks the structure and logical consistency of fusion predictions.
#'
#' @param fusions Tibble from predict_fusions()
#'
#' @return Tibble with validation results
#'
#' @export
validate_fusion_predictions <- function(fusions) {
  
  checks <- tibble::tibble()
  
  # Check required columns
  required_cols <- c("fusion_id", "gene1_symbol", "gene2_symbol", "canonical", "reading_frame")
  missing_cols <- setdiff(required_cols, names(fusions))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Required columns present",
      passed = length(missing_cols) == 0,
      details = if (length(missing_cols) == 0) 
        "All columns present" 
      else 
        paste("Missing:", paste(missing_cols, collapse = ", "))
    )
  )
  
  # Check no missing fusion IDs
  no_missing_ids <- all(!is.na(fusions$fusion_id))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "All fusion IDs populated",
      passed = no_missing_ids,
      details = if (no_missing_ids) 
        paste0(nrow(fusions), " fusions") 
      else 
        "Found missing fusion IDs"
    )
  )
  
  # Check canonical column values
  valid_canonical <- all(fusions$canonical %in% c(TRUE, FALSE, NA))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Canonical values valid",
      passed = valid_canonical,
      details = paste(
        "TRUE:", sum(fusions$canonical == TRUE, na.rm = TRUE),
        "FALSE:", sum(fusions$canonical == FALSE, na.rm = TRUE),
        "NA:", sum(is.na(fusions$canonical))
      )
    )
  )
  
  # Check reading frame values
  valid_frames <- all(fusions$reading_frame %in% c("in-frame", "out-frame", NA))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Reading frame values valid",
      passed = valid_frames,
      details = if (valid_frames) 
        "All values are 'in-frame', 'out-frame', or NA" 
      else 
        "Invalid reading frame values"
    )
  )
  
  # Check gene names are populated (unless intergenic)
  gene_names_valid <- all(!is.na(fusions$gene1_symbol)) &&
    (sum(is.na(fusions$gene2_symbol)) == 0 || all(is.na(fusions$gene2_symbol[is.na(fusions$canonical)])))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Gene names appropriately populated",
      passed = gene_names_valid,
      details = if (gene_names_valid) 
        "Gene naming consistent with fusion type" 
      else 
        "Missing or inconsistent gene names"
    )
  )
  
  checks
}


# ==========================================
# Validate Expression Analysis
# ==========================================

#' Validate Expression Analysis Output
#'
#' Checks the structure and values of expression analysis results.
#'
#' @param expression_hotspots Tibble from identify_fusion_hotspots()
#'
#' @return Tibble with validation results
#'
#' @export
validate_expression_analysis <- function(expression_hotspots) {
  
  checks <- tibble::tibble()
  
  # Check required columns
  required_cols <- c("fusion_id", "tissue", "gene1_tpm", "gene2_tpm", "log2_ratio", "is_hotspot")
  missing_cols <- setdiff(required_cols, names(expression_hotspots))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Required columns present",
      passed = length(missing_cols) == 0,
      details = if (length(missing_cols) == 0) 
        "All columns present" 
      else 
        paste("Missing:", paste(missing_cols, collapse = ", "))
    )
  )
  
  # Check numeric columns are valid
  numeric_valid <- all(!is.na(expression_hotspots$gene1_tpm)) &&
    all(!is.na(expression_hotspots$gene2_tpm)) &&
    all(expression_hotspots$gene1_tpm >= 0) &&
    all(expression_hotspots$gene2_tpm >= 0)
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "TPM values valid (non-negative)",
      passed = numeric_valid,
      details = if (numeric_valid) 
        "All TPM values >= 0" 
      else 
        "Found negative or missing TPM values"
    )
  )
  
  # Check is_hotspot values
  valid_hotspot <- all(is.logical(expression_hotspots$is_hotspot))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "is_hotspot values are logical",
      passed = valid_hotspot,
      details = paste(
        "TRUE:", sum(expression_hotspots$is_hotspot, na.rm = TRUE),
        "FALSE:", sum(!expression_hotspots$is_hotspot, na.rm = TRUE)
      )
    )
  )
  
  # Check hotspot reason populated when is_hotspot = TRUE
  hotspot_reasons_populated <- all(
    !is.na(expression_hotspots$hotspot_reason[expression_hotspots$is_hotspot == TRUE])
  )
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Hotspot reasons populated when is_hotspot=TRUE",
      passed = hotspot_reasons_populated,
      details = if (hotspot_reasons_populated) 
        "All hotspots have reasons" 
      else 
        "Some hotspots missing reasons"
    )
  )
  
  checks
}


# ==========================================
# Validate Reports
# ==========================================

#' Validate Report Output
#'
#' Checks the structure and consistency of generated reports.
#'
#' @param report List returned from generate_fusion_report()
#'
#' @return Tibble with validation results
#'
#' @export
validate_report <- function(report) {
  
  checks <- tibble::tibble()
  
  # Check structure
  required_sections <- c("fusion_summary", "breakpoint_summary", "hotspot_summary", "report_metadata")
  missing_sections <- setdiff(required_sections, names(report))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Report sections present",
      passed = length(missing_sections) == 0,
      details = if (length(missing_sections) == 0) 
        "All 4 sections present" 
      else 
        paste("Missing:", paste(missing_sections, collapse = ", "))
    )
  )
  
  # Check each section has data
  sections_nonempty <- 
    nrow(report$fusion_summary) > 0 &&
    nrow(report$breakpoint_summary) > 0 &&
    nrow(report$hotspot_summary) > 0 &&
    nrow(report$report_metadata) == 1
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "All report sections populated",
      passed = sections_nonempty,
      details = paste(
        "Fusions:", nrow(report$fusion_summary),
        "Breakpoints:", nrow(report$breakpoint_summary),
        "Hotspots:", nrow(report$hotspot_summary),
        "Metadata:", nrow(report$report_metadata)
      )
    )
  )
  
  # Check fusion summary columns
  fusion_cols <- c("fusion_id", "gene1_symbol", "gene2_symbol", "fusion_type", 
                   "reading_frame", "hotspot_tissue_count")
  fusion_cols_valid <- all(fusion_cols %in% names(report$fusion_summary))
  
  checks <- checks %>% dplyr::bind_rows(
    tibble::tibble(
      check_name = "Fusion summary has required columns",
      passed = fusion_cols_valid,
      details = if (fusion_cols_valid) 
        "All columns present" 
      else 
        paste("Missing:", paste(setdiff(fusion_cols, names(report$fusion_summary)), collapse = ", "))
    )
  )
  
  # Check metadata calculations are consistent
  if (nrow(report$report_metadata) > 0) {
    metadata <- report$report_metadata[1, ]
    
    # Check if required columns exist
    required_cols <- c("total_fusions", "canonical_fusions", "non_canonical_fusions", "no_fusion_cases")
    has_required_cols <- all(required_cols %in% names(report$report_metadata))
    
    # Initialize as FALSE
    metadata_consistent <- FALSE
    
    # Only check consistency if columns exist and are numeric
    if (has_required_cols) {
      metadata_consistent <- tryCatch({
        total_fusions_val <- as.numeric(metadata$total_fusions[[1]])
        canonical_val <- as.numeric(metadata$canonical_fusions[[1]])
        non_canonical_val <- as.numeric(metadata$non_canonical_fusions[[1]])
        no_fusion_val <- as.numeric(metadata$no_fusion_cases[[1]])
        
        total_fusions_val == nrow(report$fusion_summary) &&
          (canonical_val + non_canonical_val + no_fusion_val) == nrow(report$fusion_summary)
      }, error = function(e) FALSE)
    }
    
    checks <- checks %>% dplyr::bind_rows(
      tibble::tibble(
        check_name = "Metadata calculations consistent",
        passed = isTRUE(metadata_consistent),
        details = if (isTRUE(metadata_consistent)) 
          "Metadata matches actual data" 
        else if (!has_required_cols)
          "Required metadata columns missing"
        else 
          "Metadata inconsistent with data"
      )
    )
  }
  
  checks
}


# ==========================================
# End-to-End Pipeline Validation
# ==========================================

#' Run Complete Pipeline Validation
#'
#' Orchestrates validation across all phases.
#'
#' @param annotations List from load_gencode_annotations()
#' @param breakpoint_annotations Tibble from annotate_breakpoints()
#' @param fusions Tibble from predict_fusions()
#' @param expression_hotspots Tibble from identify_fusion_hotspots()
#' @param report List from generate_fusion_report()
#'
#' @return List containing:
#'   - validation_results: Tibble with all checks
#'   - phase_summaries: Named list with passed/total per phase
#'   - overall_passed: Boolean indicating if all checks passed
#'
#' @export
validate_complete_pipeline <- function(annotations, breakpoint_annotations, fusions, 
                                       expression_hotspots, report) {
  
  cat("\n")
  cat("========================================\n")
  cat("PIPELINE VALIDATION REPORT\n")
  cat("========================================\n\n")
  
  # Run all validations
  annotation_checks <- validate_annotations(annotations)
  breakpoint_checks <- validate_breakpoint_annotations(breakpoint_annotations)
  fusion_checks <- validate_fusion_predictions(fusions)
  expression_checks <- validate_expression_analysis(expression_hotspots)
  report_checks <- validate_report(report)
  
  # Combine all checks
  all_checks <- dplyr::bind_rows(
    annotation_checks %>% dplyr::mutate(phase = "Annotations"),
    breakpoint_checks %>% dplyr::mutate(phase = "Breakpoints"),
    fusion_checks %>% dplyr::mutate(phase = "Fusions"),
    expression_checks %>% dplyr::mutate(phase = "Expression"),
    report_checks %>% dplyr::mutate(phase = "Reports")
  )
  
  # Calculate summaries
  phase_summaries <- all_checks %>%
    dplyr::group_by(phase) %>%
    dplyr::summarise(
      passed = sum(passed),
      total = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pct_passed = round(100 * passed / total, 1),
      status = ifelse(passed == total, "✓ PASS", "✗ FAIL")
    )
  
  # Print phase summaries
  for (i in seq_len(nrow(phase_summaries))) {
    row <- phase_summaries[i, ]
    cat(sprintf("%s: %s (%d/%d checks)\n", 
                row$phase, row$status, row$passed, row$total))
  }
  
  cat("\n")
  
  # Overall result
  overall_passed <- all(phase_summaries$passed == phase_summaries$total)
  total_checks <- nrow(all_checks)
  total_passed <- sum(all_checks$passed)
  
  cat(sprintf("Overall: %d/%d checks passed\n\n", total_passed, total_checks))
  
  if (!overall_passed) {
    cat("Failed checks:\n")
    failed_checks <- all_checks %>% dplyr::filter(!passed)
    for (i in seq_len(nrow(failed_checks))) {
      row <- failed_checks[i, ]
      cat(sprintf("  %s - %s: %s\n", row$phase, row$check_name, row$details))
    }
    cat("\n")
  }
  
  list(
    validation_results = all_checks,
    phase_summaries = phase_summaries,
    overall_passed = overall_passed,
    total_passed = total_passed,
    total_checks = total_checks
  )
}


# ==========================================
# Data Quality Metrics
# ==========================================

#' Calculate Data Quality Metrics
#'
#' Computes comprehensive quality metrics for the pipeline outputs.
#'
#' @param fusions Tibble from predict_fusions()
#' @param expression_hotspots Tibble from identify_fusion_hotspots()
#' @param report List from generate_fusion_report()
#'
#' @return Tibble with quality metrics
#'
#' @export
calculate_quality_metrics <- function(fusions, expression_hotspots, report) {
  
  metrics <- tibble::tibble()
  
  # Fusion discovery metrics
  n_fusions <- nrow(fusions)
  n_canonical <- sum(fusions$canonical == TRUE, na.rm = TRUE)
  pct_canonical <- if (n_fusions > 0) 100 * n_canonical / n_fusions else 0
  
  metrics <- metrics %>% dplyr::bind_rows(
    tibble::tibble(
      category = "Fusion Discovery",
      metric = "Total fusions",
      value = n_fusions
    ),
    tibble::tibble(
      category = "Fusion Discovery",
      metric = "Canonical fusions",
      value = n_canonical
    ),
    tibble::tibble(
      category = "Fusion Discovery",
      metric = "% Canonical",
      value = round(pct_canonical, 1)
    )
  )
  
  # Reading frame metrics
  n_inframe <- sum(fusions$reading_frame == "in-frame", na.rm = TRUE)
  pct_inframe <- if (n_canonical > 0) 100 * n_inframe / n_canonical else 0
  
  metrics <- metrics %>% dplyr::bind_rows(
    tibble::tibble(
      category = "Reading Frame",
      metric = "In-frame predictions",
      value = n_inframe
    ),
    tibble::tibble(
      category = "Reading Frame",
      metric = "% In-frame (of canonical)",
      value = round(pct_inframe, 1)
    )
  )
  
  # Expression hotspot metrics
  n_hotspots <- sum(expression_hotspots$is_hotspot, na.rm = TRUE)
  n_tissues <- n_distinct(expression_hotspots$tissue)
  pct_hotspot <- if (nrow(expression_hotspots) > 0) 100 * n_hotspots / nrow(expression_hotspots) else 0
  
  metrics <- metrics %>% dplyr::bind_rows(
    tibble::tibble(
      category = "Expression Analysis",
      metric = "Total tissue-fusion results",
      value = nrow(expression_hotspots)
    ),
    tibble::tibble(
      category = "Expression Analysis",
      metric = "Hotspots identified",
      value = n_hotspots
    ),
    tibble::tibble(
      category = "Expression Analysis",
      metric = "% Results as hotspots",
      value = round(pct_hotspot, 1)
    ),
    tibble::tibble(
      category = "Expression Analysis",
      metric = "Unique tissues",
      value = n_tissues
    )
  )
  
  # Data completeness
  fusion_completeness <- round(100 * (nrow(fusions) - sum(is.na(fusions$reading_frame))) / nrow(fusions), 1)
  
  metrics <- metrics %>% dplyr::bind_rows(
    tibble::tibble(
      category = "Data Completeness",
      metric = "Fusion reading frames populated (%)",
      value = fusion_completeness
    )
  )
  
  metrics
}


# ==========================================
# Print Validation Summary
# ==========================================

#' Print Validation Summary Report
#'
#' Formats and displays validation results for user review.
#'
#' @param validation_result List returned from validate_complete_pipeline()
#'
#' @export
print_validation_summary <- function(validation_result) {
  
  cat("\n")
  cat("========================================\n")
  cat("VALIDATION SUMMARY\n")
  cat("========================================\n\n")
  
  phase_summaries <- validation_result$phase_summaries
  
  for (i in seq_len(nrow(phase_summaries))) {
    row <- phase_summaries[i, ]
    status_symbol <- if (row$passed == row$total) "✓" else "✗"
    cat(sprintf("%s Phase %s: %d/%d checks passed\n", 
                status_symbol, row$phase, row$passed, row$total))
  }
  
  cat("\n")
  
  if (validation_result$overall_passed) {
    cat("✓ ALL CHECKS PASSED\n")
    cat("Pipeline is ready for production use.\n\n")
  } else {
    cat("✗ SOME CHECKS FAILED\n")
    cat("Review failed checks above before proceeding.\n\n")
  }
  
  invisible(NULL)
}
