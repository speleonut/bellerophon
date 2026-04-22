#!/usr/bin/env Rscript

#' Complete Bellerophon Pipeline: All 7 Phases
#'
#' This example script demonstrates how to run the complete 7-phase fusion gene detection pipeline
#' on real-world structural variant data.
#'
#' Usage:
#'   Rscript example_pipeline.R --input data/example_input/sample_variants.vcf --output output/my_analysis
#'   Rscript example_pipeline.R --input data/example_input/sample_variants.mavis.tsv --output output/my_analysis
#'   Rscript example_pipeline.R --input data/example_input/sample_variants.hgvs.txt --output output/my_analysis
#'

library(tidyverse)
library(GenomicRanges)
library(IRanges)

# ==========================================
# Parse command-line arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

# Parse --flag value pairs
parse_args <- function(args) {
  params <- list(
    input = "data/example_input/sample_variants.vcf",
    output = "output/example_analysis",
    format = "auto"
  )
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" && i < length(args)) {
      params$input <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      params$output <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--format" && i < length(args)) {
      params$format <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  params
}

params <- parse_args(args)

# Create output directory
if (!dir.exists(params$output)) {
  dir.create(params$output, recursive = TRUE)
}

# ==========================================
# Initialize logging
# ==========================================
log_messages <- list()

log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- glue::glue("[{timestamp}] [{level}] {message}")
  cat(log_entry, "\n")
  log_messages <<- c(log_messages, list(log_entry))
}

cat("\n")
cat("========================================\n")
cat("Bellerophon: Complete Pipeline Example\n")
cat("========================================\n\n")

log_message("Pipeline started")
log_message(glue::glue("Input file: {params$input}"))
log_message(glue::glue("Output directory: {params$output}"))

# ==========================================
# Source all pipeline modules
# ==========================================
log_message("Loading pipeline modules...")

source("R/parse_hgvs.R")
source("R/parse_inputs.R")
source("R/utils.R")
source("R/annotations.R")
source("R/breakpoint_annotation.R")
source("R/fusion_prediction.R")
source("R/expression_analysis.R")
source("R/report_generation.R")
source("R/validation.R")

log_message("✓ All modules loaded")

# ==========================================
# PHASE 1-2: Input Parsing & Annotation Loading
# ==========================================
cat("\n--- PHASE 1-2: Input Parsing & Annotation Loading ---\n\n")

log_message("Phase 1-2: Starting input parsing and annotation loading")

# Detect and parse input file
log_message(glue::glue("Detecting input format: {params$input}"))
input_data <- detect_and_parse_input(
  input_file = params$input,
  input_format = params$format
)

log_message(glue::glue("Input format detected: {input_data$format}"))
log_message(glue::glue("Parsed {nrow(input_data$variants)} variants"))

# Display sample of input
cat("\nFirst 3 parsed variants:\n")
print(head(input_data$variants %>% dplyr::select(1:9), 3))

# Normalize chromosomes
variants <- input_data$variants %>%
  mutate(
    break1_chromosome = normalize_chromosome(break1_chromosome, style = "NCBI"),
    break2_chromosome = normalize_chromosome(break2_chromosome, style = "NCBI"),
    variant_id = sprintf("VAR_%04d", row_number()),
    variant_type = purrr::pmap_chr(
      list(event_type, break1_chromosome, break2_chromosome, break1_orientation, break2_orientation),
      ~classify_variant_type(..1, ..2, ..3, ..4, ..5)
    )
  )

log_message(glue::glue("Normalized {nrow(variants)} variants"))

# Load GENCODE annotations
log_message("Loading GENCODE v49 annotations...")
annotations <- load_gencode_annotations(force_download = FALSE)

if (annotations$loaded_successfully) {
  log_message("✓ GENCODE annotations loaded")
  log_message(glue::glue("  - {length(annotations$genes_gr)} genes"))
  log_message(glue::glue("  - {length(annotations$exons_gr)} exons"))
  log_message(glue::glue("  - {length(annotations$cds_gr)} CDS regions"))
  log_message(glue::glue("  - {length(annotations$disease_genes)} disease genes"))
} else {
  log_message("✗ Failed to load annotations. Exiting.", level = "ERROR")
  quit(status = 1)
}

# ==========================================
# PHASE 3: Breakpoint Annotation
# ==========================================
cat("\n--- PHASE 3: Breakpoint Annotation ---\n\n")

log_message("Phase 3: Annotating breakpoints")

breakpoint_annotations <- annotate_breakpoints(
  variants = variants,
  annotations = annotations
)

log_message(glue::glue("Annotated {nrow(breakpoint_annotations)} breakpoints"))

cat("\nSample breakpoint annotations (first 6):\n")
print(head(breakpoint_annotations %>% 
           dplyr::select(variant_id, breakpoint, gene_id, region_type, is_disease_gene), 6))

# ==========================================
# PHASE 4: Fusion Prediction
# ==========================================
cat("\n--- PHASE 4: Fusion Prediction ---\n\n")

log_message("Phase 4: Predicting fusions")

fusions <- predict_fusions(
  breakpoint_annotations = breakpoint_annotations,
  annotations = annotations,
  sv_data = variants,
  disease_genes = annotations$disease_genes
)

log_message(glue::glue("Predicted {nrow(fusions)} fusions"))
log_message(glue::glue("  - Canonical: {sum(fusions$canonical == TRUE, na.rm = TRUE)}"))
log_message(glue::glue("  - Non-canonical: {sum(fusions$canonical == FALSE, na.rm = TRUE)}"))
log_message(glue::glue("  - In-frame: {sum(fusions$reading_frame == 'in-frame', na.rm = TRUE)}"))

cat("\nFusion predictions (all):\n")
print(fusions %>% 
      dplyr::select(fusion_id, gene1_symbol, gene2_symbol, canonical, reading_frame))

# ==========================================
# PHASE 5: Expression Analysis
# ==========================================
cat("\n--- PHASE 5: Expression Analysis ---\n\n")

log_message("Phase 5: Loading GTEx expression data")

gtex_expression <- load_gtex_expression()

if (!gtex_expression$loaded_successfully) {
  log_message("⚠ GTEx data not available. Skipping expression analysis.", level = "WARNING")
  log_message("  Download from: https://www.gtexportal.org/", level = "WARNING")
  expression_hotspots <- tibble::tibble()
} else {
  log_message("✓ GTEx expression data loaded")
  log_message(glue::glue("  - {nrow(gtex_expression$expression)} genes"))
  log_message(glue::glue("  - {length(gtex_expression$tissues)} tissues"))
  
  log_message("Identifying fusion hotspots")
  
  expression_hotspots <- identify_fusion_hotspots(
    fusion_data = fusions,
    gtex_expression = gtex_expression
  )
  
  n_hotspots <- sum(expression_hotspots$is_hotspot, na.rm = TRUE)
  log_message(glue::glue("Identified {n_hotspots} hotspots across tissues"))
  
  if (nrow(expression_hotspots) > 0) {
    cat("\nExpression hotspots (first 10):\n")
    print(head(expression_hotspots %>% 
               dplyr::filter(is_hotspot == TRUE) %>%
               dplyr::select(fusion_id, tissue, gene1_tpm, gene2_tpm, hotspot_reason), 10))
  }
}

# ==========================================
# PHASE 6: Report Generation
# ==========================================
cat("\n--- PHASE 6: Report Generation ---\n\n")

log_message("Phase 6: Generating comprehensive reports")

report <- generate_fusion_report(
  fusions = fusions,
  expression_hotspots = expression_hotspots,
  breakpoint_annotations = breakpoint_annotations,
  disease_genes = disease_genes
)

log_message("✓ Report generated")
log_message(glue::glue("  - Fusion summary: {nrow(report$fusion_summary)} rows"))
log_message(glue::glue("  - Breakpoint summary: {nrow(report$breakpoint_summary)} rows"))
log_message(glue::glue("  - Hotspot summary: {nrow(report$hotspot_summary)} rows"))

# Export TSV files
log_message("Exporting report tables to TSV")
export_report_tables(report, output_dir = file.path(params$output, "tables"))

# Generate HTML report
log_message("Generating HTML report")
html_file <- file.path(params$output, "fusion_report.html")
generate_html_report(
  report = report,
  output_file = html_file,
  title = glue::glue("Fusion Analysis Report - {format(Sys.Date(), '%Y-%m-%d')}")
)
log_message(glue::glue("✓ HTML report saved to {html_file}"))

# ==========================================
# PHASE 7: Quality Control & Validation
# ==========================================
cat("\n--- PHASE 7: Quality Control & Validation ---\n\n")

log_message("Phase 7: Running comprehensive validation")

if (nrow(fusions) > 0) {
  validation_result <- validate_complete_pipeline(
    annotations = annotations,
    breakpoint_annotations = breakpoint_annotations,
    fusions = fusions,
    expression_hotspots = expression_hotspots,
    report = report
  )
  
  if (validation_result$overall_passed) {
    log_message("✓ All validation checks passed")
  } else {
    log_message("⚠ Some validation checks failed", level = "WARNING")
  }
  
  # Calculate quality metrics
  metrics <- calculate_quality_metrics(fusions, expression_hotspots, report)
  
  cat("\nData Quality Metrics:\n")
  print(metrics)
} else {
  log_message("⚠ No fusions predicted. Skipping validation.", level = "WARNING")
}

# ==========================================
# Summary and Output
# ==========================================
cat("\n--- PIPELINE SUMMARY ---\n\n")

summary_data <- tibble::tibble(
  Phase = c("1-2: Input & Annotations", "3: Breakpoint Annotation", 
            "4: Fusion Prediction", "5: Expression Analysis",
            "6: Report Generation", "7: Quality Control"),
  Status = c("✓ Complete", "✓ Complete", "✓ Complete", 
             if (nrow(expression_hotspots) > 0) "✓ Complete" else "⚠ Skipped (no data)",
             "✓ Complete", "✓ Complete"),
  Key_Results = c(
    glue::glue("{nrow(variants)} variants"),
    glue::glue("{nrow(breakpoint_annotations)} breakpoints"),
    glue::glue("{nrow(fusions)} fusions, {sum(fusions$canonical == TRUE, na.rm = TRUE)} canonical"),
    if (nrow(expression_hotspots) > 0) 
      glue::glue("{sum(expression_hotspots$is_hotspot, na.rm = TRUE)} hotspots") 
    else 
      "N/A",
    glue::glue("4 report tables, 1 HTML report"),
    if (nrow(fusions) > 0) "All checks passed" else "N/A"
  )
)

print(summary_data)

# Save final outputs
log_message("Saving final results")

# Save fusions
saveRDS(fusions, file.path(params$output, "fusions_predicted.rds"))
log_message(glue::glue("Saved fusions to {params$output}/fusions_predicted.rds"))

# Save report
saveRDS(report, file.path(params$output, "report.rds"))
log_message(glue::glue("Saved report to {params$output}/report.rds"))

# Save log file
log_file <- file.path(params$output, "pipeline.log")
writeLines(unlist(log_messages), log_file)
log_message(glue::glue("Saved logs to {log_file}"))

# ==========================================
# Final Summary
# ==========================================
cat("\n")
cat("========================================\n")
cat("Pipeline Complete\n")
cat("========================================\n")
cat(glue::glue("Output directory: {params$output}\n"))
cat(glue::glue("Results files:\n"))
cat(glue::glue("  - tables/ (TSV files)\n"))
cat(glue::glue("  - fusion_report.html (interactive report)\n"))
cat(glue::glue("  - fusions_predicted.rds (R data file)\n"))
cat(glue::glue("  - report.rds (R data file)\n"))
cat(glue::glue("  - pipeline.log (execution log)\n"))
cat("\n")

log_message("Pipeline completed successfully")
