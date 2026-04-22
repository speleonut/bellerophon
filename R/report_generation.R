#' Phase 6: Fusion Report Generation
#'
#' Synthesizes outputs from Phases 1-5 into comprehensive fusion gene reports.
#' Generates summary tables and HTML reports suitable for downstream analysis.
#'
#' @import tidyverse
#' @import GenomicRanges

# ==========================================
# Generate Complete Fusion Report
# ==========================================

#' Generate Comprehensive Fusion Report
#'
#' Orchestrates report generation by combining fusion predictions, 
#' breakpoint annotations, and expression analysis into unified output tables.
#'
#' @param fusions Tibble from predict_fusions() with columns:
#'   fusion_id, variant_type, bp1_gene_id, bp1_gene_symbol, bp2_gene_id, bp2_gene_symbol,
#'   bp1_gene_strand, bp2_gene_strand, is_canonical, canonical_fusion_gene_promoter_breakpoint,
#'   in_frame, frame_offset
#'
#' @param expression_hotspots Tibble from identify_fusion_hotspots() with columns:
#'   fusion_id, tissue, gene1_tpm, gene2_tpm, log2_ratio, is_hotspot, hotspot_reason
#'
#' @param breakpoint_annotations Tibble from annotate_breakpoints() with columns:
#'   variant_id, seqname, pos, breakpoint, variant_type, gene_id, gene_symbol, 
#'   gene_strand, is_exon, num_exon, is_intron, num_intron, is_cds, is_disease_gene
#'
#' @param disease_genes Optional tibble with gene_id and flag columns for disease association
#'
#' @return List containing:
#'   - fusion_summary: Tibble with one row per fusion
#'   - breakpoint_summary: Tibble with detailed breakpoint information
#'   - hotspot_summary: Tibble with tissue-level hotspot data
#'   - report_metadata: Tibble with report generation details
#'
#' @export
generate_fusion_report <- function(fusions, expression_hotspots, breakpoint_annotations, disease_genes = NULL) {
  
  # Validate inputs
  if (!is.data.frame(fusions) || nrow(fusions) == 0) {
    return(list(
      fusion_summary = tibble::tibble(),
      breakpoint_summary = tibble::tibble(),
      hotspot_summary = tibble::tibble(),
      report_metadata = tibble::tibble(error = "No fusions provided")
    ))
  }
  
  # ========== FUSION SUMMARY ==========
  # One row per fusion with key metrics
  # First, extract breakpoint positions from annotations
  breakpoint_positions <- breakpoint_annotations %>%
    dplyr::select(variant_id, breakpoint, seqname, pos) %>%
    tidyr::pivot_wider(
      names_from = breakpoint,
      values_from = c(seqname, pos),
      names_glue = "bp{breakpoint}_{.value}"
    )
  
  fusion_summary <- fusions %>%
    dplyr::mutate(
      # Determine disease gene flags by checking if gene IDs are in disease_genes list
      bp1_is_disease = bp1_gene_id %in% (disease_genes$gene_id %||% character()),
      bp2_is_disease = bp2_gene_id %in% (disease_genes$gene_id %||% character()),
      disease_genes = dplyr::case_when(
        bp1_is_disease & bp2_is_disease ~ paste(bp1_gene_symbol, bp2_gene_symbol, sep = " & "),
        bp1_is_disease ~ bp1_gene_symbol,
        bp2_is_disease ~ bp2_gene_symbol,
        TRUE ~ "None"
      ),
      disease_gene_count = as.integer(bp1_is_disease) + as.integer(bp2_is_disease),
      fusion_type = dplyr::case_when(
        is_canonical == TRUE ~ "Canonical",
        is_canonical == FALSE ~ "Non-canonical",
        TRUE ~ "No fusion"
      ),
      reading_frame = dplyr::case_when(
        in_frame == TRUE ~ "In-frame",
        in_frame == FALSE ~ "Out-of-frame",
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::left_join(
      breakpoint_positions,
      by = c("variant_id" = "variant_id")
    ) %>%
    dplyr::left_join(
      expression_hotspots %>%
        dplyr::filter(is_hotspot) %>%
        dplyr::group_by(fusion_id) %>%
        dplyr::summarise(
          hotspot_tissue_count = dplyr::n(),
          hotspot_tissues = paste(unique(tissue), collapse = "; "),
          .groups = "drop"
        ),
      by = "fusion_id"
    ) %>%
    dplyr::mutate(
      hotspot_tissue_count = ifelse(is.na(hotspot_tissue_count), 0, hotspot_tissue_count),
      hotspot_tissues = ifelse(is.na(hotspot_tissues), "None", hotspot_tissues)
    ) %>%
    dplyr::select(
      fusion_id,
      bp1_gene_symbol,
      bp2_gene_symbol,
      fusion_type,
      reading_frame,
      disease_genes,
      disease_gene_count,
      bp1_seqname,
      bp1_pos,
      bp2_seqname,
      bp2_pos,
      is_canonical,
      in_frame,
      hotspot_tissue_count,
      hotspot_tissues
    )
  
  # ========== BREAKPOINT SUMMARY ==========
  # Detailed breakpoint information with gene-level summary
  breakpoint_summary <- breakpoint_annotations %>%
    dplyr::select(
      variant_id,
      breakpoint,
      seqname,
      pos,
      gene_id,
      gene_symbol,
      gene_strand,
      is_exon,
      num_exon,
      is_intron,
      num_intron,
      is_cds,
      is_disease_gene
    ) %>%
    dplyr::mutate(
      is_disease_gene = ifelse(is.na(is_disease_gene), FALSE, is_disease_gene),
      region_type = dplyr::case_when(
        is_exon ~ "Exon",
        is_intron ~ "Intron",
        is_cds ~ "CDS",
        is.na(gene_id) ~ "Intergenic",
        TRUE ~ "Gene region"
      ),
      gene_disruption_summary = dplyr::if_else(
        is.na(gene_symbol),
        "Intergenic region",
        paste0(gene_symbol, " (", region_type, ")")
      )
    ) %>%
    dplyr::arrange(variant_id, breakpoint)
  
  # ========== HOTSPOT SUMMARY ==========
  # Tissue-level fusion permissiveness analysis
  hotspot_summary <- expression_hotspots %>%
    dplyr::mutate(
      gene1_tpm = round(gene1_tpm, 2),
      gene2_tpm = round(gene2_tpm, 2),
      log2_ratio = round(log2_ratio, 3),
      hotspot_status = dplyr::if_else(is_hotspot, "Hotspot", "Not hotspot")
    ) %>%
    dplyr::select(
      fusion_id,
      tissue,
      gene1_tpm,
      gene2_tpm,
      log2_ratio,
      hotspot_status,
      hotspot_reason
    ) %>%
    dplyr::arrange(fusion_id, tissue)
  
  # ========== REPORT METADATA ==========
  # Summary of report generation
  report_metadata <- tibble::tibble(
    generated_at = Sys.time(),
    total_fusions = nrow(fusions),
    canonical_fusions = sum(fusions$is_canonical == TRUE, na.rm = TRUE),
    non_canonical_fusions = sum(fusions$is_canonical == FALSE, na.rm = TRUE),
    no_fusion_cases = sum(is.na(fusions$is_canonical), na.rm = TRUE),
    predicted_in_frame = sum(fusions$in_frame == TRUE, na.rm = TRUE),
    total_disease_genes = sum(fusion_summary$disease_gene_count > 0, na.rm = TRUE),
    total_tissues_analyzed = n_distinct(expression_hotspots$tissue),
    total_hotspots_identified = sum(expression_hotspots$is_hotspot, na.rm = TRUE)
  )
  
  list(
    fusion_summary = fusion_summary,
    breakpoint_summary = breakpoint_summary,
    hotspot_summary = hotspot_summary,
    report_metadata = report_metadata
  )
}


# ==========================================
# Export Report Tables
# ==========================================

#' Export Report to TSV Files
#'
#' Writes report tables to tab-separated files for easy viewing and downstream use.
#'
#' @param report List returned from generate_fusion_report()
#' @param output_dir Directory to write files to (default: "output/")
#'
#' @return Invisibly returns list of written file paths
#'
#' @export
export_report_tables <- function(report, output_dir = "output/") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  files_written <- list()
  
  # Write fusion summary
  if (nrow(report$fusion_summary) > 0) {
    fusion_file <- file.path(output_dir, paste0("fusion_summary_", timestamp, ".tsv"))
    readr::write_tsv(report$fusion_summary, fusion_file)
    files_written$fusion_summary <- fusion_file
    cat("✓ Fusion summary written to:", fusion_file, "\n")
  }
  
  # Write breakpoint summary
  if (nrow(report$breakpoint_summary) > 0) {
    breakpoint_file <- file.path(output_dir, paste0("breakpoint_summary_", timestamp, ".tsv"))
    readr::write_tsv(report$breakpoint_summary, breakpoint_file)
    files_written$breakpoint_summary <- breakpoint_file
    cat("✓ Breakpoint summary written to:", breakpoint_file, "\n")
  }
  
  # Write hotspot summary
  if (nrow(report$hotspot_summary) > 0) {
    hotspot_file <- file.path(output_dir, paste0("hotspot_summary_", timestamp, ".tsv"))
    readr::write_tsv(report$hotspot_summary, hotspot_file)
    files_written$hotspot_summary <- hotspot_file
    cat("✓ Hotspot summary written to:", hotspot_file, "\n")
  }
  
  # Write metadata
  metadata_file <- file.path(output_dir, paste0("report_metadata_", timestamp, ".tsv"))
  readr::write_tsv(report$report_metadata, metadata_file)
  files_written$metadata <- metadata_file
  cat("✓ Report metadata written to:", metadata_file, "\n")
  
  invisible(files_written)
}


# ==========================================
# Generate HTML Report
# ==========================================

#' Generate HTML Report with Summary Statistics
#'
#' Creates an interactive HTML report with tables, statistics, and visualizations.
#' Useful for sharing results with collaborators.
#'
#' @param report List returned from generate_fusion_report()
#' @param output_file Path to write HTML report (default: "output/fusion_report.html")
#' @param title Report title (default: "Bellerophon Fusion Report")
#'
#' @return Invisibly returns the output file path
#'
#' @export
generate_html_report <- function(report, output_file = "output/fusion_report.html", 
                                  title = "Bellerophon Fusion Report") {
  
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  # Build HTML content
  html_content <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "  <meta charset='utf-8'>",
    paste0("  <title>", title, "</title>"),
    "  <style>",
    "    body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }",
    "    .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }",
    "    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }",
    "    h2 { color: #34495e; margin-top: 30px; }",
    "    .metadata { background-color: #ecf0f1; padding: 15px; border-radius: 4px; margin: 20px 0; }",
    "    .stat-box { display: inline-block; margin: 10px 20px 10px 0; }",
    "    .stat-box strong { color: #3498db; font-size: 18px; }",
    "    table { border-collapse: collapse; width: 100%; margin: 15px 0; }",
    "    th { background-color: #3498db; color: white; padding: 12px; text-align: left; }",
    "    td { padding: 10px; border-bottom: 1px solid #ddd; }",
    "    tr:hover { background-color: #f9f9f9; }",
    "    .hotspot { background-color: #fff9e6; }",
    "    .disease-gene { color: #b80058; font-weight: bold; }",
    "    .canonical { color: #00bbad; font-weight: bold; }",
    "    .non-canonical { color: #ebac23; font-weight: bold; }",
    "    .footer { margin-top: 30px; text-align: center; color: #7f8c8d; font-size: 12px; }",
    "  </style>",
    "</head>",
    "<body>",
    "  <div class='container'>",
    paste0("    <h1>", title, "</h1>"),
    paste0("    <p>Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>")
  )
  
  # Add metadata section
  if (nrow(report$report_metadata) > 0) {
    metadata <- report$report_metadata[1, ]
    html_content <- c(html_content,
      "    <div class='metadata'>",
      "      <h2>Report Summary</h2>",
      paste0("      <div class='stat-box'><strong>", metadata$total_fusions, "</strong> Total Fusions</div>"),
      paste0("      <div class='stat-box'><strong>", metadata$canonical_fusions, "</strong> Canonical</div>"),
      paste0("      <div class='stat-box'><strong>", metadata$non_canonical_fusions, "</strong> Non-Canonical</div>"),
      paste0("      <div class='stat-box'><strong>", metadata$predicted_in_frame, "</strong> In-Frame</div>"),
      paste0("      <div class='stat-box'><strong>", metadata$total_disease_genes, "</strong> Disease Genes</div>"),
      paste0("      <div class='stat-box'><strong>", metadata$total_hotspots_identified, "</strong> Hotspots</div>"),
      "    </div>"
    )
  }
  
  # Add fusion summary table
  if (nrow(report$fusion_summary) > 0) {
    html_content <- c(html_content,
      "    <h2>Fusion Summary</h2>",
      "    <table>",
      paste0("      <tr>",
        "<th>Fusion ID</th>",
        "<th>Gene 1</th>",
        "<th>Gene 2</th>",
        "<th>Type</th>",
        "<th>Reading Frame</th>",
        "<th>Disease Genes</th>",
        "<th>Hotspot Tissues</th>",
        "</tr>")
    )
    
    for (i in seq_len(nrow(report$fusion_summary))) {
      row <- report$fusion_summary[i, ]
      fusion_type_class <- if (grepl("Canonical", row$fusion_type[[1]])) "canonical" else "non-canonical"
      hotspot_class <- if (!is.na(row$hotspot_tissue_count[[1]]) && row$hotspot_tissue_count[[1]] > 0) "hotspot" else ""
      disease_class <- if (!is.na(row$disease_gene_count[[1]]) && row$disease_gene_count[[1]] > 0) "disease-gene" else ""
      
      html_content <- c(html_content,
        paste0("      <tr class='", hotspot_class, "'>",
          "<td>", row$fusion_id[[1]], "</td>",
          "<td>", row$gene1_symbol[[1]], "</td>",
          "<td>", row$gene2_symbol[[1]], "</td>",
          "<td class='", fusion_type_class, "'>", row$fusion_type[[1]], "</td>",
          "<td>", row$reading_frame[[1]] %||% "NA", "</td>",
          "<td class='", disease_class, "'>", row$disease_genes[[1]], "</td>",
          "<td>", row$hotspot_tissue_count[[1]], "</td>",
          "</tr>")
      )
    }
    
    html_content <- c(html_content, "    </table>")
  }
  
  # Add hotspot table (preview)
  if (nrow(report$hotspot_summary) > 0) {
    html_content <- c(html_content,
      "    <h2>Expression Hotspots (Preview)</h2>",
      "    <table>",
      paste0("      <tr>",
        "<th>Fusion ID</th>",
        "<th>Tissue</th>",
        "<th>Gene1 TPM</th>",
        "<th>Gene2 TPM</th>",
        "<th>Log2 Ratio</th>",
        "<th>Status</th>",
        "</tr>")
    )
    
    # Show only first 20 hotspots
    hotspots_preview <- report$hotspot_summary %>%
      dplyr::filter(hotspot_status == "Hotspot") %>%
      dplyr::slice_head(n = 20)
    
    for (i in seq_len(nrow(hotspots_preview))) {
      row <- hotspots_preview[i, ]
      html_content <- c(html_content,
        paste0("      <tr class='hotspot'>",
          "<td>", row$fusion_id[[1]], "</td>",
          "<td>", row$tissue[[1]], "</td>",
          "<td>", row$gene1_tpm[[1]], "</td>",
          "<td>", row$gene2_tpm[[1]], "</td>",
          "<td>", row$log2_ratio[[1]], "</td>",
          "<td>", row$hotspot_status[[1]], "</td>",
          "</tr>")
      )
    }
    
    if (nrow(hotspots_preview) > 0 && nrow(hotspots_preview) < nrow(report$hotspot_summary %>% dplyr::filter(hotspot_status == "Hotspot"))) {
      html_content <- c(html_content,
        paste0("      <tr><td colspan='6' style='text-align:center;'>... and ",
          nrow(report$hotspot_summary %>% dplyr::filter(hotspot_status == "Hotspot")) - nrow(hotspots_preview),
          " more hotspots (see full report tables)</td></tr>")
      )
    }
    
    html_content <- c(html_content, "    </table>")
  }
  
  # Close HTML
  html_content <- c(html_content,
    "    <div class='footer'>",
    "      <p>Bellerophon Fusion Gene Detection Pipeline</p>",
    "      <p>Phase 6: Report Generation</p>",
    "    </div>",
    "  </div>",
    "</body>",
    "</html>"
  )
  
  # Write HTML file
  writeLines(html_content, con = output_file)
  
  cat("✓ HTML report written to:", output_file, "\n")
  
  invisible(output_file)
}


# ==========================================
# Utilities
# ==========================================

#' Print Report Summary to Console
#'
#' Displays a formatted summary of the fusion report to the R console.
#'
#' @param report List returned from generate_fusion_report()
#'
#' @return Invisibly returns NULL
#'
#' @export
print_report_summary <- function(report) {
  
  metadata <- report$report_metadata[1, ]
  
  cat("\n")
  cat("=====================================\n")
  cat("FUSION REPORT SUMMARY\n")
  cat("=====================================\n\n")
  
  cat("Report Metadata:\n")
  cat("  Generated:", format(metadata$generated_at, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("  Total Fusions:", metadata$total_fusions, "\n")
  cat("  Canonical:", metadata$canonical_fusions, "\n")
  cat("  Non-Canonical:", metadata$non_canonical_fusions, "\n")
  cat("  No Fusion Cases:", metadata$no_fusion_cases, "\n\n")
  
  cat("Predicted Outcomes:\n")
  cat("  In-Frame:", metadata$predicted_in_frame, "\n")
  cat("  Disease Genes Involved:", metadata$total_disease_genes, "\n")
  cat("  Tissues Analyzed:", metadata$total_tissues_analyzed, "\n")
  cat("  Hotspots Identified:", metadata$total_hotspots_identified, "\n\n")
  
  cat("Fusion Summary (first 5):\n")
  if (nrow(report$fusion_summary) > 0) {
    print(report$fusion_summary %>% dplyr::slice_head(n = 5), n = Inf)
  } else {
    cat("  No fusions in report\n")
  }
  
  cat("\n")
  
  invisible(NULL)
}
