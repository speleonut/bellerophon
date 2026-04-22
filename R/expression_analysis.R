#' Expression Analysis - Phase 5
#'
#' This module queries GTEx tissue expression data and identifies tissues
#' permissive for fusion gene expression. Includes visualization functions.
#'
#' Expression hotspots are tissues where fusion expression is likely:
#' 1. Upstream gene expressed (TPM >= 1) AND downstream not (TPM < 1)
#' 2. Large expression difference: log2(upstream_tpm + 0.1) - log2(downstream_tpm + 0.1) > 1
#'

#' Load GTEx Expression Data
#'
#' Loads GTEx v11 median TPM expression from local or remote file.
#' GCT format: skip 2 header lines, first column is gene ID (ENSG with version),
#' second column is gene symbol, remaining columns are tissues.
#'
#' @param gtex_file Path to GTEx gct.gz file
#'
#' @return List containing:
#'   - expression: tibble of expression data (gene_id, gene_symbol, tissue columns with TPM values)
#'   - tissues: character vector of tissue names
#'   - loaded_successfully: logical
#'
load_gtex_expression <- function(gtex_file = "data/annotations/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz") {
  
  cat("Loading GTEx expression data...\n")
  
  tryCatch({
    # Check if local file exists, if not download
    if (!file.exists(gtex_file)) {
      cat("  Downloading GTEx data from remote server...\n")
      gtex_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
      download.file(gtex_url, gtex_file, mode = "wb")
      cat("  ✓ Downloaded GTEx data\n")
    }
    
    # Read GCT file: skip 2 header lines, read tab-separated
    gtex_data <- readr::read_tsv(
      gtex_file,
      skip = 2,  # Skip "#1.2" and "74628  68" lines
      col_types = readr::cols(.default = readr::col_double(), Name = readr::col_character(), Description = readr::col_character())
    )
    
    # Rename first column and extract tissue names
    gtex_data <- gtex_data %>%
      dplyr::rename(gene_id = Name, gene_symbol = Description) %>%
      tibble::as_tibble()
    
    # Extract tissue names (all columns except gene_id and gene_symbol)
    tissues <- names(gtex_data)[!(names(gtex_data) %in% c("gene_id", "gene_symbol"))]
    
    cat("  ✓ Loaded GTEx expression data: ", nrow(gtex_data), " genes × ", 
        length(tissues), " tissues\n", sep = "")
    
    result <- list(
      expression = gtex_data,
      tissues = tissues,
      loaded_successfully = TRUE
    )
    
    return(result)
    
  }, error = function(e) {
    cat("✗ Error loading GTEx expression data: ", as.character(e), "\n", sep = "")
    return(list(loaded_successfully = FALSE, error = as.character(e)))
  })
}

#' Query Gene Expression Across Tissues
#'
#' Returns expression data for specified genes across all tissues.
#' Strips version numbers from gene IDs for matching.
#'
#' @param gene_ids Character vector of ENSG IDs (with or without version numbers)
#' @param gtex_expression List from load_gtex_expression()
#' @param tissues Optional character vector to filter specific tissues
#'
#' @return Tibble with columns:
#'   - gene_id: ENSG ID (version stripped)
#'   - gene_symbol: gene name
#'   - tissue: tissue name
#'   - tpm: median TPM value
#'   - log2_tpm: log2(TPM + 0.1)
#'
query_gene_expression <- function(gene_ids, gtex_expression, tissues = NULL) {
  
  if (!gtex_expression$loaded_successfully) {
    cat("✗ GTEx expression data not loaded\n")
    return(tibble::tibble())
  }
  
  # Strip version numbers from input gene IDs
  gene_ids_stripped <- gsub("\\.[0-9]+$", "", gene_ids)
  
  # Filter GTEx data for requested genes
  expr_data <- gtex_expression$expression %>%
    dplyr::mutate(gene_id_stripped = gsub("\\.[0-9]+$", "", gene_id)) %>%
    dplyr::filter(gene_id_stripped %in% gene_ids_stripped) %>%
    dplyr::select(-gene_id_stripped)
  
  if (nrow(expr_data) == 0) {
    cat("Warning: No genes found in GTEx data\n")
    return(tibble::tibble())
  }
  
  # Select tissues if specified
  tissue_cols <- if (!is.null(tissues)) {
    intersect(tissues, gtex_expression$tissues)
  } else {
    gtex_expression$tissues
  }
  
  # Pivot to long format
  result <- expr_data %>%
    dplyr::select(gene_id, gene_symbol, all_of(tissue_cols)) %>%
    tidyr::pivot_longer(
      cols = -c(gene_id, gene_symbol),
      names_to = "tissue",
      values_to = "tpm"
    ) %>%
    dplyr::mutate(
      tpm = tidyr::replace_na(tpm, 0),
      log2_tpm = log2(tpm + 0.1)
    ) %>%
    tibble::as_tibble()
  
  return(result)
}

#' Identify Fusion Expression Hotspots
#'
#' For each tissue, identifies if expression profile is permissive for fusion.
#' Flags tissues where upstream gene (bp1) is expressed but downstream (bp2) is not,
#' or where expression difference is large (log2 ratio > 1).
#'
#' **Hotspot Conditions**:
#' 1. Upstream TPM >= 1 AND downstream TPM < 1
#' 2. log2(upstream_tpm + 0.1) - log2(downstream_tpm + 0.1) > 1
#'
#' @param fusion_data Tibble with columns: fusion_id, bp1_gene_id, bp1_gene_symbol, bp2_gene_id, bp2_gene_symbol
#' @param gtex_expression List from load_gtex_expression()
#'
#' @return Tibble with columns:
#'   - fusion_id: fusion identifier
#'   - gene1_id, gene1_symbol: upstream gene (bp1)
#'   - gene2_id, gene2_symbol: downstream gene (bp2)
#'   - tissue: tissue name
#'   - gene1_tpm, gene1_log2tpm: upstream expression
#'   - gene2_tpm, gene2_log2tpm: downstream expression
#'   - is_hotspot: logical (TRUE if tissue permits fusion)
#'   - hotspot_reason: reason for hotspot designation
#'
identify_fusion_hotspots <- function(fusion_data, gtex_expression) {
  
  if (!gtex_expression$loaded_successfully || nrow(fusion_data) == 0) {
    return(tibble::tibble())
  }
  
  # Extract unique gene IDs from bp1 and bp2 columns
  gene_ids <- unique(c(
    fusion_data$bp1_gene_id, fusion_data$bp2_gene_id
  )) %>% na.omit()
  
  # Query expression for all genes in fusions
  all_expr <- query_gene_expression(gene_ids, gtex_expression)
  
  if (nrow(all_expr) == 0) {
    return(tibble::tibble())
  }
  
  # Process each fusion
  hotspot_results <- list()
  result_idx <- 1
  
  for (i in seq_len(nrow(fusion_data))) {
    fusion <- fusion_data[i, ]
    fusion_id <- fusion$fusion_id[1]
    gene1_id <- gsub("\\.[0-9]+$", "", fusion$bp1_gene_id[1])
    gene1_symbol <- fusion$bp1_gene_symbol[1]
    gene2_id <- gsub("\\.[0-9]+$", "", fusion$bp2_gene_id[1])
    gene2_symbol <- fusion$bp2_gene_symbol[1]
    
    # Get expression for both genes
    gene1_expr <- all_expr %>%
      dplyr::filter(stringr::str_detect(gene_id, paste0("^", gene1_id, "(\\.\\d+)?$"))) %>%
      dplyr::rename(gene1_tpm = tpm, gene1_log2tpm = log2_tpm) %>%
      dplyr::select(tissue, gene1_tpm, gene1_log2tpm)
    
    gene2_expr <- all_expr %>%
      dplyr::filter(stringr::str_detect(gene_id, paste0("^", gene2_id, "(\\.\\d+)?$"))) %>%
      dplyr::rename(gene2_tpm = tpm, gene2_log2tpm = log2_tpm) %>%
      dplyr::select(tissue, gene2_tpm, gene2_log2tpm)
    
    # Join expression data by tissue
    if (nrow(gene1_expr) > 0 && nrow(gene2_expr) > 0) {
      tissue_expr <- dplyr::full_join(
        gene1_expr, gene2_expr,
        by = "tissue"
      ) %>%
        tidyr::replace_na(list(gene1_tpm = 0, gene1_log2tpm = log2(0.1),
                               gene2_tpm = 0, gene2_log2tpm = log2(0.1)))
      
      # Identify hotspots
      tissue_expr <- tissue_expr %>%
        dplyr::mutate(
          # Condition 1: upstream expressed, downstream not
          hotspot_cond1 = gene1_tpm >= 1 & gene2_tpm < 1,
          # Condition 2: large expression difference
          log2_ratio = gene1_log2tpm - gene2_log2tpm,
          hotspot_cond2 = log2_ratio > 1,
          is_hotspot = hotspot_cond1 | hotspot_cond2,
          hotspot_reason = dplyr::case_when(
            hotspot_cond1 & hotspot_cond2 ~ "Upstream expressed; downstream silent; large ratio",
            hotspot_cond1 ~ "Upstream expressed; downstream silent",
            hotspot_cond2 ~ "Large expression difference (log2 ratio > 1)",
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::select(-hotspot_cond1, -hotspot_cond2) %>%
        dplyr::mutate(
          fusion_id = fusion_id,
          gene1_id = gene1_id,
          gene1_symbol = gene1_symbol,
          gene2_id = gene2_id,
          gene2_symbol = gene2_symbol,
          .before = tissue
        )
      
      hotspot_results[[result_idx]] <- tissue_expr
      result_idx <- result_idx + 1
    }
  }
  
  if (length(hotspot_results) == 0) {
    return(tibble::tibble())
  }
  
  result_df <- dplyr::bind_rows(hotspot_results)
  return(result_df)
}

#' Plot Fusion Gene Expression Heatmap
#'
#' Generates a heatmap showing expression of fusion partner genes across tissues.
#' Useful for visualizing which tissues may be permissive for fusion expression.
#'
#' @param fusion_id Fusion identifier
#' @param gene1_symbol Symbol of upstream gene
#' @param gene2_symbol Symbol of downstream gene
#' @param gtex_expression List from load_gtex_expression()
#' @param output_file Optional file path to save heatmap (PDF or PNG)
#'
#' @return Invisibly returns the heatmap plot object
#'
plot_fusion_expression_heatmap <- function(gene1_symbol, gene2_symbol, 
                                          gtex_expression, 
                                          output_file = NULL,
                                          title = NULL) {
  
  if (!gtex_expression$loaded_successfully) {
    cat("✗ GTEx expression data not loaded\n")
    return(invisible(NULL))
  }
  
  # Query expression for both genes
  expr_data <- query_gene_expression(
    c(gene1_symbol, gene2_symbol),
    gtex_expression
  )
  
  if (nrow(expr_data) == 0) {
    cat("✗ No expression data found for genes\n")
    return(invisible(NULL))
  }
  
  # Prepare matrix for heatmap (log2 TPM)
  heatmap_matrix <- expr_data %>%
    dplyr::select(gene_symbol, tissue, log2_tpm) %>%
    tidyr::pivot_wider(
      names_from = tissue,
      values_from = log2_tpm,
      values_fill = log2(0.1)
    ) %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()
  
  # Create color palette
  colors <- grDevices::colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100)
  
  # Generate heatmap title
  if (is.null(title)) {
    title <- paste(gene1_symbol, "--", gene2_symbol, "Expression Across Tissues")
  }
  
  # Create heatmap
  hm <- pheatmap::pheatmap(
    heatmap_matrix,
    cluster_cols = TRUE,
    cluster_rows = FALSE,
    color = colors,
    main = title,
    breaks = seq(-4, 12, length.out = 101),
    display_numbers = FALSE,
    fontsize = 10,
    cellwidth = 8,
    cellheight = 15
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = 14, height = 6)
    print(hm)
    grDevices::dev.off()
    cat("✓ Heatmap saved to:", output_file, "\n")
  }
  
  return(invisible(hm))
}
