#' Fusion Gene Prediction - Phase 4 (Revised)
#'
#' This module predicts canonical and non-canonical gene fusions from annotated breakpoints
#' using a decision matrix based on structural variant type and strand orientation.
#'
#' **Canonical Fusion**: Two genes disrupted at different breakpoints with strand orientations
#'   that allow functional fusion genes to form, as determined by the decision matrix.
#'
#' **Non-Canonical Fusion**: Same gene disrupted twice, or one gene disrupted with other
#'   breakpoint intergenic, or strand orientations that don't support fusion.
#'

#' Canonical Fusion Decision Matrix
#'
#' Defines which combinations of variant_type and strand orientations produce canonical fusions
#'
#' @return Tibble with columns: variant_type, strand_bp1, strand_bp2, is_canonical
#'
create_canonical_decision_matrix <- function() {
  
  tibble::tibble(
    variant_type = c(
      # Inversions
      "inversion", "inversion", "inversion", "inversion",
      # Duplications
      "duplication", "duplication", "duplication", "duplication",
      # Deletions
      "deletion", "deletion", "deletion", "deletion",
      # Translocations
      "translocation", "translocation", "translocation", "translocation",
      # Translocations (inverted)
      "translocation_inverted", "translocation_inverted", "translocation_inverted", "translocation_inverted"
    ),
    strand_bp1 = c(
      # Inversions
      "+", "+", "-", "-",
      # Duplications
      "+", "+", "-", "-",
      # Deletions
      "+", "+", "-", "-",
      # Translocations
      "+", "+", "-", "-",
      # Translocations (inverted)
      "+", "+", "-", "-"
    ),
    strand_bp2 = c(
      # Inversions
      "+", "-", "+", "-",
      # Duplications
      "+", "-", "+", "-",
      # Deletions
      "+", "-", "+", "-",
      # Translocations
      "+", "-", "+", "-",
      # Translocations (inverted)
      "-", "+", "+", "-"
    ),
    is_canonical = c(
      # Inversions
      FALSE, TRUE, TRUE, FALSE,
      # Duplications
      TRUE, FALSE, FALSE, TRUE,
      # Deletions
      TRUE, FALSE, FALSE, TRUE,
      # Translocations
      TRUE, FALSE, FALSE, TRUE,
      # Translocations (inverted)
      FALSE, TRUE, TRUE, FALSE
    ),
    promoter_breakpoint = c(
      # Inversions
      NA_integer_, 1L, 1L, NA_integer_,
      # Duplications
      2L, NA_integer_, NA_integer_, 1L,
      # Deletions
      1L, NA_integer_, NA_integer_, 2L,
      # Translocations
      1L, NA_integer_, NA_integer_, 2L,
      # Translocations (inverted)
      NA_integer_, 1L, 1L, NA_integer_
    )
  )
}

#' Predict Fusion Genes from Annotated Breakpoints (Revised)
#'
#' Takes gene-level breakpoint annotations and generates fusion predictions.
#' Creates permutations of breakpoint combinations for each variant and classifies
#' canonicity using the decision matrix.
#'
#' @param breakpoint_annotations Tibble from annotate_breakpoints() with columns:
#'   variant_id, seqname, pos, breakpoint (1 or 2), variant_type, is_gene, gene_id,
#'   gene_symbol, gene_strand, is_exon, num_exon, is_intron, num_intron, is_cds, etc.
#' @param annotations List from load_gencode_annotations()
#' @param disease_genes Tibble from load_disease_genes() (optional)
#'
#' @return Tibble with fusion predictions including fusion_id, gene pair info, and canonicity
#'
predict_fusions <- function(breakpoint_annotations, annotations, disease_genes = NULL) {

  # Initialize results tibble
  results <- tibble::tibble(
    fusion_id = character(),
    variant_id = character(),
    variant_type = character(),
    bp1_gene_id = character(),
    bp1_gene_symbol = character(),
    bp2_gene_id = character(),
    bp2_gene_symbol = character(),
    bp1_gene_strand = character(),
    bp2_gene_strand = character(),
    is_canonical = logical(),
    canonical_fusion_gene_promoter_breakpoint = integer(),
    in_frame = logical(),
    upstream_end_phase = integer(),
    downstream_end_phase = integer(),
    frame_offset = integer()
  )

  # Input validation
  if (nrow(breakpoint_annotations) == 0) {
    return(results)
  }

  # Get unique variants
  unique_variants <- unique(breakpoint_annotations$variant_id)
  
  # Load decision matrix
  decision_matrix <- create_canonical_decision_matrix()

  fusion_results <- list()
  fusion_idx <- 1

  # Process each variant
  for (variant_id in unique_variants) {
    
    # Get annotations for this variant
    variant_annot <- breakpoint_annotations %>%
      dplyr::filter(variant_id == !!variant_id)
    
    # Get breakpoint 1 and 2 annotations
    bp1_annot <- variant_annot %>% dplyr::filter(breakpoint == 1)
    bp2_annot <- variant_annot %>% dplyr::filter(breakpoint == 2)
    
    if (nrow(bp1_annot) == 0 || nrow(bp2_annot) == 0) {
      next  # Skip if either breakpoint missing
    }
    
    variant_type <- bp1_annot$variant_type[1]
    
    # Determine if this variant type supports reciprocal fusions
    # Inversions and translocation_inverted can form reciprocal fusions (bp2::bp1 in addition to bp1::bp2)
    supports_reciprocal <- variant_type %in% c("inversion", "translocation", "translocation_inverted")
    
    # Generate all permutations of bp1 and bp2 genes
    for (i in seq_len(nrow(bp1_annot))) {
      for (j in seq_len(nrow(bp2_annot))) {
        
        # Process both permutations (bp1::bp2 and bp2::bp1) if supports_reciprocal
        # Otherwise only process bp1::bp2
        permutations_to_process <- list(list(bp1 = bp1_annot[i, ], bp2 = bp2_annot[j, ]))
        
        if (supports_reciprocal) {
          # Add reciprocal permutation: bp2::bp1
          permutations_to_process[[2]] <- list(bp1 = bp2_annot[j, ], bp2 = bp1_annot[i, ])
        }
        
        for (perm_idx in seq_along(permutations_to_process)) {
          perm <- permutations_to_process[[perm_idx]]
          bp1_row <- perm$bp1
          bp2_row <- perm$bp2
          
          bp1_gene_id <- bp1_row$gene_id[1]
          bp1_gene_symbol <- bp1_row$gene_symbol[1]
          bp1_gene_strand <- bp1_row$gene_strand[1]
          bp1_is_exon <- bp1_row$is_exon[1]
          bp1_num_exon <- bp1_row$num_exon[1]
          bp1_is_intron <- bp1_row$is_intron[1]
          bp1_num_intron <- bp1_row$num_intron[1]
          
          bp2_gene_id <- bp2_row$gene_id[1]
          bp2_gene_symbol <- bp2_row$gene_symbol[1]
          bp2_gene_strand <- bp2_row$gene_strand[1]
          bp2_is_exon <- bp2_row$is_exon[1]
          bp2_num_exon <- bp2_row$num_exon[1]
          bp2_is_intron <- bp2_row$is_intron[1]
          bp2_num_intron <- bp2_row$num_intron[1]
          
          # Create fusion ID
          bp1_label <- ifelse(bp1_is_exon, "e", ifelse(bp1_is_intron, "i", "0"))
          bp2_label <- ifelse(bp2_is_exon, "e", ifelse(bp2_is_intron, "i", "0"))
          
          bp1_symbol <- if (!is.na(bp1_gene_symbol)) bp1_gene_symbol else bp1_gene_id
          bp2_symbol <- if (!is.na(bp2_gene_symbol)) bp2_gene_symbol else bp2_gene_id
          
          bp1_symbol <- if (is.na(bp1_symbol)) "" else as.character(bp1_symbol)
          bp2_symbol <- if (is.na(bp2_symbol)) "" else as.character(bp2_symbol)
          
          fusion_id <- glue::glue("{bp1_symbol}::{bp2_symbol}_{bp1_label}_{bp2_label}_{variant_id}")
          
          # Determine canonicity - only if both breakpoints have genes
          is_canonical <- FALSE
          promoter_bp <- NA_integer_
          
          if (!is.na(bp1_gene_strand) && !is.na(bp2_gene_strand)) {
            # Look up in decision matrix
            matrix_match <- decision_matrix %>%
              dplyr::filter(
                variant_type == !!variant_type,
                strand_bp1 == !!bp1_gene_strand,
                strand_bp2 == !!bp2_gene_strand
              ) %>%
              dplyr::slice(1)
            
            if (nrow(matrix_match) > 0) {
              is_canonical <- matrix_match$is_canonical[1]
              promoter_bp <- matrix_match$promoter_breakpoint[1]
            }
          }
          
          # For canonical fusions, calculate reading frame
          in_frame_result <- NA
          if (is_canonical && !is.na(promoter_bp)) {
            # Determine which breakpoint is promoter and which is 3' end
            if (promoter_bp == 1) {
              # bp1 has promoter, bp2 provides 3' end
              bp_promoter <- bp1_row
              bp_3end <- bp2_row
            } else {
              # bp2 has promoter, bp1 provides 3' end
              bp_promoter <- bp2_row
              bp_3end <- bp1_row
            }
            
            # Calculate reading frame
            in_frame_result <- predict_reading_frame(
              bp_promoter_annot = bp_promoter,
              bp_3end_annot = bp_3end,
              annotations = annotations,
              variant_type = variant_type
            )
          }
          
          # Extract reading frame results if available
          in_frame_val <- NA
          upstream_end_phase_val <- NA_integer_
          downstream_end_phase_val <- NA_integer_
          frame_offset_val <- NA_integer_
          if (is.list(in_frame_result)) {
            in_frame_val <- in_frame_result$in_frame
            upstream_end_phase_val <- in_frame_result$upstream_end_phase
            downstream_end_phase_val <- in_frame_result$downstream_end_phase
            frame_offset_val <- in_frame_result$offset
          }
          
          fusion_results[[fusion_idx]] <- tibble::tibble(
            fusion_id = fusion_id,
            variant_id = variant_id,
            variant_type = variant_type,
            bp1_gene_id = bp1_gene_id,
            bp1_gene_symbol = bp1_gene_symbol,
            bp2_gene_id = bp2_gene_id,
            bp2_gene_symbol = bp2_gene_symbol,
            bp1_gene_strand = bp1_gene_strand,
            bp2_gene_strand = bp2_gene_strand,
            is_canonical = is_canonical,
            canonical_fusion_gene_promoter_breakpoint = promoter_bp,
            in_frame = in_frame_val,
            upstream_end_phase = upstream_end_phase_val,
            downstream_end_phase = downstream_end_phase_val,
            frame_offset = frame_offset_val
          )
          
          fusion_idx <- fusion_idx + 1
        }  # Close perm_idx loop
      }
    }
  }

  if (length(fusion_results) == 0) {
    return(results)
  }

  result_df <- dplyr::bind_rows(fusion_results)
  return(result_df)
}

#' Predict Reading Frame for Canonical Fusions (Revised)
#'
#' Calculates whether fusion maintains reading frame using CDS phase information.
#' Uses MANE_Select transcripts when available for more robust predictions.
#'
#' **CDS Phase Information**:
#' - Phase 0: Full codon at start of exon
#' - Phase 1: Last 2 bases of codon at start
#' - Phase 2: Last 1 base of codon at start
#' - end_phase = (phase + width) mod 3
#'
#' **In-Frame Conditions**:
#' - For promoter breakpoint: looking upstream (use follow () to find the CDS prior to the breakpoint)
#' - For 3' end breakpoint: looking downstream (use precede () to find CDS after the breakpoint)
#' - downstream_end_phase == upstream_phase (adjusted for codon phase)
#'
#' @param bp_promoter_annot Tibble row for breakpoint providing promoter
#' @param bp_3end_annot Tibble row for breakpoint providing 3' end
#' @param annotations List from load_gencode_annotations() with cds, tx_lookup, etc.
#' @param variant_type Character indicating structural variant type
#'
#' @return List with:
#'   - in_frame: logical (TRUE if reading frame maintained)
#'   - offset: integer (0, 1, or 2 representing codon phase)
#'   Returns NA if frame cannot be determined
#'
predict_reading_frame <- function(bp_promoter_annot, bp_3end_annot, annotations, variant_type = NA_character_) {
  
  # Require gene info for both breakpoints
  if (is.na(bp_promoter_annot$gene_id[1]) || is.na(bp_3end_annot$gene_id[1])) {
    return(NA)
  }
  
  # Require CDS presence at both breakpoints
  if (!bp_promoter_annot$is_cds[1] || !bp_3end_annot$is_cds[1]) {
    return(NA)
  }
  
  # Check if CDS data available with phase information
  if (is.null(annotations$cds) || length(annotations$cds) == 0) {
    return(NA)
  }
  
  cds_gr <- annotations$cds
  if (!("phase" %in% colnames(mcols(cds_gr)))) {
    return(NA)
  }
  
  # Get gene IDs for CDS selection
  promoter_gene_id <- bp_promoter_annot$gene_id[1]
  three_end_gene_id <- bp_3end_annot$gene_id[1]
  
  # Get all CDS for each gene
  cds_promoter_gene <- cds_gr[names(cds_gr) %in% (annotations$tx_lookup %>%
    dplyr::filter(gene_id == !!promoter_gene_id) %>%
    dplyr::pull(tx_id))]
  
  cds_3end_gene <- cds_gr[names(cds_gr) %in% (annotations$tx_lookup %>%
    dplyr::filter(gene_id == !!three_end_gene_id) %>%
    dplyr::pull(tx_id))]
  
  if (length(cds_promoter_gene) == 0 || length(cds_3end_gene) == 0) {
    return(NA)
  }
  
  # Create breakpoint ranges with strand information
  bp_promoter_strand <- bp_promoter_annot$gene_strand[1]
  bp_3end_strand <- bp_3end_annot$gene_strand[1]
  
  bp_promoter_range <- GenomicRanges::GRanges(
    seqnames = bp_promoter_annot$seqname[1],
    ranges = IRanges::IRanges(start = bp_promoter_annot$pos[1], end = bp_promoter_annot$pos[1]),
    strand = ifelse(is.na(bp_promoter_strand), "*", bp_promoter_strand)
  )
  
  bp_3end_range <- GenomicRanges::GRanges(
    seqnames = bp_3end_annot$seqname[1],
    ranges = IRanges::IRanges(start = bp_3end_annot$pos[1], end = bp_3end_annot$pos[1]),
    strand = ifelse(is.na(bp_3end_strand), "*", bp_3end_strand)
  )
  
  # Find upstream CDS from promoter breakpoint (use follow)
  upstream_cds_idx <- GenomicRanges::follow(bp_promoter_range, cds_promoter_gene, select = "all")
  
  if (length(upstream_cds_idx) == 0 || all(is.na(upstream_cds_idx))) {
    return(NA)
  }
  
  # Select CDS - prefer MANE_Select if available
  upstream_cds_candidates <- cds_promoter_gene[subjectHits(upstream_cds_idx)]
  upstream_cds <- select_mane_cds(upstream_cds_candidates, annotations)
  
  # Find downstream CDS from 3' end breakpoint (use precede)
  downstream_cds_idx <- GenomicRanges::precede(bp_3end_range, cds_3end_gene, select = "all")
  
  if (length(downstream_cds_idx) == 0 || all(is.na(downstream_cds_idx))) {
    return(NA)
  }
  
 # Select CDS - prefer MANE_Select if available
  downstream_cds_candidates <- cds_3end_gene[subjectHits(downstream_cds_idx)]    
  downstream_cds <- select_mane_cds(downstream_cds_candidates, annotations)
  
  # Extract phase information
  upstream_phase <- mcols(upstream_cds)$phase[1]
  downstream_phase <- mcols(downstream_cds)$phase[1]
  
  if (is.na(upstream_phase) || is.na(downstream_phase)) {
    return(NA)
  }
  
  # Calculate end phase of upstream CDS
  upstream_end_phase <- (upstream_phase + BiocGenerics::width(upstream_cds)[1]) %% 3
  
  # Check in-frame conditions based on phase relationship
  in_frame <- FALSE
  frame_offset <- NA_integer_
  
  if (upstream_end_phase == 0 && downstream_phase == 0) {
    in_frame <- TRUE
    frame_offset <- 0L
  } else if (upstream_end_phase == 1 && downstream_phase == 2) {
    in_frame <- TRUE
    frame_offset <- 0L
  } else if (upstream_end_phase == 2 && downstream_phase == 1) {
    in_frame <- TRUE
    frame_offset <- 0L
  } else {
    # Out of frame
    in_frame <- FALSE
    frame_offset <- downstream_phase
  }
  
  return(list(
    in_frame = in_frame,
    upstream_end_phase = upstream_end_phase,
    downstream_end_phase = downstream_phase,
    offset = frame_offset
  ))
}

#' Select MANE CDS for Reading Frame Analysis
#'
#' Prefers MANE_Select transcripts when available for more robust predictions.
#' Falls back to first available CDS if no MANE_Select.
#'
#' @param cds_candidates GRanges object with CDS regions
#' @param annotations List from load_gencode_annotations()
#'
#' @return Single GRanges object (best CDS choice)
#'
select_mane_cds <- function(cds_candidates, annotations) {
  
  # Handle empty input
  if (length(cds_candidates) == 0) {
    return(NA)
  }
  
  # If only one candidate, return it directly
  if (length(cds_candidates) == 1) {
    return(cds_candidates)
  }
  
  # Check for MANE_Select tag in metadata
  mc <- mcols(cds_candidates)
  if (!is.null(mc) && "tag" %in% colnames(mc)) {
    tags <- mc$tag
    
    if (!is.null(tags) && length(tags) > 0) {
      # tags is likely a CharacterList with multiple values per element
      # Use sapply to check if MANE_Select appears in any tag for each element
      idx_mane <- which(sapply(tags, function(x) {
        any(grepl("MANE_Select", x, fixed = TRUE))
      }))
      
      if (length(idx_mane) > 0) {
        # Return the first MANE_SELECT CDS
        return(cds_candidates[idx_mane[1]])
      }
    }
  }
  
  # Fallback: return first CDS
  return(cds_candidates[1])
}
