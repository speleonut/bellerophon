#' Breakpoint Annotation - Phase 3
#'
#' This module identifies genes and transcripts disrupted by structural variant breakpoints.
#' Maps genomic coordinates to gene/transcript features and classifies disruption types.
#'

#' Split Multi-Breakpoint Variants into Separate Breakpoint Sets
#'
#' Converts a tibble with break1 and break2 columns into two separate tibbles,
#' one for each breakpoint, with columns renamed for compatibility with annotation.
#'
#' @param variants Tibble with columns: variant_id, break1_chromosome, break1_position,
#'   break1_orientation, break2_chromosome, break2_position, break2_orientation
#'
#' @return List with two elements:
#'   - bp1: Tibble with columns variant_id, seqname, pos, strand, breakpoint=1
#'   - bp2: Tibble with columns variant_id, seqname, pos, strand, breakpoint=2
#'
split_breakpoints_for_annotation <- function(variants) {
  
  # Create breakpoint 1 tibble
  variants_bp1 <- variants %>%
    dplyr::select(
      variant_id,
      seqname = break1_chromosome,
      pos = break1_position,
      strand = break1_orientation
    ) %>%
    dplyr::mutate(breakpoint = 1L)
  
  # Create breakpoint 2 tibble
  variants_bp2 <- variants %>%
    dplyr::select(
      variant_id,
      seqname = break2_chromosome,
      pos = break2_position,
      strand = break2_orientation
    ) %>%
    dplyr::mutate(breakpoint = 2L)
  
  return(list(bp1 = variants_bp1, bp2 = variants_bp2))
}

#' Annotate Structural Variant Breakpoints with Gene Information
#'
#' For each breakpoint in input variants, identifies overlapping genes and transcripts
#' and classifies which genomic region (exon, intron, UTR, etc.) is disrupted.
#'
#' @param variants GRanges or data.frame with breakpoint coordinates.
#'   Can be in two formats:
#'   1. Multi-breakpoint: data.frame with columns break1_chromosome, break1_position, break1_orientation,
#'      break2_chromosome, break2_position, break2_orientation, variant_id, variant_type
#'   2. Single format: data.frame/GRanges with columns seqname/chrom, pos, strand, variant_id
#'
#' @param annotations List returned from load_gencode_annotations()
#'   Must contain: genes, transcripts, exons, introns, cds GRanges/tibbles
#'
#' @param disease_genes Tibble from load_disease_genes() with gene_id column
#'
#' @return Tibble with columns:
#'   - variant_id: identifier for the structural variant
#'   - seqname: chromosome
#'   - pos: breakpoint position
#'   - breakpoint: 1 or 2 (first or second breakpoint)
#'   - variant_type: type of structural variant
#'   - is_gene: TRUE if breakpoint overlaps any gene
#'   - gene_id: ENSEMBL gene ID (NA if intergenic)
#'   - gene_symbol: gene name (NA if intergenic)
#'   - gene_strand: strand of gene (NA if intergenic)
#'   - is_exon: TRUE if breakpoint overlaps exon(s)
#'   - num_exon: count of exons overlapping this breakpoint for this gene
#'   - is_intron: TRUE if breakpoint overlaps intron(s)
#'   - num_intron: count of introns overlapping this breakpoint for this gene
#'   - is_cds: TRUE if breakpoint overlaps CDS
#'   - is_disease_gene: TRUE/FALSE flag
#'
annotate_breakpoints <- function(variants, annotations, disease_genes = NULL) {
  
  # Check for required annotation objects
  if (!("exons" %in% names(annotations)) || is.null(annotations$exons)) {
    warning("Annotations missing exons GRanges object. Cannot annotate breakpoints.")
    return(tibble::tibble(
      variant_id = character(),
      breakpoint = integer(),
      seqname = character(),
      pos = integer(),
      gene_id = character(),
      gene_symbol = character(),
      transcript_id = character(),
      region_type = character(),
      region_details = character(),
      is_disease_gene = logical(),
      distance_to_nearest_gene = integer()
    ))
  }
  
  # Detect if input has multiple breakpoints per variant (break1/break2 format)
  is_multi_breakpoint <- is.data.frame(variants) && 
    ("break1_chromosome" %in% names(variants)) && 
    ("break2_chromosome" %in% names(variants))
  
  if (is_multi_breakpoint) {
    # Split variants into separate breakpoint sets
    breakpoint_sets <- split_breakpoints_for_annotation(variants)
    
    # Annotate each breakpoint set
    annotated_results <- list()
    for (bp_num in c(1, 2)) {
      bp_name <- paste0("bp", bp_num)
      bp_variants <- breakpoint_sets[[bp_name]]
      
      # Annotate this breakpoint set
      bp_annotations <- annotate_breakpoint_set(
        breakpoint_variants = bp_variants,
        annotations = annotations,
        disease_genes = disease_genes,
        sv_variants = variants
      )
      
      annotated_results[[bp_num]] <- bp_annotations
    }
    
    # Combine annotations from both breakpoints
    combined_results <- dplyr::bind_rows(annotated_results)
    return(combined_results)
    
  } else {
    # Single breakpoint format - use original logic
    annotated_results <- annotate_breakpoint_set(
      breakpoint_variants = variants,
      annotations = annotations,
      disease_genes = disease_genes,
      sv_variants = NULL
    )
    return(annotated_results)
  }
}

#' Annotate a Single Breakpoint Set (Revised - Gene-Level Summary)
#'
#' Performs annotation for a set of breakpoints and summarizes results at the gene level.
#' Returns one row per breakpoint-gene pair with counts of overlapping exons/introns/CDS.
#'
#' @param breakpoint_variants Data frame with columns: variant_id, seqname, pos, strand, breakpoint
#' @param annotations Annotation list from load_gencode_annotations()
#' @param disease_genes Disease genes tibble
#' @param sv_variants Original SV variants tibble (for variant_type information)
#'
#' @return Tibble with gene-level summary per breakpoint
#'
annotate_breakpoint_set <- function(breakpoint_variants, annotations, disease_genes = NULL, sv_variants = NULL) {
  
  results <- tibble::tibble(
    variant_id = character(),
    seqname = character(),
    pos = integer(),
    breakpoint = integer(),
    variant_type = character(),
    is_gene = logical(),
    gene_id = character(),
    gene_symbol = character(),
    gene_strand = character(),
    is_exon = logical(),
    num_exon = integer(),
    is_intron = logical(),
    num_intron = integer(),
    is_cds = logical(),
    is_disease_gene = logical()
  )
  
  # Convert variants to GRanges if needed
  if (is.data.frame(breakpoint_variants) && !inherits(breakpoint_variants, "GRanges")) {
    variants_gr <- convert_df_to_granges(breakpoint_variants)
  } else {
    variants_gr <- breakpoint_variants
  }
  
  if (length(variants_gr) == 0) {
    return(results)
  }
  
  # Process each breakpoint
  annotated_results <- list()
  result_idx <- 1
  
  for (i in seq_len(length(variants_gr))) {
    variant <- variants_gr[i]
    variant_id <- names(variant)[1] %||% as.character(i)
    seqname <- as.character(seqnames(variant))[1]
    pos <- start(variant)[1]
    strand_val <- as.character(strand(variant))[1]
    
    # Get breakpoint number from input data
    breakpoint_num <- 1
    variant_type_val <- NA_character_
    if (is.data.frame(breakpoint_variants)) {
      if ("breakpoint" %in% names(breakpoint_variants) && i <= nrow(breakpoint_variants)) {
        breakpoint_num <- breakpoint_variants$breakpoint[i]
      }
      # Get variant_type if available
      if ("variant_type" %in% names(breakpoint_variants) && i <= nrow(breakpoint_variants)) {
        variant_type_val <- breakpoint_variants$variant_type[i]
      }
    }
    
    # Alternatively get variant_type from sv_variants
    if (is.na(variant_type_val) && !is.null(sv_variants)) {
      sv_row <- sv_variants %>%
        dplyr::filter(variant_id == !!variant_id) %>%
        dplyr::slice(1)
      if (nrow(sv_row) > 0) {
        variant_type_val <- sv_row$variant_type[1]
      }
    }
    
    # Query for overlapping exons
    exon_overlaps <- IRanges::findOverlaps(variant, annotations$exons, ignore.strand = TRUE)
    
    # Query for overlapping introns
    intron_overlaps <- IRanges::findOverlaps(variant, annotations$introns, ignore.strand = TRUE)
    
    # Note: CDS overlap is NOT checked here because introns within CDS boundaries won't overlap directly.
    # Instead, predict_reading_frame() will independently retrieve all CDS for candidate genes via tx_lookup.
    
    # Collect all unique genes at this breakpoint
    all_genes <- character()
    
    if (length(exon_overlaps) > 0) {
      exon_gene_ids <- annotations$tx_lookup %>%
        dplyr::filter(tx_id %in% names(annotations$exons)[subjectHits(exon_overlaps)]) %>%
        dplyr::pull(gene_id) %>%
        unique()
      all_genes <- c(all_genes, exon_gene_ids)
    }
    
    if (length(intron_overlaps) > 0) {
      intron_gene_ids <- annotations$tx_lookup %>%
        dplyr::filter(tx_id %in% names(annotations$introns)[subjectHits(intron_overlaps)]) %>%
        dplyr::pull(gene_id) %>%
        unique()
      all_genes <- c(all_genes, intron_gene_ids)
    }
    
    all_genes <- unique(all_genes)
    
    # If genes found, create one row per gene
    if (length(all_genes) > 0) {
      for (gene_id in all_genes) {
        # Count exons for this gene
        exon_count <- 0
        if (length(exon_overlaps) > 0) {
          exon_tx_ids <- names(annotations$exons)[subjectHits(exon_overlaps)]
          exon_genes <- annotations$tx_lookup %>%
            dplyr::filter(tx_id %in% exon_tx_ids) %>%
            dplyr::pull(gene_id)
          exon_count <- sum(exon_genes == gene_id)
        }
        
        # Count introns for this gene
        intron_count <- 0
        if (length(intron_overlaps) > 0) {
          intron_tx_ids <- names(annotations$introns)[subjectHits(intron_overlaps)]
          intron_genes <- annotations$tx_lookup %>%
            dplyr::filter(tx_id %in% intron_tx_ids) %>%
            dplyr::pull(gene_id)
          intron_count <- sum(intron_genes == gene_id)
        }
        
        # Check if gene has CDS (will be determined by predict_reading_frame if needed)
        # Note: We check if the gene exists in any transcripts with CDS data
        has_cds <- FALSE
        if (!is.null(annotations$tx_lookup)) {
          gene_txs <- annotations$tx_lookup %>%
            dplyr::filter(gene_id == !!gene_id) %>%
            dplyr::pull(tx_id)
          # Assume CDS present if we have transcripts for this gene
          # predict_reading_frame will verify actual CDS overlap
          has_cds <- length(gene_txs) > 0
        }
        
        # Get gene info
        gene_info <- get_gene_info(gene_id, annotations, disease_genes)
        
        # Get gene strand
        gene_strand <- NA_character_
        if (!is.null(annotations$gene_lookup)) {
          gene_row <- annotations$gene_lookup %>%
            dplyr::filter(gene_id == !!gene_id) %>%
            dplyr::slice(1)
          if (nrow(gene_row) > 0) {
            gene_strand <- gene_row$strand[1]
          }
        }
        
        annotated_results[[result_idx]] <- tibble::tibble(
          variant_id = variant_id,
          seqname = seqname,
          pos = pos,
          breakpoint = breakpoint_num,
          variant_type = variant_type_val,
          is_gene = TRUE,
          gene_id = gene_id,
          gene_symbol = gene_info$gene_symbol,
          gene_strand = gene_strand,
          is_exon = exon_count > 0,
          num_exon = as.integer(exon_count),
          is_intron = intron_count > 0,
          num_intron = as.integer(intron_count),
          is_cds = has_cds,
          is_disease_gene = gene_info$is_disease_gene
        )
        
        result_idx <- result_idx + 1
      }
    } else {
      # Intergenic breakpoint - create one row with NA gene values
      annotated_results[[result_idx]] <- tibble::tibble(
        variant_id = variant_id,
        seqname = seqname,
        pos = pos,
        breakpoint = breakpoint_num,
        variant_type = variant_type_val,
        is_gene = FALSE,
        gene_id = NA_character_,
        gene_symbol = NA_character_,
        gene_strand = NA_character_,
        is_exon = FALSE,
        num_exon = 0L,
        is_intron = FALSE,
        num_intron = 0L,
        is_cds = FALSE,
        is_disease_gene = FALSE
      )
      
      result_idx <- result_idx + 1
    }
  }
  
  if (length(annotated_results) == 0) {
    return(results)
  }
  
  result_df <- dplyr::bind_rows(annotated_results)
  return(result_df)
}

#' Convert data.frame to GRanges
#'
#' @param df Data frame with columns: variant_id, chrom (or seqname), start (or pos), end (optional), strand (optional)
#'
#' @return GRanges object
#'
convert_df_to_granges <- function(df) {
  
  # Identify chromosome column name
  chrom_col <- intersect(names(df), c("seqname", "chrom", "chromosome"))[1]
  if (is.na(chrom_col)) stop("Cannot find chromosome column (seqname, chrom, or chromosome)")
  
  # Identify position column name
  pos_col <- intersect(names(df), c("pos", "start", "position"))[1]
  if (is.na(pos_col)) stop("Cannot find position column")
  
  # Check for end column
  if ("end" %in% names(df)) {
    end_col <- "end"
  } else {
    # Use position as both start and end
    end_col <- pos_col
  }
  
  # Convert to GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = df[[chrom_col]],
    ranges = IRanges::IRanges(start = df[[pos_col]], end = df[[end_col]]),
    strand = if ("strand" %in% names(df)) df$strand else "*"
  )
  
  # Add the variant id as the feature name
  names(gr) <- df$variant_id
  return(gr)
}

#' Get Gene Information
#'
#' Retrieves gene symbol and disease gene status for a given gene ID
#' Uses the cached gene_lookup table for fast queries
#'
#' @param gene_id ENSEMBL gene ID
#' @param annotations Annotation list (must include gene_lookup tibble)
#' @param disease_genes Disease gene tibble
#'
#' @return List with gene_symbol and is_disease_gene
#'
get_gene_info <- function(gene_id, annotations, disease_genes = NULL) {
  
  gene_symbol <- NA_character_
  is_disease_gene <- FALSE
  
  # Try to get gene symbol from gene_lookup table
  if (!is.na(gene_id) && !is.null(annotations$gene_lookup)) {
    lookup_row <- annotations$gene_lookup %>%
      dplyr::filter(gene_id == !!gene_id) %>%
      dplyr::slice(1)
    
    if (nrow(lookup_row) > 0) {
      gene_symbol <- lookup_row$gene_symbol[1]
    }
  }
  
  # Check disease gene list
  if (!is.null(disease_genes) && !is.na(gene_id)) {
    is_disease_gene <- gsub("\\.[0-9]+$", "", gene_id) %in% disease_genes$gene_id
  }
  
  return(list(
    gene_symbol = gene_symbol,
    is_disease_gene = is_disease_gene
  ))
}

#' Find Nearest Gene to a Position
#'
#' For intergenic breakpoints, find the nearest gene
#'
#' @param variant GRanges with single breakpoint position
#' @param annotations Annotation list
#'
#' @return List with gene_id and distance
#'
find_nearest_gene <- function(variant, annotations) {
  
  if (!is.data.frame(annotations$genes)) {
    return(list(gene_id = NA_character_, distance = NA_integer_))
  }
  
  genes_gr <- convert_df_to_granges(annotations$genes)
  
  # Find nearest gene
  nearest <- IRanges::nearest(variant, genes_gr)
  
  if (is.na(nearest)) {
    return(list(gene_id = NA_character_, distance = NA_integer_))
  }
  
  nearest_gene <- genes_gr[nearest]
  gene_id <- annotations$genes$gene_id[nearest]
  
  # Calculate distance
  distance <- as.integer(abs(
    ifelse(start(variant) < start(nearest_gene), 
           start(nearest_gene) - start(variant),
           start(variant) - end(nearest_gene))
  ))
  
  return(list(gene_id = gene_id, distance = distance))
}

#' Classify Breakpoint Region Type (Advanced)
#'
#' Determines if breakpoint falls in 5'UTR, 3'UTR, CDS, etc.
#' Requires detailed transcript structure information
#'
#' @param breakpoint_pos Breakpoint position (integer)
#' @param transcript GenomicFeatures transcript object
#' @param cds_ranges CDS ranges for transcript
#' @param utrs_5 5' UTR ranges
#' @param utrs_3 3' UTR ranges
#'
#' @return Character string indicating region type
#'
classify_region_type <- function(breakpoint_pos, transcript, cds_ranges, utrs_5, utrs_3) {
  
  # Check 5' UTR
  if (!is.null(utrs_5) && any(breakpoint_pos >= start(utrs_5) & breakpoint_pos <= end(utrs_5))) {
    return("5_UTR")
  }
  
  # Check CDS
  if (!is.null(cds_ranges) && any(breakpoint_pos >= start(cds_ranges) & breakpoint_pos <= end(cds_ranges))) {
    return("CDS")
  }
  
  # Check 3' UTR
  if (!is.null(utrs_3) && any(breakpoint_pos >= start(utrs_3) & breakpoint_pos <= end(utrs_3))) {
    return("3_UTR")
  }
  
  return("intron")
}
