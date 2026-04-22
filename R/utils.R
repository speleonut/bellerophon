#' Utility Functions for Chromosomal Analysis
#'
#' This module provides utility functions for:
#' - Chromosome normalization (1/chr1 formats)
#' - Handling uncertain breakpoint positions
#' - Variant type classification
#'

#' Normalize Chromosome Names
#'
#' Converts chromosome names to a standard format.
#' Supports input formats: 1-22, X, Y, MT, chr1-chr22, chrX, chrY, chrM, chrMT
#'
#' @param chr Character vector of chromosome names
#' @param style Character string specifying output style: "NCBI" (1,2,...,X,Y,MT) or "UCSC" (chr1,chr2,...,chrX,chrY,chrM)
#'
normalize_chromosome <- function(chr, style = "UCSC") {
  # Handle NULL and return as-is if empty
  if (is.null(chr) || length(chr) == 0) {
    return(chr)
  }
  
  # Convert to character
  chr <- as.character(chr)
  chr <- stringr::str_trim(chr)
  
  # Remove 'chr' prefix if present
  chr_normalized <- stringr::str_replace(chr, "^chr", "")
  
  # Handle MT/M variants
  chr_normalized <- dplyr::case_when(
    chr_normalized %in% c("MT", "M") ~ "MT",
    TRUE ~ chr_normalized
  )
  
  # Apply output style
  if (style == "UCSC") {
    chr_normalized <- dplyr::case_when(
      chr_normalized == "MT" ~ "chrM",
      !stringr::str_detect(chr_normalized, "^chr") ~ paste0("chr", chr_normalized),
      TRUE ~ chr_normalized
    )
  } else if (style == "NCBI") {
    chr_normalized <- stringr::str_replace(chr_normalized, "^chr", "")
    chr_normalized <- dplyr::case_when(
      chr_normalized == "M" ~ "MT",
      is.na(chr_normalized) ~ NA_character_,
      TRUE ~ chr_normalized
    )
  }
  
  return(chr_normalized)
}

#' Handle Uncertain Breakpoint Positions
#'
#' When a breakpoint is represented as a range (e.g., break1_start=100, break1_end=200),
#' this function selects the outermost position relative to the strand orientation.
#'
#' For forward strand (+): use the maximum end position (rightmost)
#' For reverse strand (-): use the minimum start position (leftmost)
#'
#' This approach preserves the widest impact margin for fusion prediction.
#'
#' @param position_start Numeric; start position of uncertain breakpoint
#' @param position_end Numeric; end position of uncertain breakpoint
#' @param orientation Character; "+" (forward) or "-" (reverse) strand
#'
#' @return Numeric; the selected outermost position
#'
handle_uncertain_breakpoint <- function(position_start, position_end, orientation) {
  # If start and end are the same, it's a point breakpoint
  if (is.na(position_start) || is.na(position_end)) {
    return(NA_integer_)
  }
  
  if (position_start == position_end) {
    return(as.integer(position_start))
  }
  
  # Ensure start <= end
  if (position_start > position_end) {
    temp <- position_start
    position_start <- position_end
    position_end <- temp
  }
  
  # Select outermost position based on orientation
  if (orientation == "-") {
    # Reverse strand: use minimum (leftmost) position
    selected_pos <- as.integer(position_start)
  } else {
    # Forward strand or unknown: use maximum (rightmost) position
    selected_pos <- as.integer(position_end)
  }
  
  return(selected_pos)
}

#' Classify Variant Type
#'
#' Classifies a structural variant into one of the standard types:
#' - deletion: segment removed
#' - duplication: tandem copy of segment
#' - insertion: new sequence inserted
#' - inversion: segment inverted
#' - translocation: segments swapped between chromosomes
#' - unclassified: could not determine type
#'
#' @param event_type Character; the reported event type (may be from VCF ALT, MAVIS, or HGVS)
#' @param break1_chromosome Character; chromosome of breakpoint 1
#' @param break2_chromosome Character; chromosome of breakpoint 2
#' @param break1_orientation Character; "+" or "-"
#' @param break2_orientation Character; "+" or "-"
#'
classify_variant_type <- function(event_type,
                                  break1_chromosome,
                                  break2_chromosome,
                                  break1_orientation = NA_character_,
                                  break2_orientation = NA_character_) {
  
  if (is.na(event_type)) {
    event_type <- ""
  }
  
  event_type_lower <- tolower(event_type)
  
  # Explicit event type classification
  if (stringr::str_detect(event_type_lower, "del")) {
    return("deletion")
  } else if (stringr::str_detect(event_type_lower, "dup")) {
    return("duplication")
  } else if (stringr::str_detect(event_type_lower, "inv")) {
    return("inversion")
  } else if (stringr::str_detect(event_type_lower, "ins")) {
    return("insertion")
  } else if (stringr::str_detect(event_type_lower, "bnd|tra")) {
    return("translocation")
  }
  
  # If not explicit, infer from chromosome and orientation
  break1_chromosome <- normalize_chromosome(break1_chromosome)
  break2_chromosome <- normalize_chromosome(break2_chromosome)
  
  if (break1_chromosome != break2_chromosome) {
    return("translocation")
  } else if (break1_orientation == "-" && break2_orientation == "-") {
    return("inversion")
  } else if (!is.na(event_type) && event_type != "") {
    # Try to parse remaining event types
    if (stringr::str_detect(event_type_lower, "unk|complex")) {
      return("unclassified")
    }
  }
  
  return("unclassified")
}

#' Validate Variant Coordinates
#'
#' Checks that variant coordinates are valid (integers, non-negative)
#'
#' @param variant_df Data frame with columns: break1_position_start, break1_position_end,
#'                                            break2_position_start, break2_position_end
#'
#' @return Data frame with additional column: coordinates_valid (logical)
#'
validate_coordinates <- function(variant_df) {
  variant_df %>%
    mutate(
      coordinates_valid = (
        !is.na(break1_position_start) & !is.na(break1_position_end) &
        !is.na(break2_position_start) & !is.na(break2_position_end) &
        break1_position_start > 0 & break1_position_end > 0 &
        break2_position_start > 0 & break2_position_end > 0
      )
    )
}

#' Generate Unique Variant ID
#'
#' Creates a 6-digit hexadecimal unique identifier
#'
#' @param n Number of IDs to generate
#' @param seed Optional seed for reproducibility
#'
generate_variant_ids <- function(n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate random 6-digit hex values
  sprintf("%06X", sample(0:16777215, n, replace = FALSE))
}

#' Check if All Chromosomes Valid
#'
#' @param chr Vector of chromosome names
#'
#' @return Data frame with chromosome and is_valid columns
#'
check_valid_chromosomes <- function(chr) {
  valid_chrs <- as.character(c(1:22, "X", "Y", "MT"))
  
  tibble::tibble(
    chromosome = normalize_chromosome(chr),
    is_valid = normalize_chromosome(chr) %in% valid_chrs
  )
}
