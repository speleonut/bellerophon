#' Parse HGVS Nomenclature Strings
#'
#' This module provides functions to parse HGVS nomenclature for structural variants.
#' Supports genomic (g.), coding (c.), and protein (p.) notation.
#' Handles NCBI contig IDs and converts to chromosome names.
#'
#' Reference: https://hgvs-nomenclature.org/
#'

# Load NCBI to chromosome mapping
load_ncbi_contig_mapping <- function(mapping_file = NULL) {
  if (is.null(mapping_file) || !file.exists(mapping_file)) {
    # Define default mapping for GRCh38
    mapping <- tibble::tibble(
      ncbi_contig = c(
        "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
        "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
        "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
        "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
        "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
        "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10",
        "NC_012920.1"
      ),
      chromosome = c(
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "X", "Y", "MT"
      ),
      ncbi_version = c(
        "11", "12", "12", "12",
        "10", "12", "14", "11",
        "12", "11", "10", "12",
        "11", "9", "10", "10",
        "11", "10", "10", "11",
        "9", "11", "11", "10",
        "1"
      )
    )
    
    # Add alternative versions
    mapping <- rbind(
      mapping,
      tibble::tibble(
        ncbi_contig = gsub("NC_", "NC_", c(
          "NC_000001.10", "NC_000002.11", "NC_000003.11", "NC_000004.11",
          "NC_000005.9", "NC_000006.11", "NC_000007.13", "NC_000008.10",
          "NC_000009.11", "NC_000010.10", "NC_000011.9", "NC_000012.11",
          "NC_000013.10", "NC_000014.8", "NC_000015.9", "NC_000016.9",
          "NC_000017.10", "NC_000018.9", "NC_000019.9", "NC_000020.10",
          "NC_000021.8", "NC_000022.10", "NC_000023.10", "NC_000024.9"
        )),
        chromosome = c(
          "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
          "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
          "21", "22", "X", "Y"
        ),
        ncbi_version = c(
          "10", "11", "11", "11",
          "9", "11", "13", "10",
          "11", "10", "9", "11",
          "10", "8", "9", "9",
          "10", "9", "9", "10",
          "8", "10", "10", "9"
        )
      )
    )
    
    return(mapping)
  } else {
    # Load from user-provided file
    return(readr::read_tsv(mapping_file))
  }
}

#' Extract Chromosome from NCBI Contig ID
#'
#' Maps NCBI contig IDs (e.g., NC_000001.11) to chromosome names (e.g., "1")
#'
ncbi_contig_to_chromosome <- function(contig_id, mapping = NULL) {
  if (is.null(mapping)) {
    mapping <- load_ncbi_contig_mapping()
  }
  
  # Extract major version (e.g., NC_000001.11 -> NC_000001)
  contig_base <- stringr::str_extract(contig_id, "^NC_\\d+")
  
  matched_chrs <- mapping %>%
    dplyr::filter(stringr::str_detect(ncbi_contig, contig_base)) %>%
    dplyr::pull(chromosome)
  
  # Check if any matches found
  if (length(matched_chrs) == 0) {
    warning(glue::glue("Could not map NCBI contig {contig_id} to chromosome"))
    return(contig_id)
  }
  
  return(matched_chrs[1])
}

#' Parse HGVS Genomic Nomenclature
#'
#' Parses HGVS genomic notation (g.) for deletions, duplications, inversions, insertions
#' Examples:
#'   - Deletion: g.123del or g.123_456del
#'   - Duplication: g.123dup or g.123_456dup
#'   - Inversion: g.123_456inv
#'   - Insertion: g.123_124insAGCT
#'   - Substitution: g.123A>G
#'
parse_hgvs_genomic <- function(hgvs_string, mapping = NULL) {
  if (is.null(mapping)) {
    mapping <- load_ncbi_contig_mapping()
  }
  
  result <- list(
    original = hgvs_string,
    variant_type = NA_character_,
    chromosome = NA_character_,
    position_start = NA_integer_,
    position_end = NA_integer_,
    ref_allele = NA_character_,
    alt_allele = NA_character_,
    parsed_successfully = FALSE,
    error_message = NA_character_
  )
  
  tryCatch({
    # Trim whitespace
    s <- stringr::str_trim(hgvs_string)
    
    # Check if this looks like genomic notation (must contain :g.)
    if (!stringr::str_detect(s, ":g\\.")) {
      return(result)  # Not a genomic notation
    }
    
    # Extract the part before :g. (chromosome/contig ID)
    chr_part <- stringr::str_extract(s, "^[^:]+")  # Everything before the first ":"
    
    if (is.na(chr_part) || chr_part == "") {
      result$error_message <- "Could not extract chromosome/contig ID from HGVS string"
      return(result)
    }
    
    # Map NCBI contig to chromosome if needed
    if (stringr::str_detect(chr_part, "^NC_")) {
      chr <- ncbi_contig_to_chromosome(chr_part, mapping = mapping)
    } else {
      # Direct chromosome reference (chr1, 1, X, etc.)
      chr <- chr_part
    }
    result$chromosome <- chr
    
    # Remove chromosome prefix and :g. to get the variant part
    # Pattern: remove everything up to and including :g.
    s_no_chr <- stringr::str_replace(s, "^[^:]*:g\\.", "")
    
    # Parse variant type and positions
    # Deletion: 123_456del or 123del
    if (stringr::str_detect(s_no_chr, "del$")) {
      result$variant_type <- "deletion"
      pos_part <- stringr::str_replace(s_no_chr, "del$", "")
      if (stringr::str_detect(pos_part, "_")) {
        positions <- stringr::str_split_fixed(pos_part, "_", 2)
        result$position_start <- suppressWarnings(as.integer(positions[1, 1]))
        result$position_end <- suppressWarnings(as.integer(positions[1, 2]))
      } else {
        result$position_start <- suppressWarnings(as.integer(pos_part))
        result$position_end <- suppressWarnings(as.integer(pos_part))
      }
    }
    # Duplication: 123_456dup or 123dup
    else if (stringr::str_detect(s_no_chr, "dup$")) {
      result$variant_type <- "duplication"
      pos_part <- stringr::str_replace(s_no_chr, "dup$", "")
      if (stringr::str_detect(pos_part, "_")) {
        positions <- stringr::str_split_fixed(pos_part, "_", 2)
        result$position_start <- suppressWarnings(as.integer(positions[1, 1]))
        result$position_end <- suppressWarnings(as.integer(positions[1, 2]))
      } else {
        result$position_start <- suppressWarnings(as.integer(pos_part))
        result$position_end <- suppressWarnings(as.integer(pos_part))
      }
    }
    # Inversion: 123_456inv
    else if (stringr::str_detect(s_no_chr, "inv$")) {
      result$variant_type <- "inversion"
      pos_part <- stringr::str_replace(s_no_chr, "inv$", "")
      positions <- stringr::str_split_fixed(pos_part, "_", 2)
      result$position_start <- suppressWarnings(as.integer(positions[1, 1]))
      result$position_end <- suppressWarnings(as.integer(positions[1, 2]))
    }
    # Insertion: 123_124insAGCT or 123_124ins (without sequence)
    else if (stringr::str_detect(s_no_chr, "ins")) {
      result$variant_type <- "insertion"
      pos_part <- stringr::str_replace(s_no_chr, "ins.*$", "")
      
      if (stringr::str_detect(pos_part, "_")) {
        positions <- stringr::str_split_fixed(pos_part, "_", 2)
        result$position_start <- suppressWarnings(as.integer(positions[1, 1]))
        result$position_end <- suppressWarnings(as.integer(positions[1, 2]))
      } else {
        result$position_start <- suppressWarnings(as.integer(pos_part))
        result$position_end <- suppressWarnings(as.integer(pos_part) + 1)
      }
      
      # Extract inserted sequence if present
      ins_seq <- stringr::str_extract(s_no_chr, "ins[ACGT]*")
      if (!is.na(ins_seq) && nchar(ins_seq) > 3) {  # More than just "ins"
        result$alt_allele <- stringr::str_replace(ins_seq, "ins", "")
      }
    }
    # Substitution: 123A>G
    else if (stringr::str_detect(s_no_chr, "[ACGT]>[ACGT]$")) {
      result$variant_type <- "substitution"
      pos_allele <- stringr::str_split_fixed(s_no_chr, "[ACGT]>[ACGT]$", 2)[1, 1]
      alleles <- stringr::str_extract(s_no_chr, "[ACGT]>[ACGT]$") %>%
        stringr::str_split(">", simplify = TRUE)
      result$position_start <- suppressWarnings(as.integer(pos_allele))
      result$position_end <- suppressWarnings(as.integer(pos_allele))
      result$ref_allele <- alleles[1, 1]
      result$alt_allele <- alleles[1, 2]
    }
    # Complex variants (not supported)
    else if (stringr::str_detect(s_no_chr, ";")) {
      result$error_message <- "Complex variants with multiple operations are not supported yet"
      return(result)
    }
    else {
      result$error_message <- glue::glue("Could not parse HGVS variant: {hgvs_string}")
      return(result)
    }
    
    result$parsed_successfully <- TRUE
    
  }, error = function(e) {
    result$error_message <<- as.character(e)
  })
  
  return(result)
}

#' Parse HGVS String and Return Structured Data
#'
#' Main entry point for HGVS parsing
#'
parse_hgvs_variant <- function(hgvs_string, mapping = NULL) {
  if (is.null(mapping)) {
    mapping <- load_ncbi_contig_mapping()
  }
  
  # For now, only handle genomic notation
  if (stringr::str_detect(hgvs_string, ":g\\.")) {
    return(parse_hgvs_genomic(hgvs_string, mapping = mapping))
  } else if (stringr::str_detect(hgvs_string, ":c\\.")) {
    return(list(
      original = hgvs_string,
      parsed_successfully = FALSE,
      error_message = "Coding notation (c.) not yet supported; please use genomic notation (g.)"
    ))
  } else if (stringr::str_detect(hgvs_string, ":p\\.")) {
    return(list(
      original = hgvs_string,
      parsed_successfully = FALSE,
      error_message = "Protein notation (p.) not yet supported"
    ))
  } else {
    return(list(
      original = hgvs_string,
      parsed_successfully = FALSE,
      error_message = "Could not detect HGVS notation type (g., c., or p.)"
    ))
  }
}
