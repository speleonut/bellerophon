#' Input Detection and Parsing
#'
#' This module detects the format of input files (VCF, MAVIS, HGVS)
#' and parses them into a standardized internal format.
#'

#' Detect File Compression Format
#'
#' @param input_file Path to input file
#'
#' @return Character; one of "gzip", "bgzip", "uncompressed"
#'

detect_compression <- function(input_file) {
  if (!file.exists(input_file)) {
    stop(glue::glue("File not found: {input_file}"))
  }
  
  # Check magic bytes
  bytes <- readBin(input_file, "raw", n = 4)
  
  # gzip: 1f 8b
  if (length(bytes) >= 2 && bytes[1] == 0x1f && bytes[2] == 0x8b) {
    return("gzip")
  }
  
  # bgzip is a variant of gzip, will be detected as gzip
  # bgzip additionally has an empty gzip block at the end
  
  return("uncompressed")
}

#' Decompress File if Needed
#'
#' @param input_file Path to input file
#'
#' @return Character; path to uncompressed file (may be original file)
#'
decompress_if_needed <- function(input_file) {
  compression <- detect_compression(input_file)

  if (compression == "uncompressed") {
    return(input_file)
  }

  # Decompress to temporary file
  temp_file <- tempfile()
  R.utils::gunzip(input_file, destname = temp_file, remove = FALSE)

  return(temp_file)
}

#' Detect Input Format
#'
#' Attempts to detect the format of the input file by examining the first few lines.
#'
#' @param input_file Path to input file
#' @param ncbi_contig_mapping_file Path to NCBI contig mapping file (optional)
#'
#' @return List with elements:
#'   - format: Character; "vcf", "mavis", or "hgvs"
#'   - confidence: Numeric; confidence score (0-1)
#'   - first_lines: Character vector; first few lines of file
#'
detect_input_format <- function(input_file, ncbi_contig_mapping_file = NULL) {
  # Decompress if needed
  work_file <- decompress_if_needed(input_file)
  
  # Read first 20 lines
  first_lines <- readLines(work_file, n = 20)
  
  result <- list(
    format = "unknown",
    confidence = 0,
    first_lines = first_lines
  )
  
  # Look for VCF signature
  vcf_header_count <- sum(stringr::str_detect(first_lines, "^##fileformat=VCF"))
  vcf_column_header <- sum(stringr::str_detect(first_lines, "^#CHROM"))
  
  if (vcf_header_count > 0 || vcf_column_header > 0) {
    result$format <- "vcf"
    result$confidence <- 0.95
    return(result)
  }
  
  # Look for MAVIS format (has required columns)
  mavis_required <- c("break1_chromosome", "break1_position_start", "break1_position_end",
                      "break2_chromosome", "break2_position_start", "break2_position_end",
                      "break1_orientation", "break2_orientation", "event_type")
  
  # Check header row
  if (length(first_lines) > 0) {
    header <- first_lines[1]
    
    # Check if all MAVIS required columns are present
    mavis_match_count <- sum(sapply(mavis_required, function(x) stringr::str_detect(header, x)))
    
    if (mavis_match_count >= 8) {  # At least 8 of 9 required columns
      result$format <- "mavis"
      result$confidence <- 0.9
      return(result)
    }
  }
  
  # Look for HGVS format (lines with colon and variant notation)
  hgvs_pattern <- "^[A-Za-z0-9_:.]+(:g\\.|:c\\.|:p\\.)"
  hgvs_count <- sum(stringr::str_detect(first_lines, hgvs_pattern))
  
  if (hgvs_count > 0) {
    result$format <- "hgvs"
    result$confidence <- 0.85
    return(result)
  }
  
  return(result)
}

#' Parse VCF Input
#'
#' Parses VCF format input into standardized variant format.
#' Currently handles simple variants; complex variants are flagged as unclassified.
#'
#' @param input_file Path to VCF file
#'
#' @return Tibble with columns:
#'   - break1_chromosome, break1_position_start, break1_position_end
#'   - break2_chromosome, break2_position_start, break2_position_end
#'   - break1_orientation, break2_orientation
#'   - event_type, source_format, source_line
#'
parse_vcf_input <- function(input_file) {
  # Read VCF using VariantAnnotation
  tryCatch({
    vcf <- VariantAnnotation::readVcf(input_file)
    
    # Extract basic information
    chrom <- as.character(seqnames(vcf))
    pos <- start(vcf)
    ref <- as.character(ref(vcf))
    alt <- as.character(alt(vcf))
    
    
    # So far can process, DEL, DUP, INV.
    # Cannot deal with single BND and translocations yet
    
    variants <- tibble::tibble(
      break1_chromosome = chrom,
      break1_position_start = pos,
      break1_position_end = pos,
      break2_chromosome = chrom,
      break2_position_start = as.numeric(pos) + as.numeric(abs(info(vcf)$SVLEN)-1),
      break2_position_end = as.numeric(pos) + as.numeric(abs(info(vcf)$SVLEN)-1),
      break1_orientation = "+",
      break2_orientation = ifelse(info(vcf)$SVTYPE=="INV", "-", "+"),
      event_type = info(vcf)$SVTYPE,
      source_format = "VCF",
      source_line = seq_along(chrom)
    )
    
    return(variants)
    
  }, error = function(e) {
    stop(glue::glue("Error parsing VCF: {e}"))
  })
}

#' Parse MAVIS Format Input
#'
#' Parses MAVIS TSV/CSV format into standardized variant format.
#'
#' @param input_file Path to MAVIS file (TSV or CSV)
#'
#' @return Tibble with standardized columns
#'
parse_mavis_input <- function(input_file) {
  # Detect delimiter (tab or comma)
  first_line <- readLines(input_file, n = 1)
  delimiter <- if (stringr::str_detect(first_line, "\t")) "\t" else ","
  
  # Read MAVIS format
  mavis_data <- readr::read_delim(
    input_file,
    delim = delimiter,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  # Check for required columns
  required_cols <- c("break1_chromosome", "break1_position_start", "break1_position_end",
                     "break2_chromosome", "break2_position_start", "break2_position_end",
                     "break1_orientation", "break2_orientation", "event_type")
  
  missing_cols <- setdiff(required_cols, names(mavis_data))
  if (length(missing_cols) > 0) {
    stop(glue::glue("MAVIS input missing required columns: {paste(missing_cols, collapse=', ')}"))
  }
  
  # Convert to standardized format
  variants <- mavis_data %>%
    mutate(
      break1_position_start = as.integer(break1_position_start),
      break1_position_end = as.integer(break1_position_end),
      break2_position_start = as.integer(break2_position_start),
      break2_position_end = as.integer(break2_position_end),
      source_format = "MAVIS",
      source_line = row_number()
    ) %>%
    dplyr::select(dplyr::all_of(c(required_cols, "source_format", "source_line")))
  
  return(variants)
}

#' Parse HGVS Format Input
#'
#' Parses HGVS format input (one variant per line).
#' Each line should contain a single variant in HGVS notation.
#'
#' @param input_file Path to HGVS file
#' @param ncbi_contig_mapping_file Path to NCBI contig mapping (optional)
#'
#' @return Tibble with standardized columns
#'
parse_hgvs_input <- function(input_file, ncbi_contig_mapping_file = NULL) {
  # Load NCBI mapping if provided
  if (is.null(ncbi_contig_mapping_file) || !file.exists(ncbi_contig_mapping_file)) {
    mapping <- load_ncbi_contig_mapping()
  } else {
    mapping <- readr::read_tsv(ncbi_contig_mapping_file)
  }
  
  # Read HGVS lines
  hgvs_lines <- readLines(input_file)
  hgvs_lines <- hgvs_lines[stringr::str_trim(hgvs_lines) != ""]  # Remove empty lines
  
  # Parse each HGVS string
  parsed_variants <- list()
  
  for (i in seq_along(hgvs_lines)) {
    hgvs_str <- hgvs_lines[i]
    
    parsed <- parse_hgvs_variant(hgvs_str, mapping = mapping)
    
    # For point variants (both breakpoints at same position),
    # infer second breakpoint as genomically adjacent
    if (parsed$parsed_successfully) {
      variant <- tibble::tibble(
        break1_chromosome = parsed$chromosome,
        break1_position_start = parsed$position_start,
        break1_position_end = parsed$position_end,
        break2_chromosome = parsed$chromosome,
        break2_position_start = parsed$position_start + 1,
        break2_position_end = parsed$position_end + 1,
        break1_orientation = "+",
        break2_orientation = "+",
        event_type = parsed$variant_type,
        source_format = "HGVS",
        source_line = i,
        hgvs_original = hgvs_str
      )
      parsed_variants[[i]] <- variant
    } else {
      # Log parsing error but continue
      warning(glue::glue("Line {i}: {parsed$error_message} (original: {hgvs_str})"))
      parsed_variants[[i]] <- NULL
    }
  }
  
  # Combine all parsed variants
  variants <- dplyr::bind_rows(parsed_variants)
  
  if (nrow(variants) == 0) {
    stop("No valid HGVS variants could be parsed from input file")
  }
  
  return(variants)
}

#' Main Entry Point: Detect and Parse Input
#'
#' @param input_file Path to input file
#' @param input_format Character; "auto", "vcf", "mavis", or "hgvs"
#' @param ncbi_contig_mapping_file Path to NCBI contig mapping (optional)
#'
#' @return List with elements:
#'   - format: Detected format
#'   - variants: Tibble of parsed variants
#'   - warnings: Character vector of warnings
#'
detect_and_parse_input <- function(input_file,
                                   input_format = "auto",
                                   ncbi_contig_mapping_file = NULL) {
  
  if (!file.exists(input_file)) {
    stop(glue::glue("Input file not found: {input_file}"))
  }
  
  warnings <- c()
  
  # Auto-detect format if not specified
  if (input_format == "auto") {
    detection <- detect_input_format(input_file, ncbi_contig_mapping_file)
    input_format <- detection$format
    if (detection$confidence < 0.7) {
      warnings <- c(warnings, glue::glue("Low confidence in format detection (score: {detection$confidence})"))
    }
  }
  
  # Parse based on detected format
  if (input_format == "vcf") {
    variants <- parse_vcf_input(input_file)
  } else if (input_format == "mavis") {
    variants <- parse_mavis_input(input_file)
  } else if (input_format == "hgvs") {
    variants <- parse_hgvs_input(input_file, ncbi_contig_mapping_file)
  } else {
    stop(glue::glue("Unknown input format: {input_format}"))
  }
  
  return(list(
    format = input_format,
    variants = variants,
    warnings = warnings
  ))
}
