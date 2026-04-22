#' Annotation Loading - Phase 2
#'
#' This module handles loading:
#' 1. GENCODE transcript annotations download directly using rtracklayer
#' 2. Disease gene list
#' 3. GTEx tissue metadata and expression data
#'

#' Load GENCODE Annotations via rtracklayer
#'
#' Downloads or loads cached GENCODE GFF as specified in fusion_pipeline.qmd for hg38 and converts to TxDb format.
#' Creates GRanges objects for fast exon/intron/CDS queries.
#'
#' @param force_download Logical; if TRUE, ignore cached files and re-download
#'
#' @return List containing:
#'   - txdb: GenomicFeatures::TxDb object
#'   - genes: GRanges object of genes
#'   - transcripts: Granges object of transcripts
#'   - exons: GRanges object of exon locations
#'   - introns: GRanges object of intron locations
#'   - cds: GRanges object of coding sequence regions
#'

load_gencode_annotations <- function(force_download = FALSE, gencode_release = 49) {
  
  cat("Loading GENCODE annotations...\n")
  
  tryCatch({
    # Check for a local copy of the gff3 in data/annotation
    # This is perhaps not the best way to handle this but it is more reliable than AnnotationHub or bioMaRt
    gtf_file <- paste("data/annotations/gencode.v",gencode_release,".basic.annotation.gff3.gz", sep = "")
    
    if (file_test("-f", gtf_file)) {
      cat("  ✓ Found local GENCODE GTF file: ",gtf_file, "\n")
      gencode_gtf <- rtracklayer::import(gtf_file) # cache the file locally if you like
    } else {
      gtf_url <- paste("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",gencode_release,"/gencode.v",gencode_release,".basic.annotation.gff3.gz", sep = "")
      gencode_gtf <- rtracklayer::import(gtf_url)
    }
    
    cat("  ✓ Imported GENCODE GTF file\n")
    
    # Convert GFF to TxDb using txdbmaker (makeTxDbFromGRanges moved to txdbmaker package)
    txdb <- txdbmaker::makeTxDbFromGRanges(gencode_gtf)
    cat("  ✓ Converted to TxDb format\n")
    
    # Extract gene lookup table from original GTF (before TxDb loses gene_name metadata)
    # This is more memory efficient than storing full gene tibbles
    gene_mask <- !is.na(gencode_gtf$type) & gencode_gtf$type == "gene"
    genes_from_gtf <- gencode_gtf[gene_mask]
    
    gene_lookup <- tibble::tibble(
      gene_id = genes_from_gtf$gene_id,
      gene_symbol = genes_from_gtf$gene_name,
      seqname = as.character(GenomicRanges::seqnames(genes_from_gtf)),
      start = BiocGenerics::start(genes_from_gtf),
      end = BiocGenerics::end(genes_from_gtf),
      strand = as.character(GenomicRanges::strand(genes_from_gtf))
    ) %>%
      distinct(gene_id, .keep_all = TRUE)
    
    cat("  ✓ Created gene lookup table: ", nrow(gene_lookup), " genes\n", sep = "")
    
    genes_gr <- genes(txdb)
    
    # Note: GENCODE GFF includes gene_name in metadata
    cat("  ✓ Extracted gene coordinates: ", length(genes_gr), " genes\n", sep = "")
    
    # Extract transcript information from txdb and make a lookup table mapping "ENSG" gene_id to "ENST" tx_id
    tx_gr <- GenomicFeatures::transcriptsBy(txdb, by = "gene") %>% 
      unlist(use.names = TRUE) # use.names = TRUE uses the gene_id as names
    
    # To make the lookup simply extract names and tx_name from tx_gr
    tx_lookup <- tibble::tibble(
      tx_id = tx_gr$tx_name,
      gene_id = names(tx_gr))

    cat("  ✓ Created transcript lookup table and transcript coordinates for: ", nrow(tx_lookup), " transcripts\n", sep = "")
    

    # Get exons from txdb and use ENST transcript IDs in names for easy lookup when annotating breakpoints
    exons_gr <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE) %>%
      unlist(use.names = TRUE) # With use.names = TRUE, names(exons_gr) is the ENST id
    
    # e.g.
    # GRanges object with 2525461 ranges and 3 metadata columns:
    #                  seqnames      ranges strand |   exon_id           exon_name exon_rank
    #                      <Rle>   <IRanges>  <Rle> | <integer>         <character> <integer>
    # ENST00000832828.1     chr1 11426-11671      + |         1 ENST00000832828.1:1         1
    # ENST00000832828.1     chr1 12010-12227      + |         3 ENST00000832828.1:2         2

       
    cat("  ✓ Extracted exons as GRanges: ", length(exons_gr), " exons\n", sep = "")
    
    # Get introns from txdb and use ENST transcript IDs in names for easy lookup when annotating breakpoints
    introns_gr <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE) %>%
      unlist(use.names = TRUE) # Same logic as above transcript id are in names()
    
    cat("  ✓ Extracted introns as GRanges\n")
    
    # Get CDS regions directly from gencode_gtf (with phase information)
    # CDS regions from txdb lose phase information, so extract from original GTF
    cds_mask <- !is.na(gencode_gtf$type) & gencode_gtf$type == "CDS"
    cds_from_gtf <- gencode_gtf[cds_mask]
    cds_gr <- cds_from_gtf
    
    # Add transcript_id to mcols for easy lookup
    mcols(cds_gr)$tx_id <- cds_from_gtf$transcript_id
    
    # Calculate end_phase for reading frame analysis
    # end_phase = (phase + width) mod 3
    if ("phase" %in% colnames(mcols(cds_gr))) {
      mcols(cds_gr)$end_phase <- (mcols(cds_gr)$phase + BiocGenerics::width(cds_gr)) %% 3
    } else {
      # Fallback if phase not available
      cat("  Warning: phase column not found in CDS; reading frame analysis may be limited\n")
    }
    
    # Use transcript IDs in names for easy lookup
    names(cds_gr) <- cds_from_gtf$transcript_id

    cat("  ✓ Extracted CDS regions as GRanges (with phase information)\n")
    
    result <- list(
      txdb = txdb,
      genes = genes_gr,
      gene_lookup = gene_lookup,
      transcripts = tx_gr,
      tx_lookup = tx_lookup,
      exons = exons_gr,
      introns = introns_gr,
      cds = cds_gr,
      loaded_successfully = TRUE
    )
    
    cat("✓ GENCODE annotations loaded successfully\n")
    return(result)
    
  }, error = function(e) {
    cat("✗ Error loading GENCODE annotations: ", as.character(e), "\n", sep = "")
    return(list(loaded_successfully = FALSE, error = as.character(e)))
  })
}

#' Load Disease Gene List
#'
#' Reads a file containing ENSEMBL gene IDs of genes associated with genetic diseases.
#' Expected file format: one ENSG ID per line (e.g., ENSG00000000003)
#'
#' @param disease_gene_file Path to disease gene list file
#'
#' @return Tibble with columns: gene_id, is_disease_gene (TRUE)
#'
load_disease_genes <- function(disease_gene_file = "data/annotations/Nijmegen.DG.ENSG.list.txt") {
  
  if (!file.exists(disease_gene_file)) {
    cat("Disease gene list not found at: ", disease_gene_file, "\n", sep = "")
    cat("Creating empty disease gene list\n")
    return(tibble::tibble(gene_id = character()))
  }
  
  tryCatch({
    # Read file - expect one ENSG ID per line
    disease_genes <- readLines(disease_gene_file) %>%
      stringr::str_trim() %>%
      # Filter out empty lines and comments
      .[!grepl("^#|^$", .)] %>%
      # Extract ENSG IDs (format: ENSGxxxxxxxxxx)
      stringr::str_extract("ENSG\\d+") %>%
      unique()
    
    result <- tibble::tibble(
      gene_id = disease_genes,
      is_disease_gene = TRUE
    )
    
    cat("✓ Loaded ", nrow(result), " disease genes\n", sep = "")
    return(result)
    
  }, error = function(e) {
    cat("✗ Error loading disease genes: ", as.character(e), "\n", sep = "")
    return(tibble::tibble(gene_id = character()))
  })
}

#' Initialize GTEx Metadata
#'
#' Retrieves available tissues and samples from GTEx portal.
#' Builds a mapping of tissue IDs to tissue names for later expression queries.
#'
#' @return List containing:
#'   - tissues: tibble of tissues (tissue_id, tissue_name, tissue_category)
#'   - gtex_gene_median_tpm: tibble of gtex medians (from file or download as needed)
#'   - loaded_successfully: logical
#'
initialize_gtex_metadata <- function(gtex_file = "data/annotations/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz") {
  
  cat("Initializing GTEx metadata...\n")
  
  tryCatch({
    # Note: gtexr package requires specific Gencode ID with a version number which is difficult to know without downloading all of the data.
    # Might as well just download the data and then strip version numbers locally when needed to compare to the annotation data
    
    gtex_tissues <- tibble::tibble(
      tissue_id = c(
        "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
        "Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
        "Bladder", "Blood", "Blood_Vessel",
        "Brain_Amygdala", "Brain_Anterior_cingulate_cortex", "Brain_Caudate",
        "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex",
        "Brain_Frontal_Cortex", "Brain_Hippocampus", "Brain_Hypothalamus",
        "Brain_Nucleus_accumbens", "Brain_Putamen", "Brain_Spinal_cord",
        "Brain_Substantia_nigra", "Breast", "Cells_EBV",
        "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse",
        "Cornea", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa",
        "Esophagus_Muscularis", "Fallopian_Tube", "Heart_Atrial_Appendage",
        "Heart_Left_Ventricle", "Kidney_Cortex", "Kidney_Medulla",
        "Liver", "Lung", "Minor_Salivary_Gland",
        "Muscle_Skeletal", "Nerve_Tibial", "Ovary",
        "Pancreas", "Pituitary", "Prostate",
        "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
        "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach",
        "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"
      ),
      tissue_name = c(
        "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland",
        "Artery - Aorta", "Artery - Coronary", "Artery - Tibial",
        "Bladder", "Blood", "Blood Vessel",
        "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate",
        "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex",
        "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus",
        "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)",
        "Brain - Substantia nigra", "Breast - Mammary Tissue", "Cells - EBV-transformed lymphocytes",
        "Cells - Transformed fibroblasts", "Colon - Sigmoid", "Colon - Transverse",
        "Cornea", "Esophagus - Gastroesophageal Junction", "Esophagus - Mucosa",
        "Esophagus - Muscularis", "Fallopian Tube", "Heart - Atrial Appendage",
        "Heart - Left Ventricle", "Kidney - Cortex", "Kidney - Medulla",
        "Liver", "Lung", "Minor Salivary Gland",
        "Muscle - Skeletal", "Nerve - Tibial", "Ovary",
        "Pancreas", "Pituitary", "Prostate",
        "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
        "Small Intestine - Terminal Ileum", "Spleen", "Stomach",
        "Testis", "Thyroid", "Uterus", "Vagina", "Whole Blood"
      ),
      tissue_category = c(
        "Adipose", "Adipose", "Endocrine",
        "Blood Vessel", "Blood Vessel", "Blood Vessel",
        "Urinary", "Blood", "Blood Vessel",
        "Brain", "Brain", "Brain",
        "Brain", "Brain", "Brain",
        "Brain", "Brain", "Brain",
        "Brain", "Brain", "Brain",
        "Brain", "Breast", "Cells",
        "Cells", "Colon", "Colon",
        "Eye", "Esophagus", "Esophagus",
        "Esophagus", "Reproductive", "Heart",
        "Heart", "Kidney", "Kidney",
        "Liver", "Lung", "Salivary Gland",
        "Muscle", "Nerve", "Reproductive",
        "Pancreas", "Pituitary", "Reproductive",
        "Skin", "Skin",
        "Small Intestine", "Spleen", "Stomach",
        "Reproductive", "Thyroid", "Reproductive", "Reproductive", "Blood"
      )
    )
    
    cat("  ✓ Initialized ", nrow(gtex_tissues), " GTEx tissues\n", sep = "")
    
    if (file_test(gtex_file)) {
      cat("  ✓ Local copy of GTEx median TPMs located")
      gtex_median_tpm <- read_tsv(gtex_file, skip = 2, lazy = TRUE)
    } else {
      gtex_url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
      gtex_median_tpm <- read_tsv(gtex_url, skip = 2)
    }
    
    result <- list(
      tissues = gtex_tissues,
      gtex_median_tpm = gtex_median_tpm,
      loaded_successfully = TRUE
    )
    
    return(result)
    
  }, error = function(e) {
    cat("✗ Error initializing GTEx metadata: ", as.character(e), "\n", sep = "")
    return(list(loaded_successfully = FALSE, error = as.character(e)))
  })
}


