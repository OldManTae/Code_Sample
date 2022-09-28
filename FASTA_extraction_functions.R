# Script for extracting fasta sequences as well as introducing variants
####################################################################################################
# Based on motif analysis script
####################################################################################################
##########################
# Download libraries
##########################
library(tidyverse)
library(Biostrings)
# library('BSgenome.Hsapiens.UCSC.hg19')
# library('BSgenome.Hsapiens.UCSC.hg38')

####################################################
# Extract Flanking Sequence of Wt and Variants
# Note: Updated to allow for specification of assembly
####################################################
# given alt sequence, determine offset
# take in list of annotated variants
# note that BSgenome is 1 indexed
get_flanking_sequence <- function(df, flank_dist, reference) {
  # note: expects df to be of the format: chrom, start, end, ref, alt, allele_id
  All_variants <- df %>%
  # Check to see if there is a "chr" - add if absent
  dplyr::mutate(chr = case_when(str_detect(chr, "chr") ~ chr,
                                TRUE ~ paste0("chr", chr)))
  #extract just variants (and add in length of ref and alt)
  variant_sequences_df <- All_variants %>%
    dplyr::select(chr, start, end, ref, alt, allele_id) %>%
    dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    dplyr::mutate(region_length = abs(end-start) + 1, ref_length = str_length(ref), alt_length = ifelse((alt == 0 | alt == "-"), 0, str_length(alt))) %>%
    dplyr::mutate(wt_sequence = paste(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = start-1),
                                      ref,
                                      BSgenome::getSeq(reference, names = chr, start = end + 1,end = end + flank_dist),
                                       sep='')) %>%
    dplyr::mutate(mt_sequence = paste(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = start - 1),
                                      ifelse((alt == 0 | alt == "-"), "", alt),
                                      BSgenome::getSeq(reference, names = chr, start = end + 1, end = end + flank_dist),
                                      sep='')) %>%
    dplyr::mutate(wt_seq_length = str_length(wt_sequence), mut_seq_length = str_length(mt_sequence))
  return(variant_sequences_df)
}

# shortened (one less at the end to allow for even number fragment length)
extract_ref_mut_sequences <- function(df, flank_dist, reference) {
  # note: expects df to be of the format: chrom, start, end, ref, alt, allele_id
  All_variants <- df %>%
    # Check to see if there is a "chr" - add if absent
    dplyr::mutate(chr = case_when(str_detect(chr, "chr") ~ chr,
                                  TRUE ~ paste0("chr", chr)))
  #extract just variants (and add in length of ref and alt)
  variant_sequences_df <- All_variants %>%
    dplyr::select(chr, start, end, ref, alt, allele_id) %>%
    dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    dplyr::mutate(region_length = abs(end-start) + 1, ref_length = str_length(ref), alt_length = ifelse((alt == 0 | alt == "-"), 0, str_length(alt))) %>%
    dplyr::mutate(wt_sequence = paste(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = start - 1),
                                      ref,
                                      BSgenome::getSeq(reference, names = chr, start = end + 1, end = end + flank_dist - 1),
                                      sep='')) %>%
    dplyr::mutate(mt_sequence = paste(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = start - 1),
                                      ifelse((alt == 0 | alt == "-"), "", alt),
                                      BSgenome::getSeq(reference, names = chr, start = end + 1, end = end + flank_dist - 1),
                                      sep='')) %>%
    dplyr::mutate(wt_seq_length = str_length(wt_sequence), mut_seq_length = str_length(mt_sequence)) %>%
    # portion of function to make sure there aren't issues indexing (mutation in wrong location)
    dplyr::mutate(check_sequence = as.character(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = end + flank_dist - 1)),
                  check = wt_sequence == check_sequence)
  return(variant_sequences_df)
}
# extract ref sequences
# shortened (one less at the end to allow for even number fragment length)
extract_ref_sequences <- function(df, flank_dist, reference) {
  # note: expects df to be of the format: chrom, start, end, ref, alt, allele_id
  All_variants <- df %>%
    # Check to see if there is a "chr" - add if absent
    dplyr::mutate(chr = case_when(str_detect(chr, "chr") ~ chr,
                                  TRUE ~ paste0("chr", chr)))
  #extract just variants (and add in length of ref and alt)
  variant_sequences_df <- All_variants %>%
    dplyr::select(chr, start, end, allele_id) %>%
    dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    dplyr::mutate(wt_sequence = as.character(BSgenome::getSeq(reference, names = chr, start = start - flank_dist, end = end + flank_dist - 1))) %>%
    dplyr::mutate(wt_seq_length = str_length(wt_sequence))
    return(variant_sequences_df)
}

