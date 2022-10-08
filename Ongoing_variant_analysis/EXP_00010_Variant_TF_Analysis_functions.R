##################################################################################################v
# function to format annovar output properly... (Annovar output is terrible...)
annovar_output_format <- function(annovar_bed, genome, test = FALSE){
  ## extract SNVs
  # given that annovar output for SNVs are 0 based while we BSgenomes is 1-based
  snvs <- annovar_bed %>%
    dplyr::filter(str_detect(ref, "^[:alpha:]+$") & 
                    str_detect(alt, "^[:alpha:]+$")) 
  
  # get granges for sequence just for ref/alt
  snvs_granges <- as_granges(snvs %>%
                               dplyr::mutate(start = end), seqnames = chr)
  
  snvs_corrected_formatted <- snvs %>%
    # add in ref from reference to check for differences
    dplyr::mutate(ref_sequence = BSgenome::getSeq(x = genome, snvs_granges) %>%
                    as.character(),
                  alt_sequence = alt
    ) %>%
    dplyr::rename(ref_origin = ref,
                  alt_origin = alt) %>%
    dplyr::rename(ref = ref_sequence,
                  alt = alt_sequence)
  
  # as snvs are fine in terms of sequence leave as is
  
  
  ## extract insertions
  ins <- annovar_bed %>%
    dplyr::filter(str_detect(alt, "^[:alpha:]+$") & 
                    (ref == "0" | ref == "-"))
  # for insertions, leave as 0 - based to extract sequence just before insertion
  # only precede if  ins is not empty
  if (plyr::empty(ins) == FALSE){
    # get granges for sequence just before ref/alt
    ins_pre_one_granges <- as_granges(ins %>%
                                        dplyr::mutate(end = start), seqnames = chr)
    
    ins_corrected <- ins %>%
      # add in sequence just before insertion
      dplyr::mutate(ref_sequence = BSgenome::getSeq(x = genome, ins_pre_one_granges) %>%
                      as.character(),
                    alt_sequence = paste0(
                      BSgenome::getSeq(x = genome, ins_pre_one_granges) %>%
                        as.character(), alt
                    )
      ) %>%
      # fix start and end annotation (0-based)
      dplyr::mutate(start_fixed = start - 1)
    
    ins_corrected_formatted <- ins_corrected %>%
      dplyr::rename(ref_origin = ref,
                    alt_origin = alt) %>%
      dplyr::select(chr, start_fixed, end, ref_sequence, alt_sequence, everything(), -c(start)) %>%
      dplyr::rename(start = start_fixed,
                    ref = ref_sequence,
                    alt = alt_sequence)
  } else{
    ins_corrected <- tibble()
    ins_corrected_formatted <- tibble()
  }
  # extract deletions
  del <- annovar_bed %>%
    dplyr::filter(str_detect(ref, "^[:alpha:]+$") & 
                    (alt == "0"| alt == "-"))
  # for deletions, leave as 0 - based to extract sequence just before insertion
  # only continue if deletion tibble not empty
  if (plyr::empty(del) == FALSE){
    # get granges for sequence just before ref/alt
    del_pre_one_granges <- as_granges(del %>%
                                        dplyr::mutate(end = start), seqnames = chr)
    
    del_corrected <- del %>%
      # add in sequence just before deletion
      dplyr::mutate(alt_sequence = BSgenome::getSeq(x = genome, del_pre_one_granges) %>%
                      as.character(),
                    ref_sequence = paste0(
                      BSgenome::getSeq(x = genome, del_pre_one_granges) %>% as.character(), ref
                    )
      ) %>%
      # fix start and end annotation (0-based)
      dplyr::mutate(start_fixed = start - 1)
    
    del_corrected_formatted <- del_corrected %>%
      dplyr::rename(ref_origin = ref,
                    alt_origin = alt) %>%
      dplyr::select(chr, start_fixed, end, ref_sequence, alt_sequence, everything(), -c(start)) %>%
      dplyr::rename(start = start_fixed,
                    ref = ref_sequence,
                    alt = alt_sequence)
  } else{
    del_corrected <- tibble()
    del_corrected_formatted <- tibble()
  }
  # combine into one tibble
  all_variants <- bind_rows(snvs_corrected_formatted,
                            ins_corrected_formatted,
                            del_corrected_formatted)
  
  # return output
  if (test == TRUE) {
    return(list(ins_corrected = ins_corrected, del_corrected = del_corrected, snvs_corrected = snvs_corrected_formatted))
  } else {
    return(all_variants)
  }
}

###################################################################################################################################################
# General function to access enrichment of TFBS impacted by variants
# 3 different analyses
# 1. TF enrichment based on variants that are enriched in cases v controls (nominal excess, not accounting for number affected/unaffected) - must use only for # case varints >= # control variants
# 2. TF enrichment based on affected only variants and not affected only
# 3. TF enrichment based on affected only variants and unaffected only variants (will lose many variants but more stringent?)
# 4. Combined TF enrichment based on counts of strong/weak TF hits across all TFs and aggregating them

# 1. TF enrichment based on variants that are enriched in cases v controls (nominal excess, not accounting for number affected/unaffected) - must use only for # case varints >= # control variants
TF_enrichment_test <- function(motifbreakR_results, variants_with_affected_status){
  # input expectations
  # expects regular motifbreakR results with SNP_id metacolumn
  # expects variants with affected status to info to have columns: name (snp id info), 
  # affected_status (whether the variant is from affected or not affected individual
  # note that filtering of variants should be done prior to inputting into function
  
  # note that the function is assessing whether certain TFBS are disproportionately affected by SNVs 
  # from affected individuals. Each allele is classified as being enriched (nominally) 
  # for patient variants vs unaffected variants. This numerical excess does not take into consideration
  # total number of patient or unaffected individual variants (so need to be careful about imbalance).
  # This function is meant to be used when we have equal or less affected patient variants vs unaffected variants.
  
  # count number of each variants in affected and unaffected individuals. CLassify enrichment
  variants_with_affected_status_counts <- variants_with_affected_status %>%
    dplyr::group_by(name, affected_status) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    pivot_wider(names_from = affected_status, names_prefix = "affected_stat_", values_from = count, values_fill = 0) %>%
    dplyr::mutate(affected_only = case_when(affected_stat_2 > 0 & affected_stat_1 == 0 ~ "affected_only",
                                            TRUE ~ "not_affected_only"),
                  affected_enrichment = affected_stat_2 - affected_stat_1) %>%
    dplyr::select(name, affected_only, affected_enrichment) %>%
    dplyr::mutate(affected_enrichment_status = case_when(affected_enrichment > 0 ~ "affected_enriched",
                                                         TRUE ~ "not_affected_enriched"))
  
  # annotate motifbreakR results
  motifbreakR_results_annotated <- motifbreakR_results %>%
    as_tibble() %>%
    left_join(x = .,
              y = variants_with_affected_status %>%
                dplyr::mutate(name_adjusted = case_when(
                  str_length(ref) > 1 | str_length(alt) > 1 ~ paste0(chr, ":", start + 1,"-", end, ":", ref, ":", alt),
                  TRUE ~ name)),
              by = c("SNP_id" = "name_adjusted")) %>%
    # remove snps that are not annotated in the patient variant list 
    # (useful motifbreakR results was general and needs filtering)
    dplyr::filter(name %in% variants_with_affected_status$name) %>%
    left_join(x = .,
              y = variants_with_affected_status_counts,
              by = c("name" = "name"))
  
  
  # count total number of unique variants (need to count non-motif containing elements too)
  variants_total_counts <- variants_with_affected_status_counts %>%
    dplyr::select(name, affected_enrichment_status) %>%
    distinct() %>%
    dplyr::group_by(affected_enrichment_status) %>%
    dplyr::summarise(total_count =  dplyr::n())
  
  
  # for each TF, make aggregation/count
  motifbreakR_TF_counts <- motifbreakR_results_annotated %>%
    dplyr::group_by(geneSymbol, affected_enrichment_status) %>%
    dplyr::summarise(count =  dplyr::n())
  
  # add totals
  motifbreakR_TF_counts_w_totals <- motifbreakR_TF_counts %>%
    # complete missing categories (if there are no variants, then fill with 0)
    dplyr::ungroup() %>%
    tidyr::complete(geneSymbol, affected_enrichment_status, fill = list(count = 0)) %>%
    left_join(.,
              y = variants_total_counts,
              by = c("affected_enrichment_status" = "affected_enrichment_status")) %>%
    dplyr::mutate(proportion = count/total_count) %>%
    pivot_wider(names_from = affected_enrichment_status,
                values_from = c(proportion, count, total_count),
                values_fill = 0)
  
  # perform hypergeometric test for enrichment - BH adjusted p_values
  motifbreakR_TF_enrichment_test <- motifbreakR_TF_counts_w_totals %>%
    ungroup() %>%
    rowwise() %>%
    dplyr::mutate(
      enrichment_score = log2(
        (count_affected_enriched / total_count_affected_enriched) / 
          ((count_not_affected_enriched + count_affected_enriched) / 
             (total_count_affected_enriched + total_count_not_affected_enriched)))
    ) %>%
    dplyr::mutate(
      p_value = phyper(
        q = count_affected_enriched - 1,
        m = total_count_affected_enriched,
        n = total_count_not_affected_enriched,
        k = count_affected_enriched + count_not_affected_enriched,
        lower.tail = FALSE
      )
    ) %>%
    ungroup() %>%
    adjust_pvalue(p = "p_value", method = "BH")
  return(motifbreakR_TF_enrichment_test)
}


# 2. TF enrichment based on affected only variants and not affected only
TF_enrichment_affected_only_test <- function(motifbreakR_results, variants_with_affected_status){
  # input expectations
  # expects regular motifbreakR results with SNP_id metacolumn
  # expects variants with affected status to info to have columns: name (snp id info), 
  # affected_status (whether the variant is from affected or not affected individual
  # note that filtering of variants should be done prior to inputting into function
  
  # note that the function is assessing whether certain TFBS are disproportionately affected by SNVs 
  # from affected individuals. Each allele is classified as being enriched (nominally) 
  # for patient variants vs unaffected variants. This numerical excess does not take into consideration
  # total number of patient or unaffected individual variants (so need to be careful about imbalance).
  # This function is meant to be used when we have equal or less affected patient variants vs unaffected variants.
  
  # count number of each variants in affected and unaffected individuals. CLassify enrichment
  variants_with_affected_status_counts <- variants_with_affected_status %>%
    dplyr::group_by(name, affected_status) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    pivot_wider(names_from = affected_status, names_prefix = "affected_stat_", values_from = count, values_fill = 0) %>%
    dplyr::mutate(affected_only = case_when(affected_stat_2 > 0 & affected_stat_1 == 0 ~ "affected_only",
                                            TRUE ~ "not_affected_only"),
                  affected_enrichment = affected_stat_2 - affected_stat_1) %>%
    dplyr::select(name, affected_only, affected_enrichment) %>%
    dplyr::mutate(affected_enrichment_status = case_when(affected_enrichment > 0 ~ "affected_enriched",
                                                         TRUE ~ "not_affected_enriched"))
  
  # annotate motifbreakR results
  motifbreakR_results_annotated <- motifbreakR_results %>%
    as_tibble() %>%
    left_join(x = .,
              y = variants_with_affected_status %>%
                dplyr::mutate(name_adjusted = case_when(
                  str_length(ref) > 1 | str_length(alt) > 1 ~ paste0(chr, ":", start + 1,"-", end, ":", ref, ":", alt),
                  TRUE ~ name)),
              by = c("SNP_id" = "name_adjusted")) %>%
    # remove snps that are not annotated in the patient variant list 
    # (useful motifbreakR results was general and needs filtering)
    dplyr::filter(name %in% variants_with_affected_status$name) %>%
    left_join(x = .,
              y = variants_with_affected_status_counts,
              by = c("name" = "name"))
  
  
  # count total number of unique variants (need to count non-motif containing elements too)
  variants_total_counts <- variants_with_affected_status_counts %>%
    dplyr::select(name, affected_only) %>%
    distinct() %>%
    dplyr::group_by(affected_only) %>%
    dplyr::summarise(total_count =  dplyr::n())
  
  
  # for each TF, make aggregation/count
  motifbreakR_TF_counts <- motifbreakR_results_annotated %>%
    dplyr::group_by(geneSymbol, affected_only) %>%
    dplyr::summarise(count =  dplyr::n())
  
  # add totals
  motifbreakR_TF_counts_w_totals <- motifbreakR_TF_counts %>%
    # complete missing categories (if there are no variants, then fill with 0)
    dplyr::ungroup() %>%
    tidyr::complete(geneSymbol, affected_only, fill = list(count = 0)) %>%
    left_join(.,
              y = variants_total_counts,
              by = c("affected_only" = "affected_only")) %>%
    dplyr::mutate(proportion = count/total_count) %>%
    pivot_wider(names_from = affected_only,
                values_from = c(proportion, count, total_count),
                values_fill = 0)
  
  # perform hypergeometric test for enrichment - BH adjusted p_values
  motifbreakR_TF_enrichment_test <- motifbreakR_TF_counts_w_totals %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      enrichment_score = log2(
        (count_affected_only / total_count_affected_only) / 
          ((count_not_affected_only + count_affected_only) /
             (total_count_affected_only + total_count_not_affected_only)))
    ) %>%
    dplyr::mutate(
      p_value = phyper(
        q = count_affected_only - 1,
        m = total_count_affected_only,
        n = total_count_not_affected_only,
        k = count_affected_only + count_not_affected_only,
        lower.tail = FALSE
      )
    ) %>%
    ungroup() %>%
    adjust_pvalue(p = "p_value", method = "BH")
  return(motifbreakR_TF_enrichment_test)
}


# 3. TF enrichment based on affected only variants and unaffected only variants (will lose many variants but more stringent?)
TF_enrichment_affected_unaffected_only_test <- function(motifbreakR_results, variants_with_affected_status){
  # input expectations
  # expects regular motifbreakR results with SNP_id metacolumn
  # expects variants with affected status to info to have columns: name (snp id info), 
  # affected_status (whether the variant is from affected or not affected individual
  # note that filtering of variants should be done prior to inputting into function
  
  # note that the function is assessing whether certain TFBS are disproportionately affected by SNVs 
  # from affected individuals. Each allele is classified as being enriched (nominally) 
  # for patient variants vs unaffected variants. This numerical excess does not take into consideration
  # total number of patient or unaffected individual variants (so need to be careful about imbalance).
  # This function is meant to be used when we have equal or less affected patient variants vs unaffected variants.
  
  # count number of each variants in affected and unaffected individuals
  variants_with_affected_status_counts <- variants_with_affected_status %>%
    dplyr::group_by(name, affected_status) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    pivot_wider(names_from = affected_status, names_prefix = "affected_stat_", values_from = count, values_fill = 0) %>%
    dplyr::mutate(affected_only = case_when(affected_stat_2 > 0 & affected_stat_1 == 0 ~ "affected_only",
                                            affected_stat_1 > 0 & affected_stat_2 == 0 ~ "unaffected_only",
                                            TRUE ~ "shared"),
                  affected_enrichment = affected_stat_2 - affected_stat_1) %>%
    # remove variants that are in both affected and unaffected individuals
    dplyr::filter(affected_only != "shared") %>%
    dplyr::select(name, affected_only, affected_enrichment) %>%
    dplyr::mutate(affected_enrichment_status = case_when(affected_enrichment > 0 ~ "affected_enriched",
                                                         TRUE ~ "not_affected_enriched"))
  
  # annotate motifbreakR results
  motifbreakR_results_annotated <- motifbreakR_results %>%
    as_tibble() %>%
    left_join(x = .,
              y = variants_with_affected_status %>%
                dplyr::mutate(name_adjusted = case_when(
                  str_length(ref) > 1 | str_length(alt) > 1 ~ paste0(chr, ":", start + 1,"-", end, ":", ref, ":", alt),
                  TRUE ~ name)),
              by = c("SNP_id" = "name_adjusted")) %>%
    # remove snps that are not annotated in the patient variant list 
    # (useful motifbreakR results was general and needs filtering)
    dplyr::filter(name %in% variants_with_affected_status$name) %>%
    left_join(x = .,
              y = variants_with_affected_status_counts,
              by = c("name" = "name"))
  
  
  # count total number of unique variants (need to count non-motif containing elements too)
  variants_total_counts <- variants_with_affected_status_counts %>%
    dplyr::select(name, affected_only) %>%
    distinct() %>%
    dplyr::group_by(affected_only) %>%
    dplyr::summarise(total_count =  dplyr::n())
  
  
  # for each TF, make aggregation/count
  motifbreakR_TF_counts <- motifbreakR_results_annotated %>%
    dplyr::group_by(geneSymbol, affected_only) %>%
    dplyr::summarise(count =  dplyr::n())
  
  # add totals
  motifbreakR_TF_counts_w_totals <- motifbreakR_TF_counts %>%
    # complete missing categories (if there are no variants, then fill with 0)
    dplyr::ungroup() %>%
    tidyr::complete(geneSymbol, affected_only, fill = list(count = 0)) %>%
    left_join(.,
              y = variants_total_counts,
              by = c("affected_only" = "affected_only")) %>%
    dplyr::mutate(proportion = count/total_count) %>%
    pivot_wider(names_from = affected_only,
                values_from = c(proportion, count, total_count),
                values_fill = 0)
  
  # perform hypergeometric test for enrichment - BH adjusted p_values
  motifbreakR_TF_enrichment_test <- motifbreakR_TF_counts_w_totals %>%
    ungroup() %>%
    rowwise() %>%
    dplyr::mutate(
      enrichment_score = log2(
        (count_affected_only / total_count_affected_only) / 
          ((count_unaffected_only + count_affected_only) / 
             (total_count_affected_only + total_count_unaffected_only)))
    ) %>%
    dplyr::mutate(
      p_value = phyper(
        q = count_affected_only - 1,
        m = total_count_affected_only,
        n = total_count_unaffected_only,
        k = count_affected_only + count_unaffected_only,
        lower.tail = FALSE
      )
    ) %>%
    ungroup() %>%
    adjust_pvalue(p = "p_value", method = "BH")
  return(motifbreakR_TF_enrichment_test)
}


# 4. Combined TF enrichment based on counts of strong/weak TF hits across all TFs and aggregating them
TF_aggregate_enrichment_test <- function(motifbreakR_results, variants_with_affected_status){
  # input expectations
  # expects regular motifbreakR results with SNP_id metacolumn
  # expects variants with affected status to info to have columns: name (snp id info), 
  # affected_status (whether the variant is from affected or not affected individual
  # note that filtering of variants should be done prior to inputting into function
  
  # note that the function is assessing whether certain TFBS are disproportionately affected by SNVs 
  # from affected individuals. Each allele is classified as being enriched (nominally) 
  # for patient variants vs unaffected variants. This numerical excess does not take into consideration
  # total number of patient or unaffected individual variants (so need to be careful about imbalance).
  # This function is meant to be used when we have equal or less affected patient variants vs unaffected variants.
  
  # count number of each variants in affected and unaffected individuals. CLassify enrichment
  variants_with_affected_status_counts <- variants_with_affected_status %>%
    dplyr::group_by(name, affected_status) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    pivot_wider(names_from = affected_status, names_prefix = "affected_stat_", values_from = count, values_fill = 0) %>%
    dplyr::mutate(affected_only = case_when(affected_stat_2 > 0 & affected_stat_1 == 0 ~ "affected_only",
                                            TRUE ~ "not_affected_only"),
                  affected_enrichment = affected_stat_2 - affected_stat_1) %>%
    dplyr::select(name, affected_only, affected_enrichment) %>%
    dplyr::mutate(affected_enrichment_status = case_when(affected_enrichment > 0 ~ "affected_enriched",
                                                         TRUE ~ "not_affected_enriched"))
  
  # annotate motifbreakR results
  motifbreakR_results_annotated <- motifbreakR_results %>%
    as_tibble() %>%
    left_join(x = .,
              y = variants_with_affected_status %>%
                dplyr::mutate(name_adjusted = case_when(
                  str_length(ref) > 1 | str_length(alt) > 1 ~ paste0(chr, ":", start + 1,"-", end, ":", ref, ":", alt),
                  TRUE ~ name)),
              by = c("SNP_id" = "name_adjusted")) %>%
    # remove snps that are not annotated in the patient variant list 
    # (useful motifbreakR results was general and needs filtering)
    dplyr::filter(name %in% variants_with_affected_status$name) %>%
    left_join(x = .,
              y = variants_with_affected_status_counts,
              by = c("name" = "name"))
  
  
  # count total number of unique variants (need to count non-motif containing elements too)
  variants_total_counts <- variants_with_affected_status_counts %>%
    dplyr::select(name, affected_enrichment_status) %>%
    distinct() %>%
    dplyr::group_by(affected_enrichment_status) %>%
    dplyr::summarise(total_count =  dplyr::n())
  
  
  # for each TF, make aggregation/count (only 1 TF)
  motifbreakR_TF_aggregate_counts <- motifbreakR_results_annotated %>%
    dplyr::group_by(affected_enrichment_status, name) %>%
    dplyr::summarise(count =  dplyr::n()) %>%
    dplyr::group_by(affected_enrichment_status) %>%
    dplyr::summarise(count = dplyr::n())
  
  # add totals
  motifbreakR_TF_counts_w_totals <- motifbreakR_TF_aggregate_counts %>%
    # complete missing categories (if there are no variants, then fill with 0)
    dplyr::ungroup() %>%
    tidyr::complete(affected_enrichment_status, fill = list(count = 0)) %>%
    left_join(.,
              y = variants_total_counts,
              by = c("affected_enrichment_status" = "affected_enrichment_status")) %>%
    dplyr::mutate(proportion = count/total_count) %>%
    pivot_wider(names_from = affected_enrichment_status,
                values_from = c(proportion, count, total_count),
                values_fill = 0)
  
  # perform hypergeometric test for enrichment - BH adjusted p_values
  motifbreakR_TF_enrichment_test <- motifbreakR_TF_counts_w_totals %>%
    ungroup() %>%
    rowwise() %>%
    dplyr::mutate(
      enrichment_score = log2(
        (count_affected_enriched / total_count_affected_enriched) / 
          ((count_not_affected_enriched + count_affected_enriched) / 
             (total_count_affected_enriched + total_count_not_affected_enriched)))
    ) %>%
    dplyr::mutate(
      p_value = phyper(
        q = count_affected_enriched - 1,
        m = total_count_affected_enriched,
        n = total_count_not_affected_enriched,
        k = count_affected_enriched + count_not_affected_enriched,
        lower.tail = FALSE
      )
    ) %>%
    ungroup() %>%
    adjust_pvalue(p = "p_value", method = "BH")
  return(motifbreakR_TF_enrichment_test)
}

