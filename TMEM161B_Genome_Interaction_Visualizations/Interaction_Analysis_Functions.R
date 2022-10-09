# function to extract interactions that overlap with a gene's TSS - with wiggle room due to imprecision of interaction methods
split_interactions_by_geneTSS <- function(interactions, gene_symbol, padding = 0){
  # extract TSS for gene
  print("loading hg19 genome...")
  # import hg19 genome
  ensembl <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  print("extracting gene TSS...")
  # extract gene TSS
  gene_TSS <- biomaRt::select(ensembl, 
                              keys = gene_symbol,
                              columns = c('external_gene_name','chromosome_name', 'transcription_start_site'),
                              keytype = "external_gene_name") %>%
    dplyr::mutate(transcription_start_site_w_padding = transcription_start_site + padding,
                  chromosome_name = paste("chr", chromosome_name, sep = ""))
  
  ## check read 1
  # find intersection between gene_TSS and interactions
  print("checking read 1 for TSS overlap...")
  interactions_TSS_intersect_R1 <- genome_intersect(x = interactions, y = gene_TSS, by = c("chr1" = "chromosome_name", "start1" = "transcription_start_site", "end1" = "transcription_start_site_w_padding"), mode = "left")
  
  ## check read 2
  # find intersection between gene_TSS and interactions
  print("checking read 2 for TSS overlap...")
  interactions_TSS_intersect_R2 <- genome_intersect(x = interactions, y = gene_TSS, by = c("chr2" = "chromosome_name", "start2" = "transcription_start_site", "end2" = "transcription_start_site_w_padding"), mode = "left")
  # combine data and format data
  interactions_TSS_intersect <- bind_rows(interactions_TSS_intersect_R1,interactions_TSS_intersect_R2) %>%
    dplyr::mutate(name = row_number(),
                  score = 1,
                  strand1 = "*",
                  strand2= "*") %>%
    dplyr::select(chr1, start1, end1, chr2, start2, end2, name, score, strand1, strand2, everything())
  
  return(interactions_TSS_intersect)
}

# this function generates bedpe from interactions annotated by gene symbol
#eg SCN3A	2	165619546	165624430
gene_symbol_interaction_to_bedpe <- function(df){
  print("loading hg19 genome...")
  # import hg19 genome
  ensembl <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  print("extracting gene symbols: assumes column name is `gene`")
  # extracting gene symbols
  gene_symbols <- as.character(unique(df$gene))
  print("generating gene TSS coordinates")
  # convert to gene TSS - picks transcript that has TSS closest to gene start site (largest transcript typically)
  gene_TSS <- biomaRt::select(ensembl, 
                              keys = gene_symbols,
                              columns = c('external_gene_name','chromosome_name', 'start_position', 'transcription_start_site'),
                              keytype = "external_gene_name") %>%
    dplyr::mutate(seq_diff = abs(as.numeric(start_position) - as.numeric(transcription_start_site))) %>%
    dplyr::group_by(external_gene_name) %>%
    dplyr::arrange(desc(seq_diff)) %>%
    dplyr::top_n(n=1) %>%
    dplyr::select(-c(seq_diff, start_position))
  print("adding gene TSS information")
  # add gene TSS information
  interaction_bedpe <- full_join(df, gene_TSS, by = c("gene" = "external_gene_name")) %>%
    dplyr::mutate(transcription_start_site_50bp = transcription_start_site + 50) %>%
    dplyr::mutate(chromosome_name = paste("chr",chromosome_name, sep = ""),
                  chromosome = paste("chr",chromosome, sep = ""),
                  name = row_number(),
                  score = 1,
                  strand1 = "*",
                  strand2= "*") %>%
    dplyr::select(chromosome_name, transcription_start_site, transcription_start_site_50bp, chromosome, start, end, name, score, strand1, strand2, everything())
  return(interaction_bedpe)
}


# take in Won et al interaction data and replace ENSGID for gene_name
replace_ENSGID_gene_name <- function(df){
  
  # add gene names
  ensembl <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  ENGS_ids <- as.character(df$ENSGID_for_TSS)
  Won_interactions_gene_symbols <- biomaRt::select(ensembl, 
                                                        keys = ENGS_ids,
                                                        columns = c("ensembl_gene_id","external_gene_name"),
                                                        keytype = "ensembl_gene_id")
  # merge symbols and reformat data
  Won_interactions <- full_join(df, Won_interactions_gene_symbols, by = c("ENSGID_for_TSS" = "ensembl_gene_id")) %>%
    dplyr::rename(gene = external_gene_name)
  
  return(Won_interactions)
}

# liftover
hg19_to_hg38 <- function(df, chain){
  converted_df <- liftOver(makeGRangesFromDataFrame(df, keep.extra.columns = TRUE), chain)
  converted_df <- as(unlist(converted_df), "data.frame")
  return(converted_df)
}
