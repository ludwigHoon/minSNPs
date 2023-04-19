snp_fragment_mapping <- list(
    
)

fragment_to_snps <- function(fragments, snp_fragment_mapping, existing_snps){
    lapply(fragments, function(frag))

}



#' \code{process_individual_read}
#'
#' @description
#' \code{process_individual_read} is used to read intermediate 
#' output files generated by minSNPs.
#' @param read_id the filepath
#' @param sequence_id the type of data contained in intermediate file
#' @param search_table the search table containing list of fragments to search
#' @return a dataframe of the matches and metadata
process_individual_read <- function(read_file, sequence_id, search_table){
    result <- data.frame(sequence_id = c(), read_id = c(), matched_string = c(),
                         snp_id = c(), extracted_snps = c(), target_pos = c(),
                         match_count = c(), match_indexes = c(), stringsAsFactors = F)

    fragments <- read_fasta(read_file)
    read_id <- strsplit(read_file, split=".fa")[[1]][1]
    if (length(fragments) > 1){
        for (i in seq_len(length(fragments))){

        }
    } else{

    }

    return(result)
}

check_process_valid <- function(){

}

check_

result <- data.frame(fragment_id = c(), matched_string = c(), snp_id = c(), extracted_snps = c(), target_pos = c(), match_count = c(), match_indexes = c(),
  stringsAsFactors = F)
for (i in seq_len(length(test))) {
  c_count <- BiocParallel::bplapply(s_strings, function(string) {
    match <- unlist(gregexpr(string, paste(test[[i]], collapse = "")))
    return(match[which(match > 0)])
  }, BPPARAM = BiocParallel::MulticoreParam(workers = 8))
  s_m <- lapply(c_count, length)
  matched_string <- s_strings[which(s_m > 0)]
  match_count <- unlist(s_m[which(s_m > 0)])
  match_indexes <- unlist(lapply(c_count[which(s_m > 0)], paste, collapse = ", "))
  snp_id <- searches[match(matched_string, searches$search_string), "fasta_position"]
  target_pos <- searches[match(matched_string, searches$search_string), "target_position"]
  ext_snps <- searches[match(matched_string, searches$search_string), "snp_string"]
  stopifnot(length(matched_string) == length(match_count))
  stopifnot(length(matched_string) == length(match_indexes))
  new_data <- data.frame(fragment_id = rep(i, length(matched_string)), snp_id = snp_id, target_pos = target_pos, extracted_snps = ext_snps, matched_string = matched_string,
    match_count = match_count, match_indexes = match_indexes, stringsAsFactors = F)
  result <- rbind(result, new_data)
}