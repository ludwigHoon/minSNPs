#' \code{search_from_fastq_reads}
#'
#' @description
#' \code{search_from_fastq_reads} identify the matches
#' from a list of search strings
#' @param fastq_file fastq file containing the runs to search from
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param bp BiocParallel backend to use for parallelization
#' @param output_temp_result whether to output the temporary results
#' @param output_temp_result_dir directory to output the temporary results
#' @return will return a dataframe containing: -
#' file_path, read_id, matched_string_id, match_strand, match_count
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
search_from_fastq_reads <- function(fastq_file, search_tables,
    output_temp_result = TRUE, temp_result_folder = "./temp_results",
    bp = MulticoreParam()) {
    
    reads <- read_sequences_from_fastq(fastq_file, bp = bp)
    read_ids <- names(reads)
    names(reads) <- paste(read_ids, "_+", sep = "")
    reads_rc <- bplapply(reads, function(seq) {
        return (reverse_complement(seq))
        }, BPPARAM = bp)
    names(reads_rc) <- paste(read_ids, "_-", sep = "")
    all_reads <- c(reads, reads_rc)

    if (output_temp_result) {
        if (!dir.exists(temp_result_folder)) {
            dir.create(temp_result_folder)
        }
        temp_result_file <- file.path(temp_result_folder,
            paste0(basename(fastq_file), ".csv"))
    }

    temp_result <- bplapply(seq_len(nrow(search_tables)), function(i, search_tables, reads, output_temp, result_file){
        search_sequence <- search_tables[i, "sequence"]
        result <- sequence_reads_match_count(search_sequence, reads)
        result_df <- data.frame(reads = names(reads), sequence = search_sequence, count = unlist(result))
        if (output_temp) {
            write.table(result_df, result_file, col.names = FALSE, row.names = FALSE, append = TRUE)    
        }
        return(result_df)
    }, output_temp = output_temp_result, search_tables = search_tables, result_file = temp_result_file, reads = all_reads, BPPARAM = bp)
    
    return(do.call(rbind, temp_result))
}

sequence_reads_match_count <- function(search_sequence, reads) {
    matches <- gregexpr(paste(search_sequence, collapse = ""),
        reads)
    sequence_found_in_reads <- bplapply(matches, function(match) {
        found <- attr(match, "match.length")
        if (length(found) == 1) {
            if (found[1] == -1) {
                return(0)
            }
        }
        return(length(found))
    }, BPPARAM = BiocParallel::SerialParam())
    names(sequence_found_in_reads) <- names(reads)
    return(sequence_found_in_reads)
}
