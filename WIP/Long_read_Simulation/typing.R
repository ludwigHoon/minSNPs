#' \code{search_from_fastq_reads}
#'
#' @description
#' \code{search_from_fastq_reads} identify the matches
#' from a list of search strings
#' @param fastq_file fastq file containing the runs to search from
#' @param search_table a dataframe with the following columns:
#' - "type","sequence","result"
#' @param bp BiocParallel backend to use for parallelization (per file)
#' @param bp_per_read BiocParallel backend to use for parallelization (per read)
#' @param output_temp_result whether to output the temporary results
#' @param output_temp_result_dir directory to output the temporary results
#' @return will return a dataframe containing: -
#' file_path, read_id, matched_string_id, strand, match_count
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
search_from_fastq_reads <- function(fastq_file, search_tables,
    output_temp_result = TRUE, temp_result_folder = "./temp_results",
    bp = MulticoreParam(), bp_per_read = SerialParam()) {
    
    reads <- read_sequences_from_fastq(fastq_file, bp = bp)
    read_ids <- names(reads)
}