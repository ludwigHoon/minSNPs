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
#' @importFrom data.table rbindlist
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
    }

    temp_result <- bplapply(seq_len(nrow(search_tables)), function(i, search_tables, reads, output_temp, result_folder, id){
        search_sequence <- search_tables[i, "sequence"]
        search_id <- search_tables[i, "id"]
        result <- sequence_reads_match_count(search_sequence, reads)
        result_df <- data.frame(reads = names(reads), sequence = search_sequence, count = unlist(result))
        if (output_temp) {
            result_file <- file.path(result_folder,
            paste0(id, "_", search_id, ".csv"))
            write.table(result_df, result_file, col.names = FALSE, row.names = FALSE, append = TRUE)    
        }
        return(result_df)
    }, output_temp = output_temp_result, search_tables = search_tables, result_folder = temp_result_folder,
    id = fastq_file, reads = all_reads, BPPARAM = bp)
    
    return(rbindlist(temp_result))
}

#' \code{sequence_reads_match_count}
#'
#' @description
#' \code{sequence_reads_match_count} find the number of matches
#' from a list of reads
#' @param search_sequence search string (single)
#' @param reads reads (multiple) from fastq file to search from (forward and backward)
#' @return will return a list of counts
#' @importFrom BiocParallel bplapply SerialParam
#' @export
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

scan_directory_and_search <- function() {
    
}

update_typing_result <- function() {

}

combine_and_remove_conflict <- function(fastq_file, search_table, output_temp_result_dir,
    count_measure = "n_reads", bp = MulticoreParam()) {
    temp_result_files <- list.files(output_temp_result_dir, pattern = fastq_file, full.names = TRUE)
    temp_result <- bplapply(temp_result_files, function(file) {
        temp <- read.csv(file, header = FALSE, row.names=NULL, sep= " " )
        temp <- temp[temp$V3 > 0, ]
        if (nrow(temp) == 0) {
            return(NULL)
        }
        result <- data.frame(sequence = unique(temp$V2), raw_match = sum(temp$V3), n_reads = nrow(temp),
            reads = paste(temp$V1, collapse = "; "))
        return(result)
    }, BPPARAM = bp)
    temp_result <- rbindlist(temp_result)
    t_result <- merge(temp_result, search_table, by.x = "sequence", by.y = "sequence", all.x = TRUE)
    snp_result <- remove_snp_conflict(t_result[t_result$type == "SNP",], count_measure)
    no_conflict_result <- rbindlist(list(t_result[t_result$type == "KMER",], snp_result))
    return(no_conflict_result)
}

#' @import data.table
remove_snp_conflict <- function(result, count_measure) {
    snp_info <- strsplit(result$id, split = "_")
    snp_id <- unlist(bplapply(snp_info, `[`, 1, BPPARAM = SerialParam()))
    snp_allele <- unlist(bplapply(snp_info, `[`, 2, BPPARAM = SerialParam()))
    n_result <- cbind(result, snp_id = snp_id, snp_allele = snp_allele)
    setDT(n_result)
    ### Table is ordered
    setorderv(n_result, c("snp_id", "raw_match", "n_reads"), c(1,-1, -1))
    
    ### Drop rows if the maximum count are the same
    ### 1. Identifying all the max value for each snp_id
    temp_result <- n_result[n_result[, .I[n_reads == max(get(count_measure))], by = snp_id]$V1]
    ### 2. Drop all row with duplicated snp_id
    to_drop <- temp_result$snp_id[duplicated(temp_result$snp_id)]
    temp_result <- temp_result[!temp_result$snp_id %in% to_drop,]
    return(temp_result[, -c("snp_id", "snp_allele")])
}


# "145702" "28049"  "6495"

# n_result[n_result[, .I[n_reads == max(n_reads)], by = snp_id]$V1]

# mydf[mydf[, value == max(value), .(A, B)]$V1, ]

# group[group[, .I[pt == max(pt)], by=Subject]$V1]

# [order(snp_id, -group_sum),][J(unique(snp_id)), on = "snp_id", mult = "first"]

# n_result[J(unique(snp_id)), on = "snp_id", mult = "first"]

# bp <- MulticoreParam(workers = 64, progress = TRUE)
# temp_result_files <- list.files(".", pattern = "SRR14933395_aa.fastq.*", full.names = TRUE)
#     temp_result <- bplapply(temp_result_files, function(file) {
#         temp <- read.csv(file, header = FALSE, row.names=NULL, sep= " " )
#         temp <- temp[temp$V3 > 0, ]
#         if (nrow(temp) == 0) {
#             return(NULL)
#         }
#         result <- data.frame(sequence = unique(temp$V2), raw_match = sum(temp$V3), n_reads = nrow(temp),
#             reads = paste(temp$V1, collapse = "; "))
#         return(result)
#     }, BPPARAM = bp)

# file <- temp_result_files[1]