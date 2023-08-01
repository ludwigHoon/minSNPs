############# INDIVIDUAL SEARCHES #############

#' \code{sequence_reads_match_count}
#' 
#' @description
#' \code{sequence_reads_match_count} look for the search sequences in reads and return the matches indexes and mean qualities
#' @param search_sequence the search sequence to look for where `.` stands for any character.
#' @param reads the sequences reads to search for. 
#' @param qualities the qualities of each bases in the reads.
#' @return will return a list containing for each read: -
#' count, mean_quality, indexes
#' @export
sequence_reads_match_count <- function(search_sequence, reads, qualities) {
    matches <- gregexpr(paste(search_sequence, collapse = ""),
        reads)
    names(matches) <- names(reads)
    sequence_found_in_reads <- lapply(names(matches), function(match_read, matches, qualities) {
        match <- matches[[match_read]]
        if (match[1] == -1) {
                return(list(count=0, mean_quality=NA, indexes=NA))
        }
        count <- length(match)
        mean_quality <- rep(0, length(match))
        indexes <- as.numeric(match)
        for (i in seq_len(length(match))) {
            ind <- match[i]
            range <- ind:(attr(match, "match.length")[i] + ind - 1)
            mean_quality[i] <- round(mean(qualities[[match_read]][range]), 2)
        }
        return(list(count=count, mean_quality=mean_quality, indexes=indexes))
    }, matches = matches, qualities = qualities)
    names(sequence_found_in_reads) <- names(reads)
    return(sequence_found_in_reads)
}

#' \code{search_from_fastq_reads}
#'
#' @description
#' \code{search_from_fastq_reads} identify the matches
#' from a list of search strings
#' @param fastq_file fastq file containing the runs to search from
#' @param search_tables a dataframe with the following columns:
#' - ["id"],"type",["sequence"],"strand","result","extra","match_ref_seq"
#' @param bp BiocParallel backend to use for parallelization
#' @param skip_n_reads number of reads to skip, default is 0
#' @param max_n_reads maximum number of reads to read, default to -1 (all)
#' @param quality_offset the quality offset to use, default to 33
#' @param progress whether to show the progress bar
#' @param output_read_length whether to output the read length, NULL - do not output; csv - output to csv file; data - output to result
#' @param output_temp_result whether to output the temporary results
#' @param temp_result_folder directory to output the temporary results
#' @param simplify_id simplify and shorten the read id to the first part
#' @return will return a list of dataframe containing: -
#' `search_id`, `sequence`, `reads`, `raw_match`, `mean_qualities`, `indexes`.
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
search_from_fastq_reads <- function(fastq_file, search_tables, skip_n_reads = 0, progress = TRUE, max_n_reads = -1,
    quality_offset = 33, output_temp_result = TRUE, temp_result_folder = "./temp_results",
    simplify_id = TRUE, output_read_length = TRUE, bp = MulticoreParam()) {
    
    reads <- read_sequences_from_fastq(fastq_file, quality_offset = quality_offset, max_n_reads = max_n_reads, skip_n_reads = skip_n_reads, bp = bp)

    read_ids <- names(reads)
    if (simplify_id){
        read_ids <- sapply(strsplit(read_ids, split = " "), `[`, 1)
    }

    qualities <- attr(reads, "qualities")
    if (is.null(qualities)) {
        qualities <- lapply(reads, function(seq) {
            return(rep(quality_offset, nchar(seq)))
        })
    }
    if (output_temp_result) {
        if (!dir.exists(temp_result_folder)) {
            dir.create(temp_result_folder)
        }
    }

    read_length_data <- NULL
    if (output_read_length) {
        reads_length <- nchar(reads)
        read_length_data <- data.frame(reads_id = read_ids, reads_length = as.numeric(reads_length))
        
        if (output_temp_result){
            write.csv(read_length_data,
                paste0(temp_result_folder, "/", basename(fastq_file), "_read_lengths.csv"),
                row.names = FALSE)
        }
    }

    names(reads) <- paste(read_ids, "_+", sep = "")
    reads_rc <- bplapply(reads, function(seq) {
        return(reverse_complement(seq))
    }, BPPARAM = bp)
    qualities_rc <- bplapply(qualities, function(seq) {
        return(rev(seq))
    }, BPPARAM = bp)
 
    names(reads_rc) <- paste(read_ids, "_-", sep = "")
    all_reads <- c(reads, reads_rc)
    all_qualities <- c(qualities, qualities_rc)
    names(all_qualities) <- names(all_reads)
    old_pr <- bp$progressbar
    if (progress) {
        bp$progressbar <- progress
    }
    print("SEARCHING ...")
    temp_result <- bplapply(seq_len(nrow(search_tables)), function(i, search_tables, reads, qualities, output_temp, result_folder, id){
        search_sequence <- search_tables[i, "sequence"]
        search_id <- search_tables[i, "id"]
        result <- sequence_reads_match_count(search_sequence, reads, qualities)
        counts <- sapply(result, function(x) {
            return(x$count)
        })
        mean_qualities <- sapply(result, function(x) {
            return(paste(x$mean_quality, collapse = ","))
        })
        indexes <- sapply(result, function(x){
            return(paste(x$indexes, collapse = ","))
        })
        result_df <- data.frame(
            search_id = search_id,
            sequence = search_sequence,
            reads = names(reads),
            raw_match = counts,
            mean_qualities = mean_qualities,
            indexes = indexes
        )
        if (output_temp) {
            result_file <- file.path(result_folder,
            paste0(id, "_", search_id, ".csv"))
            write.csv(result_df, result_file, sep = "\t", row.names = FALSE)
        }
        return(result_df)
    }, output_temp = output_temp_result, search_tables = search_tables, result_folder = temp_result_folder,
    id = basename(fastq_file), reads = all_reads, qualities = all_qualities, BPPARAM = bp)
    if (progress) {
        bp$progressbar <- old_pr
    }
    names(temp_result) <- search_tables$sequence
    
    return(
        list(result = temp_result,
        read_length = read_length_data)
    )
}

######### COMBINED FUNCTIONS #########

#' \code{combine_search_string_result}
#'
#' @description
#' \code{combine_search_string_result} combines the search results from \code{search_from_fastq_reads}
#' @param results the dataframes to collapse.
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param bp BiocParallel backend to use for parallelization
#' @param append_to_current_result the dataframe of previous result to append to
#' @return will return a dataframe containing: -
#' `sequence`, `search_id`, `reads`, `raw_match`, `mean_qualities`, `indexes`, `id`, `type`,
#' `strand`, `result`, `extra`, `match_ref_seq`, `n_reads`
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist
#' @export
combine_search_string_result <- function(results, search_table,
    append_to_current_result = data.frame(), bp = MulticoreParam()) {
    
    t_result <- merge(results, search_table, by.x = "sequence", by.y = "sequence", all.x = TRUE)
    t_result <- rbindlist(list(t_result, append_to_current_result))

    if (nrow(t_result) == length(unique(t_result$sequence))) {
        return(t_result)
    } else{
        all_sequences_count <- table(t_result$sequence)
        to_retain <- names(which(all_sequences_count == 1))
        result <- t_result[t_result$sequence %in% to_retain, ]
        
        to_combine <- names(which(all_sequences_count > 1))
        result_to_combine <- bplapply(to_combine, function(seq, temp_result){
            temps <- temp_result[temp_result$sequence == seq, ]
            res <- temps[1, ]
            res[, "raw_match"] <- sum(temps$raw_match)
            res[, "n_reads"] <- sum(temps$n_reads)

            reads <- unlist(strsplit(temps$reads, split = ";"))
            mean_qualities <- unlist(strsplit(temps$mean_qualities, split = ";"))
            indexes <- unlist(strsplit(temps$indexes, split = ";"))

            res[, "reads"] <- paste(reads, collapse = ";")
            res[, "mean_qualities"] <- paste(mean_qualities, collapse = ";")
            res[, "indexes"] <- paste(indexes, collapse = ";")
            return(res)
        }, temp_result = t_result, BPPARAM = bp)
        result_to_combine[[length(result_to_combine) + 1]] <- result
        fin_result <- rbindlist(result_to_combine)
        return(fin_result)   
    }
}

#' \code{combine_search_string_result_from_list}
#'
#' @description
#' \code{combine_search_string_result_from_list} combines the search results from \code{search_from_fastq_reads}
#' @param results the dataframes from \code{search_from_fastq_reads} to combine.
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param bp BiocParallel backend to use for parallelization
#' @param append_to_current_result the dataframe of previous result to append to
#' @return will return a dataframe containing: -
#' `sequence`, `search_id`, `reads`, `raw_match`, `mean_qualities`, `indexes`, `id`, `type`,
#' `strand`, `result`, `extra`, `match_ref_seq`, `n_reads`
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist
#' @export
combine_search_string_result_from_list <- function(results, search_table,
    append_to_current_result = data.frame(), bp = MulticoreParam()) {
    
    temp_result <- bplapply(results, function(temp) {
        temp <- temp[temp$raw_match > 0, ]
        if (nrow(temp) == 0) {
            return(
                setNames(data.frame(matrix(ncol = 7, nrow = 0)),
                    c("search_id", "sequence", "raw_match", "n_reads", "reads", "mean_qualities", "indexes"))
            )
        }
        result <- data.frame(
            search_id = unique(temp$search_id),
            sequence = unique(temp$sequence),
            raw_match = sum(temp$raw_match),
            n_reads = nrow(temp),
            reads = paste(temp$reads, collapse = ";"),
            mean_qualities = paste(temp$mean_qualities, collapse=";"),
            indexes = paste(temp$indexes, collapse = ";"))
        return(result)
    }, BPPARAM = bp)
    
    return(
        combine_search_string_result(
            results = rbindlist(temp_result),
            search_table = search_table,
            append_to_current_result = append_to_current_result,
            bp  = bp
        )
    )

}      


#' \code{combine_search_string_result_from_files}
#'
#' \code{combine_search_string_result_from_files}
#' \code{combine_search_string_result} combines the search results from temp file generated from \code{search_from_fastq_reads}
#' @param result_files the output files from \code{search_from_fastq_reads} to combine
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param bp BiocParallel backend to use for parallelization
#' @param append_to_current_result the dataframe of result to append to
#' @return will return a dataframe containing: -
#' `sequence`, `search_id`, `reads`, `raw_match`, `mean_qualities`, `indexes`, `id`, `type`,
#' `strand`, `result`, `extra`, `match_ref_seq`, `n_reads`
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist
#' @export
combine_search_string_result_from_files <- function(result_files, search_table,
    append_to_current_result = data.frame(), bp = MulticoreParam()) {
    temp_result <- bplapply(result_files, function(file) {
        temp <- read.csv(file, header = TRUE, row.names = NULL, sep = ",")
        return(temp)
    }, BPPARAM = bp)

    return(
        combine_search_string_result_from_list(
            results = temp_result,
            search_table = search_table,
            append_to_current_result = append_to_current_result,
            bp  = bp
        )
    )
}

#' \code{infer_from_combined}
#' @description
#' \code{infer_from_combined} infers the results (presence/absense of genes & CC) from the combined result
#' @param combined_result the combined result from \code{combine_search_string_result}
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @return a dataframe containing the following columns:
#' - type, rank, result, reads_count, proportion_matched, pass_filter
infer_from_combined <- function(combined_result, search_table, ...) {
    result <- list()
    arguments <- list(...)
    analysis_types <- unique(combined_result$type)
    arguments[["search_table"]] <- search_table

    for (type in analysis_types){
        process_method <- get_all_process_methods(type)
        if (is.null(process_method)){
            warning(paste0("No process method for type: ", type))
            next
        }
        arguments[["partial_result"]] <- combined_result[combined_result$type == "KMER", ]
        result[[type]] <- do.call(process_method, args = arguments)
    }
    result_df <- rbindlist(result)
    
    return(list(
        result = result_df,
        snps_found = p_snp_result$snps_found,
        proportion_snps_found = p_snp_result$proportion_snps_found
        )
    )
}

#' \code{get_all_process_methods}
#'
#' @description
#' \code{get_all_process_methods} is used to get the metrics function
#' and required parameters. Additional metric may be set by
#' assigning it to `MinSNPs_process_methods` variable.
#' @param process_name name of the metric, "" to return all,
#' `SNP` or `KMER` are provided as default.
#' @return a list, including the function to process the search sequence result
#' @export
get_all_process_methods <- function(process_name = ""){
    if (! exists("MinSNPs_process_methods")) {
        MinSNPs_process_methods <- list( #nolint
            "SNP" = process_snp_result,
            "KMER" = process_kmer_result,
        )
    }
    if (process_name == "") {
        return(MinSNPs_process_methods)
    }
    return(MinSNPs_process_methods[[process_name]])
}

#' \code{process_kmer_result}
#' 
#' @description
#' \code{process_kmer_result} processes the KMER result from \code{infer_from_combined}
#' @param partial_result the result from \code{infer_from_combined} with only KMER
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param min_match_per_read the minimum number of kmer matches in a read, discarding reads with less than this number
#' @param ... ignored
#' @return a dataframe containing the following columns:
#' - type, rank, result, reads_count, proportion_matched, pass_filter
#' @importFrom data.table rbindlist
#' @export
process_kmer_result <- function(partial_result, search_table, min_match_per_read = 1, ...) {
    partial_result <- kmer_only
    result <- list()
    for (gene in unique(search_table[search_table$type == "KMER","result"])) {
        total_searched <- nrow(search_table[search_table$type == "KMER" & search_table$result == gene,])
        
        # Discard reads with less than min_match_per_read
        all_reads <- unlist(strsplit(kmer_only[kmer_only$result == gene, "reads"], split = ";|,"))
        n_match_in_read <- table(kmer_only[kmer_only$result == gene, "reads"])
        filtered_reads <- names(n_match_in_read[n_match_in_read >= min_match_per_read])
        filtered_index <- apply(sapply(filtered_reads, grepl, x= kmer_only$reads), 1, sum) == length(filtered_reads)
        kmer_only2 <- kmer_only[kmer_only$result == gene & filtered_index, ]
        
        prop <- length(unique(kmer_only2[kmer_only2$result == gene, "sequence"])) / total_searched
        r_count <- length(unique(unlist(strsplit(kmer_only2[kmer_only2$result == gene, "reads"], split = ";|,"))))
        result[[gene]] <- data.frame(type = "KMER", rank = 1, result = gene, reads_count = r_count,
            proportion_matched = prop, pass_filter = prop >= 0.8)
    }
    result_df <- rbindlist(result)
    return(result_df)
}


#' \code{remove_snp_conflic}
#' 
#' @description
#' \code{remove_snp_conflic} removes the reads with SNPs conflicts from the result
#' @param result the result from \code{infer_from_combined}
#' @param count_measure the column name of the count measure to use for removing the conflicts
#' @return a dataframe containing the same columns as the input result with row containing conflicts removed
remove_snp_conflict <- function(result, count_measure = "n_reads") {
    snp_info <- strsplit(result$id, split = "_")
    snp_id <- sapply(snp_info, `[`, 1)
    snp_allele <- sapply(snp_info, `[`, 2)
    n_result <- cbind(result, snp_id = snp_id, snp_allele = snp_allele)
    data.table::setDT(n_result)
    ### Table is ordered
    data.table::setorderv(n_result, c("snp_id", count_measure), c(1, -1))
    
    ### Drop rows if the maximum count are the same
    ### 1. Identifying all the max value for each snp_id
    temp_result <- n_result[n_result[, .I[n_reads == max(get(count_measure))], by = snp_id]$V1]
    ### 2. Drop all row with duplicated snp_id
    to_drop <- temp_result$snp_id[duplicated(temp_result$snp_id)]
    temp_result <- temp_result[!temp_result$snp_id %in% to_drop,]
    return(temp_result[, -c("snp_id", "snp_allele")])
}

#' \code{process_snp_result}
#' 
#' @description
#' \code{process_snp_result} processes the SNP result from \code{infer_from_combined}
#' @param partial_result the result from \code{infer_from_combined} with only SNP
#' @param search_table a dataframe with the following columns:
#' - "id","type","sequence","strand","result","extra","match_ref_seq"
#' @param count_measure the column name of the count measure to use for removing the conflicts
#' @param ... ignored
#' @return a list containing:
#' - result: a dataframe containing the following columns:
#'  - type, rank, result, reads_count, proportion_matched, pass_filter
#' - snps_found: a vector containing the SNPs ID that have been identified without conflict
#' - proportion_snps_found: the proportion of SNPs found without conflict
#' @export
process_snp_result <- function(partial_result, search_table, count_measure = "n_reads", ...) {
    snp_only <- partial_result
    snp_only <- remove_snp_conflict(snp_only, count_measure = count_measure)

    snp_info <- strsplit(snp_only$id, split = "_")
    snp_id <- sapply(snp_info, `[`, 1)
    
    searched_snps <- search_table[search_table$type == "SNP", "id"]
    searched_snps_id <- sapply(strsplit(searched_snps, split = "_"), `[`, 1)

    cc_result <- table(unlist(sapply(snp_only[, "result"], strsplit, split = ";|,")))
    cc_result <- cc_result[order(cc_result, decreasing = TRUE)][1:10]
    
    r_count <- sapply(names(cc_result), function(cc){
        temp <- length(
            which(grepl(cc, unlist(snp_only$result))
                == TRUE)
        )
        return(temp)
    })

    stopifnot(names(r_count) == names(cc_result))
    setDF(snp_only)
    prop <- sapply(names(cc_result), function(cc){
        temp <- length(
            unique(snp_only[grepl(cc, unlist(snp_only$result)),"sequence"])    
        ) / length(
            which(grepl(cc, unlist(search_table$result))
                == TRUE)
        )
        return(temp)
    })
    stopifnot(names(prop) == names(cc_result))

    read_count <- sapply(names(cc_result), function(cc){
        sum_n_reads <- sum(snp_only[grepl(cc, unlist(snp_only$result)),"n_reads"])
        return(sum_n_reads)
    })

    result_df <- data.frame(type = rep("SNP", 10), rank = c(1:10), result = names(cc_result), reads_count = unname(r_count),
        proportion_matched = unname(prop), 
        pass_filter = (
            prop >= 0.8 &
            read_count >= 10 * prop * length(unique(snp_id)) # average read-depth of the reads with respective SNPs = 10
        )
    )
    
    return(list(
        result = result_df,
        snps_found = snp_id,
        proportion_snps_found = length(unique(snp_id)) / length(unique(searched_snps_id)) 
        )
    )
}