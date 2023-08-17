#Utility functions

#' \code{reverse_complement}
#'
#' @description
#' \code{reverse_complement} returns the reverse complement of the
#' given sequence
#' @param seq the sequence to reverse complement
#' @return reverse complemented sequence
#' @export
reverse_complement <- function(seq) {
    seq <- rev(strsplit(seq, split = "")[[1]])
    result <- lapply(seq, function(x) {
        if (toupper(x) == "A") {
            return("T")
        }
        if (toupper(x) == "C") {
            return("G")
        }
        if (toupper(x) == "T") {
            return("A")
        }
        if (toupper(x) == "G") {
            return("C")
        }
        return(x)
    })
    return(paste(result, collapse = ""))
}

#' \code{match_count}
#'
#' @description
#' \code{match_count} return the number of matches of the target string in the
#' given sequence
#' @param target the search target
#' @param search_from the sequence to search from
#' @return number of matches
#' @export
match_count <- function(target, search_from) {
    matches <- gregexpr(paste(target, collapse = ""),
        paste(search_from, collapse = ""))
    found <- attr(matches[[1]], "match.length")
    if (length(found) == 1) {
        if (found[1] == -1) {
           return(0)
        }
    }
    return(length(found))
}

#' \code{read_sequences_from_fastq}
#'
#' @description
#' \code{read_sequences_from_fastq} get the sequences
#' from a fastq file, it completely ignores the quality scores
#' @param fastq_file location of the fastq file
#' @param force_to_upper whether to transform sequences
#' to upper case, default to TRUE
#' @param max_n_reads maximum number of reads to read, default to -1 (all)
#' @param skip_n_reads number of reads to skip, default to 0
#' @param output_quality whether to output the quality scores, default to TRUE
#' @param bp BiocParallel backend to use for parallelization
#' @param quality_offset the quality offset to use, default to 33
#' @return will return a list of sequences, with qualities as attribute
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
read_sequences_from_fastq <- function(fastq_file, force_to_upper = TRUE, skip_n_reads = 0,
    max_n_reads = -1, output_quality = TRUE, quality_offset = 33, bp = MulticoreParam()) {

    if (max_n_reads != -1 && skip_n_reads != 0) {
        n_lines <- (max_n_reads * 4) + (skip_n_reads * 4)
    } else {
        n_lines <- (max_n_reads * 4)
    }

    lines <- readLines(fastq_file, n = n_lines)
    if (skip_n_reads != 0){
        if ((4 * skip_n_reads + 1) > length(lines)){
            warning("skip_n_reads is too large: No reads beyond that")
            return(NULL)
        }
        lines <- lines[seq(4 * skip_n_reads + 1, length(lines), 1)]
    }
    seqs_id <- lines[seq(1, length(lines), 4)]
    seqs <- lines[seq(2, length(lines), 4)]
    

    seqs_id <- unlist(bplapply(seqs_id, function(id) {
        id <- strsplit(id, split = "")[[1]]
        return(paste0(id[2:length(id)], collapse = ""))
    }, BPPARAM = bp))

    if (force_to_upper) {
        seqs <- unlist(bplapply(seqs, function(seq) {
            return(toupper(seq))
        }, BPPARAM = bp))
    }

    seqs <- as.list(seqs)
    names(seqs) <- seqs_id
    if (output_quality) {
        qualities_data <- lines[seq(4, length(lines), 4)]
        qualities <- bplapply(qualities_data, function(x) {
            return(utf8ToInt(x) - quality_offset)
        }, BPPARAM = bp)
        attr(seqs, "qualities") <- qualities
    }
    return(seqs)
}

#' \code{get_snps_set}
#'
#' @description
#' \code{get_snps_set} extract the SNP sets from the output of
#' `find_optimised_snps`.
#' @param results output from `find_optimised_snps`
#' @param as output format, either `data.frame` or `list`.
#' @importFrom utils tail
#' @return will return either 
#' 1. a dataframe containing SNPs_set (SNP position separated by ",") and score
#' 2. a list containing SNPs_set (SNP position as numeric vector) and score (attr of the list)
#' @export
get_snps_set <- function(results, as = "data.frame") {
    if (!as %in% c("data.frame", "list"))
        stop("as must be either data.frame or list")

    snp_sets <- sapply(results$results, function(x){
        return(
            tail(names(x), n = 1)
        )
    })
    score <- sapply(results$results, function(x){
        return(
            tail(x, n = 1)[[1]]
        )
    })
    if (as == "data.frame"){
        result <- data.frame(SNPs_set = snp_sets, score = score)
    } else {
        result <- unname(lapply(strsplit(snp_sets, split = ", "), as.numeric))
        attr(result, "score") <- unname(score)
    }
    return(result)
}

#' \code{process_result_file}
#'
#' @description
#' \code{process_result_file} extract the SNP sets from the saved output file.
#' @param result_filepath is the path of the saved output file.
#' @return will return a list containing SNPs_set (SNP position as numeric vector).
#' @export
process_result_file <- function(result_filepath) {
    f <- file(result_filepath, "r")
    preparse <- readLines(f)
    result <- list()
    cur <- 1
    is_result <- FALSE
    for (i in seq_len(length(preparse))) {
        if (length(grep("^Result", preparse[[i]])) >= 1) {
            is_result <- TRUE
        }
        if (is_result) {
            if (gsub("\\s+", "", preparse[[i + 1]]) == "") {
                result[[cur]] <- as.numeric(
                    strsplit(
                        gsub(
                            "\"", "",
                            strsplit(preparse[[i]], split = "\t")[[1]][1]
                        ),
                        split = ", "
                    )[[1]]
                )
                cur <- cur + 1
                is_result <- FALSE
            }
        }
    }
    close(f)
    final_selected <- (result)
    return(final_selected)
}