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
    found <- attr(match[[1]], "match.length")
    if (length(found) == 1) {
        if (found[1] == -1) {
            found <- 0
        }
    } else{
        found <- length(found)
    }
    return(found)
}

#' \code{read_sequences_from_fastq}
#'
#' @description
#' \code{read_sequences_from_fastq} get the sequences
#' from a fastq file, it completely ignores the quality scores
#' @param fastq_file location of the fastq file
#' @param force_to_upper whether to transform sequences
#' to upper case, default to TRUE
#' @param bp BiocParallel backend to use for parallelization
#' @return will return a list of sequences
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
read_sequences_from_fastq <- function(fastq_file, force_to_upper = TRUE,
    bp = MulticoreParam()) {
    lines <- readLines(fastq_file)
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
    return(seqs)
}