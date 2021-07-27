#' \code{read.fasta} is used to read fasta file, implementation
#' similar to seqinr, but much simpler and allow for spaces in sample name.
#' @param file file path
#' @param force_to_lower whether to transform sequences
#' to lower case, default to TRUE
#' @return Will return fastaDNA.
#' @export
read.fasta <- function(file, force_to_lower=TRUE) { # nolint
    lines <- readLines(file)
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseqc <- length(ind)
    if (nseqc == 0) {
        stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- lapply(seq_len(nseqc), function(i) {
        paste(lines[start[i]:end[i]], collapse = "")
        }
    )

    if (force_to_lower) {
        sequences <- tolower(sequences)
    }
    sequences <- as.list(sequences)
    sequences <- lapply(sequences, function(seqc) {
        return(unlist(strsplit(seqc, split = "")))
    })

    names_sequences <- lapply(seq_len(nseqc), function(i) {
        firstword <- lines[ind[i]]
        return(trimws(substr(firstword, 2, nchar(firstword))))
    })

    sequences <- lapply(seq_len(length(sequences)),
        function(i, seqs, names_seqs) {
            attr(seqs[[i]], "name") <- names_seqs[[i]]
            return(seqs[[i]])
        },
    names_seqs = names_sequences, seqs = sequences)

    names(sequences) <- names_sequences

    return(sequences)
}

#' \code{get_usual_length} is used to find out
#' the length of the sequences (W/O deletion).
#' @param seqc a list containing list of nucleotides. To keep it simple,
#' use provided read.fasta to import the fasta file.
#' @param bp is the biocparallel backend, default to serialParam,
#' most likely sufficient in most scenario
#' @import BiocParallel
#' @keywords internal
#' @return Will return the length of most of the samples,
#' the higher number is taken as tie breaker
get_usual_length <- function(seqc, bp=BiocParallel::SerialParam()) {
    sequences_length <-
        unlist(BiocParallel::bplapply(seqc, length, BPPARAM = bp))
    uniqv <- sort(unique(sequences_length), decreasing = TRUE)
    return(unlist(uniqv[which.max(tabulate(match(sequences_length, uniqv)))]))
}

#' \code{flag_allele} is used to find out a list of samples that has been
#' flagged with length different from most and will not be included in
#' computation of minimum SNPs.
#' @inheritParams get_usual_length
#' @keywords internal
#' @return Will return a list of ignored samples.
flag_allele <- function(seqc, bp=BiocParallel::SerialParam()) {
    normal_length <- get_usual_length(seqc)
    ignore_list <- BiocParallel::bplapply(seqc, function(x, norm_length) {
        return(length(x) != norm_length)
    }, norm_length = normal_length, BPPARAM = bp)
    names(which(ignore_list == TRUE))

    return(names(which(ignore_list == TRUE)))
}

#' \code{flag_position} is used to find out positions that will be ignored in
#' calculation (either not A,C,G,T or '-'), can be case sensitive
#' or insensitive.
#' @param pro_seqc Sequences after processed, i.e. all with the same length
#' @param dash_ignore whether to treat '-' as another type
#' @param accepted_char character to accept, default to c("A", "C", "T", "G")
#' @param ignore_case whether to be case insensitive, default to TRUE
#' @keywords internal
#' @return Will return a list of positions that need to be ignored.
flag_position <- function(pro_seqc, dash_ignore=TRUE,
    accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE) {

    if (dash_ignore == FALSE) {
        accepted_char <- c(accepted_char, "-")
    }

    if (ignore_case) {
        accepted_char <- tolower(accepted_char)
        pro_seqc <- lapply(pro_seqc, tolower)
    }

    ignored_position <- lapply(pro_seqc, function(iso, charset) {
        return(which(! iso %in% charset))
    }, charset = accepted_char)

    result <- sort(unique(unlist(ignored_position)))

    return(result)
}

#' \code{remove_dup_allele} is used to remove allele with the same name.
#' @inheritParams get_usual_length
#' @keywords internal
#' @return return the list of samples where only the first instance of
#' the allele with the duplicate name is kept
remove_dup_allele <- function(seqc) {
    all_sequences <- names(seqc)
    unique_sequences <- which(unique(all_sequences) %in% all_sequences)
    return(seqc[unique_sequences])
}

#' \code{process_allele} is used to returned the processed allelic profiles, by
#' removing the allele profile with duplicate name and length different
#' from most. 1st allele profile with the duplicated name is returned,
#' the longer length is taken as normal should there be 2 modes.
#' @inheritParams get_usual_length
#' @inheritParams flag_position
#' @return Will return the processed allelic profiles.
#' @export
process_allele <- function(seqc, bp=BiocParallel::SerialParam(),
    dash_ignore=TRUE, accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE) {

    processed <- list()
    seqc <- remove_dup_allele(seqc)
    ignored <- flag_allele(seqc, bp)
    processed$seqc <- seqc
    processed$ignored_allele <- ignored

    cat("Ignored samples:", "\n")
    cat(paste(ignored, collapse = ", "), "\n")

    processed$seqc <- within(processed$seqc, rm(list = ignored))

    ignored_position <- flag_position(processed$seqc,
        dash_ignore = dash_ignore, accepted_char = accepted_char,
        ignore_case = ignore_case)

    cat("Ignored ", length(ignored_position), " positions", "\n")

    processed$ignored_position <- ignored_position

    return(processed)
}