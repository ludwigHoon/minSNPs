#' \code{usual_length} is used to find out
#' the length of the sequences (W/O deletion).
#' @param seq a list containing list of nucleotides. To keep it simple,
#' use provided read.fasta to import the fasta file.
#' @return Will return the length of most of the samples.
#' @example tests/integrated.R
#' @export
usual_length <- function(seq) {
    sequences_length <- unlist(lapply(seq, length))
    uniqv <- unique(sequences_length)
    return(unlist(uniqv[which.max(tabulate(match(seq, uniqv)))]))
}

#' \code{flag_allele} is used to find out a list of samples
#' that has been flagged and will not be included in computation
#' of minimum SNPs.
#' @inheritParams usualLength
#' @return Will return a list of ignored samples.
#' @export
flag_allele <- function(seq) {
    normal_length <- usual_length(seq)
    ignore_list <- lapply(seq, function(x, norm_length) {
        return(length(x) != norm_length)
    }, norm_length = normal_length)
    names(which(ignore_list == TRUE))

    return(names(which(ignore_list == TRUE)))
}

#' \code{process_allele} is used to returned the processed allelic
#' profiles.
#' @inheritParams usualLength
#' @return Will return the processed allelic profiles.
#' @export
process_allele <- function(seq) {
    ignored <- flag_allele(seq)
    processed <- seq
    cat("Ignored samples:", "\n")
    for (a in ignored) {
        cat(a, "\n")
        processed[[a]] <- NULL
    }
    return(processed)
}

#' \code{flag_position} is used to find out positions that will be ignored in
#' calculation (either not A,C,G,T or '-'), can be case sensitive or insensitive
#' @param pro_seq Sequences after processed, i.e. all with the same length
#' @param dash_ignore whether to treat '-' as another type
#' @param accepted_char character to accept, default to c("A", "C", "T", "G")
#' @param ignore_case whether to be case insensitive, default to TRUE
#' @return Will return a list of positions that need to be ignored.
#' @export
flag_position <- function(pro_seq, dash_ignore=TRUE,
    accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE) {

    if (dash_ignore == FALSE) {
        accepted_char <- c(accepted_char, "-")
    }

    if (ignore_case) {
        accepted_char <- tolower(accepted_char)
        seq <- lapply(seq, tolower)
    }

    ignored_position <- lapply(seq, function(iso, charset) {
        return(which(! iso %in% charset))
    }, charset = accepted_char)

    result <- unique(unlist(ignored_position))

    cat("Ignored ", length(result), " positions", "\n")

    return(result)
}

#' \code{read.fasta} is used to read fasta file, implementation
#' similar to seqinr, but much simpler and allow for spaces in sample name
#' @param file file path
#' @param force_to_lower whether to transform sequences
#' to lower case, default to TRUE
#' @return Will return fastaDNA.
#' @export
read.fasta <- function(file, force_to_lower=TRUE) { # nolint
    lines <- readLines(file)
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseq <- length(ind)
    if (nseq == 0) {
        stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- lapply(seq_len(nseq), function(i) {
        paste(lines[start[i]:end[i]], collapse = "")
        }
    )

    if (force_to_lower) {
        sequences <- tolower(sequences)
    }
    sequences <- as.list(sequences)
    sequences <- lapply(sequences, function(seq) {
        return(unlist(strsplit(seq, split = "")))
    })

    names_sequences <- lapply(seq_len(nseq), function(i) {
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