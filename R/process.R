#' \code{read_fasta}
#'
#' @description
#' \code{read_fasta} is used to read fasta file, implementation
#' similar to seqinr, but much simpler and allow for spaces in sample name.
#' @param file file path
#' @param force_to_upper whether to transform sequences
#' to upper case, default to TRUE
#' @param bp is the biocparallel backend, default to serialParam,
#' most likely sufficient in most scenario
#' @return Will return list of named character vectors.
#' @export
read_fasta <- function(file, force_to_upper=TRUE, bp = SerialParam()) {
    lines <- readLines(file)
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseqc <- length(ind)
    if (nseqc == 0) {
        stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- bplapply(seq_len(nseqc), function(i, force_to_upper) {
        seq <- unlist(strsplit(lines[start[i]:end[i]], split = ""))
        if (force_to_upper) {
            seq <- toupper(seq)
        }
        return(seq)
    }, force_to_upper = force_to_upper, BPPARAM = bp)

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

#' \code{write_fasta}
#'
#' @description
#' \code{write_fasta} is used to write
#' the named character vectors to fasta file.
#' @inheritParams get_usual_length
#' @param filename filename of the output file
#' @return will write the alignments to file
#' @export
write_fasta <- function(seqc, filename) {
    outfile <- file(description = filename, open = "w")

    all_seq_names <- names(seqc)

    write_one_seq <- function(name, seq) {
        writeLines(paste(">", name, sep = ""), outfile)
        this_seq <- as.vector(as.list(seqc)[[name]])
        l <- length(this_seq)
        q <- floor(l / 60)
        r <- l - (60 * q)
        if (q > 0) {
            sapply(seq_len(q), function(x) {
                writeLines(paste(this_seq[(60 * (x - 1) + 1):(60 * x)],
                    collapse = ""), outfile)
            })
        }
        if (r > 0) {
            writeLines(paste(this_seq[(60 * q + 1):l], collapse = ""), outfile)
        }
    }

    lapply(all_seq_names, write_one_seq, seq = seqc)

    close(outfile)
}

#' \code{get_usual_length}
#'
#' @description
#' \code{get_usual_length} is used to find out
#' the length of the sequences (W/O deletion).
#' @param seqc a list containing list of nucleotides. To keep it simple,
#' use provided read_fasta to import the fasta file.
#' @param bp is the biocparallel backend, default to serialParam,
#' most likely sufficient in most scenario
#' @importFrom BiocParallel bplapply SerialParam
#' @keywords internal
#' @return Will return the length of most of the samples,
#' the higher number is taken as tie breaker
get_usual_length <- function(seqc, bp=SerialParam()) {
    sequences_length <-
        unlist(bplapply(seqc, length, BPPARAM = bp))
    uniqv <- sort(unique(sequences_length), decreasing = TRUE)
    return(unlist(uniqv[which.max(tabulate(match(sequences_length, uniqv)))]))
}

#' \code{flag_allele}
#'
#' @description
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

#' \code{flag_position}
#'
#' @description
#' \code{flag_position} is used to find out positions that will be ignored in
#' calculation (either not A,C,G,T or '-'), can be case sensitive
#' or insensitive.
#' @param pro_seqc Sequences after processed, i.e. all with the same length
#' @param dash_ignore whether to treat '-' as another type
#' @param accepted_char character to accept, default to c("A", "C", "T", "G")
#' @param ignore_case whether to be case insensitive, default to TRUE
#' @param remove_invariant whether to remove invariant positions, default to FALSE
#' @param biallelic_only whether to remove positions with more than 2 alleles, default to FALSE
#' @keywords internal
#' @return Will return a list of positions that need to be ignored.
flag_position <- function(pro_seqc, dash_ignore=TRUE,
    accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE, remove_invariant=FALSE, biallelic_only=FALSE, bp = SerialParam()) {

    if (dash_ignore == FALSE) {
        accepted_char <- c(accepted_char, "-")
    }

    if (ignore_case) {
        accepted_char <- toupper(accepted_char)
        pro_seqc <- lapply(pro_seqc, function(iso) {
            toupper(as.vector(iso))
        })
    }

    seq <- as.list(pro_seqc)
    ignored_position <- bplapply(seq, function(iso, charset) {
        return(which(! as.vector(iso) %in% charset))
    }, charset = accepted_char, BPPARAM = bp)

    more_pos <- bplapply(
        seq_len(length(seq[[1]])),
            function(pos, remove_invariant, biallelic_only) {
            l_pat <- length(unique(generate_pattern(seq, pos)))
            if (remove_invariant & biallelic_only & (l_pat != 2)) {
                return(pos)
            }
            if (remove_invariant & l_pat < 2) {
                return(pos)
            }
            if (biallelic_only & l_pat > 2) {
                return(pos)
            } else {
                return(NULL)
            }
        }, remove_invariant = remove_invariant, biallelic_only = biallelic_only,
    BPPARAM = bp)

    ignored_position <- c(ignored_position, unlist(more_pos))
    result <- sort(unique(unlist(ignored_position)))

    return(result)
}

#' \code{remove_dup_isolate}
#'
#' @description
#' \code{remove_dup_isolate} is used to remove isolate with the same name.
#' @inheritParams get_usual_length
#' @keywords internal
#' @return return the list of samples where only the first instance of
#' the isolate with the duplicate name is kept
remove_dup_isolate <- function(seqc) {
    all_sequences <- names(seqc)
    unique_sequences <- unique(all_sequences)

    if (any(duplicated(all_sequences))) {
        dups <- unique(all_sequences[duplicated(all_sequences)])
        cat("Found multiple isolates with same names, only first is taken:\n",
            paste(dups, collapse = ", "), "\n\n", sep = "")
    }
    return(seqc[unique_sequences])
}


#' \code{process_allele}
#'
#' @description
#' \code{process_allele} is used to returned the processed allelic profiles, by
#' removing the allele profile with duplicate name and length different
#' from most. 1st allele profile with the duplicated name is returned,
#' the longer length is taken as normal should there be 2 modes.
#' @param check_length Check the length of each sample in the matrix, default to TRUE
#' @param check_bases Check the bases of each sample in the matrix, default to TRUE
#' @inheritParams get_usual_length
#' @inheritParams flag_position
#' @return Will return the processed allelic profiles.
#' @export
process_allele <- function(seqc, bp=BiocParallel::SerialParam(), check_length=TRUE, check_bases=TRUE,
    dash_ignore=TRUE, accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE,
    remove_invariant=FALSE, biallelic_only=FALSE) {

    processed <- list()
    seqc <- remove_dup_isolate(seqc)
    
    ignored <- NULL
    if (check_length){
        ignored <- flag_allele(seqc, bp)
    }
       
    processed$seqc <- as.list(seqc)
    processed$ignored_allele <- ignored

    cat("Ignored samples:", "\n")
    cat(paste(ignored, collapse = ", "), "\n")

    processed$seqc <- within(processed$seqc, rm(list = ignored))
    names_seqs <-  names(processed$seqc)

    ignored_position <- c()
    if (check_bases) {
        ignored_position <- flag_position(processed$seqc,
            dash_ignore = dash_ignore, accepted_char = accepted_char,
            ignore_case = ignore_case, bp = bp)    
    }
    
    cat("Ignored ", length(ignored_position), " positions", "\n")

    processed$ignored_position <- ignored_position
    processed$check_length <- check_length
    processed$length <- length(processed$seqc[[1]])
    class(processed) <- "processed_seqs"
    return(processed)
}

translation <- list(
    M = c("A", "C"),
    R = c("A", "G"),
    W = c("A", "T"),
    S = c("C", "G"),
    Y = c("C", "T"),
    K = c("G", "T"),
    V = c("A", "C", "G"),
    H = c("A", "C", "T"),
    D = c("A", "G", "T"),
    B = c("C", "G", "T")
)

#' \code{resolve_IUPAC_missing}
#'
#' @description
#' \code{resolve_IUPAC_missing} is used to replace the
#' ambiguity codes found in the sequences.
#' @param seqc the sequences to be processed
#' @param log_operation whether to log the operation
#' @param log_file log file to write the operations
#' @param max_ambiguity proportion of ambiguity codes to tolerate,
#' -1 = ignore. Default to -1
#' @param replace_method how to substitute the ambiguity codes,
#' current supported methods:random and most_common, default to "random".
#' @param N_is_any_base whether to treat N as any base or substitute it
#' with one of the alleles found at the position.
#' @param output_progress whether to output progress
#' @param bp the BiocParallel backend
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom BiocParallel bplapply MulticoreParam
#' @import data.table
#' @return Will return the processed sequences.
#' @export
    resolve_IUPAC_missing <- function(seqc, log_operation = TRUE, # nolint
        log_file = "replace.log", max_ambiguity = -1,
        replace_method = "random", N_is_any_base = FALSE, output_progress = TRUE, bp = MulticoreParam()) { #nolint

        accepted_char <- c("A", "C", "T", "G")

        # Whether N should be randomly resolved to any of the accepted bases
        # or if it should be 1 of the bases present in any of the isolates
        if (N_is_any_base) {
            translation[["N"]] <- accepted_char
        }

        # Logging all operation
        if (log_operation) {
            cat("Position\tIsolate\tOriginal\tReplaced\n",
                file = log_file, append = FALSE)
        }

        if ((max_ambiguity == -1 || max_ambiguity == 1) &&
            (replace_method == "most_common" || replace_method == "most_common_parallel")) {
            stop("most_common must have max_ambiguity < 1")
        }

        if (replace_method == "most_common") {
            if (output_progress) {
                pb <- txtProgressBar(min = 0, max = length(seqc[[1]]),
                    initial = 0, style = 3)
            }
            for (p in seq_len(length(seqc[[1]]))) {
                # All the nucleotides at this position
                nucleotides <- lapply(seqc, `[[`, p)
                # These need to be modified
                to_modify <- which(!nucleotides %in% accepted_char)

                # These are valid
                valids <- which(nucleotides %in% accepted_char)
                # Most common base
                valid_bases <- as.vector(unlist(nucleotides[valids]))
                most_common <- names(sort(table(valid_bases), decreasing = TRUE)[1])

                ambiguity_ratio <- (length(to_modify) / length(seqc))
                if (ambiguity_ratio >= max_ambiguity) {
                    if (log_operation) {
                        for (t in to_modify) {
                            cat(p, "\t", names(seqc)[t], "\t", seqc[[t]][p],
                                "\t", "Skipped - >= max_ambiguity", "\n",
                                file = log_file, append = TRUE)
                        }
                    }
                    if (output_progress) {
                        setTxtProgressBar(pb, p)
                    }
                    next
                }
                if (log_operation) {
                    for (t in to_modify) {
                        cat(p, "\t", names(seqc)[t], "\t", seqc[[t]][p],
                            "\t", most_common, "\n",
                            file = log_file, append = TRUE)
                    }
                }
                for (t in to_modify) {
                    seqc[[t]][p] <- most_common
                }
            }
            if (output_progress) {
                close(pb)
            }
        } else if (replace_method == "most_common_parallel") {
            bp$progressbar <- output_progress

            to_modify <- bplapply(seq_len(length(seqc[[1]])), function(pos, seqc, accepted_char){
                all(unlist(lapply(seqc, function(x) x[pos])) %in% accepted_char)
                
                
                
            }, seqc = seqc, accepted_char = accepted_char, BPPARAM = bp)

            ind_modify <- which(to_modify == FALSE)

            replacements <- bplapply(ind_modify, function(ind, seqc, accepted_char, max_ambiguity){
                bases <- lapply(seqc, `[`, ind)
                valids <- which(bases %in% accepted_char)
                to_modify <- which(!bases %in% accepted_char)

                ambiguity_ratio <- (length(to_modify) / length(seqc))
                if (ambiguity_ratio >= max_ambiguity) {
                    return("-")
                }
                b <- data.table(table(as.vector(unlist(bases[valids]))))
                colnames(b) <- c("base", "count")
                most_common <- b[order(-"count", "base")][1,]$base
                return(most_common)
            }, seqc = seqc, accepted_char = accepted_char, max_ambiguity = max_ambiguity, BPPARAM = bp)

            modified_seqc <- bplapply(names(seqc), function(seq_name, seqc, ind_modify, replacements, accepted_char, log_operation){
                modified_seq <- seqc[[seq_name]]
                for (ind in seq_len(length(ind_modify))) {
                    if (!modified_seq[ind_modify[[ind]]] %in% accepted_char){
                        if (replacements[[ind]] != "-"){
                            if (log_operation) {
                                cat(ind_modify[[ind]], "\t", seq_name, "\t", seqc[[seq_name]][ind_modify[[ind]]],
                                    "\t", replacements[[ind]], "\n",
                                    file = paste0("replace_", seq_name, ".tmp"), append = TRUE)
                            }
                            modified_seq[ind_modify[[ind]]] <- replacements[[ind]]
                        } else {
                            cat(ind_modify[[ind]], "\t", seq_name, "\t", seqc[[seq_name]][ind_modify[[ind]]],
                                "\t", "Skipped - >= max_ambiguity", "\n",
                                file = paste0("replace_", seq_name, ".tmp"), append = TRUE)
                        }
                    }
                }
                return(modified_seq)
            }, seqc = seqc, ind_modify = ind_modify, replacements = replacements,
            accepted_char = accepted_char, log_operation = log_operation, BPPARAM = bp)

            names(modified_seqc) <- names(seqc)
            if (log_operation) {
                tmp_logs <- list.files(pattern = "replace.*.tmp")
                lapply(tmp_logs, function(x) {
                    cat(paste(paste(readLines(x), collapse="\n"), "\n"), file = log_file, append = TRUE)
                    file.remove(x)
                })
            }

            seqc <- modified_seqc
        } else if (replace_method == "random") {
            if (output_progress) {
                pb <- txtProgressBar(min = 0, max = length(seqc[[1]]),
                    initial = 0, style = 3)
            }
            for (p in seq_len(length(seqc[[1]]))) {
                # All the nucleotides at this position
                nucleotides <- lapply(seqc, `[[`, p)
                # These need to be modified
                to_modify <- which(!nucleotides %in% accepted_char)

                if (max_ambiguity > 0 && max_ambiguity < 1) {
                    ambiguity_ratio <- (length(to_modify) / length(seqc))
                    if (ambiguity_ratio >= max_ambiguity) {
                        if (log_operation) {
                            for (t in to_modify) {
                                cat(p, "\t", names(seqc)[t], "\t", seqc[[t]][p],
                                    "\t", "Skipped - >= max_ambiguity", "\n",
                                    file = log_file, append = TRUE)
                            }
                        }
                        if (output_progress) {
                            setTxtProgressBar(pb, p)
                        }
                        next
                    }
                }

                # All other valid uncleotides
                valid_replacements <- unlist(
                    unique(
                        nucleotides[which(nucleotides %in% accepted_char)]
                    )
                )

                # Iterate through isolate that needs to be modified at position p
                for (t in to_modify) {

                    if (seqc[[t]][p] == "N") {
                        # If N is any bases, it can be any of the accepted characters #nolint
                        if (N_is_any_base) {
                            candidates <- unique(
                                c(translation[[seqc[[t]][p]]],
                                valid_replacements)
                            )
                        } else {
                            # Otherwise, it is one of the valid bases
                            # shown in at least one of the isolates
                            candidates <- valid_replacements
                        }
                    } else {
                        # If it's not N, just refer to the translation list
                        candidates <- unique(translation[[seqc[[t]][p]]])
                    }
                    replacement <-
                        sample(candidates, 1, replace = TRUE)
                    if (log_operation) {
                        cat(p, "\t", names(seqc)[t], "\t", seqc[[t]][p],
                            "\t", replacement, "\n",
                            file = log_file, append = TRUE)
                    }
                    seqc[[t]][p] <- replacement
                }
                if (output_progress) {
                    setTxtProgressBar(pb, p)
                }
            }
            if (output_progress) {
                close(pb)
            }
        } else {
            stop("Unknown replace method")
        }
        return(seqc)
    }