# seed <- as.numeric(Sys.time())

devtools::load_all("./MinSNPs")
vivax <- read_fasta("Pv_Nature_communications_2018_VQSLOD3.min0.fasta")
seed <- 2021
set.seed(seed)

# //TODO add on other IUPAC.
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


calculate_percent_w_amb_fail <- function(pattern, goi) {
    target_seqs <- character()
    translation[["N"]] <- c("A", "C", "T", "G")
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
    }
    all_unique_pattern <- unlist(unique(pattern))
    goi_patterns <- unlist(unique(pattern[goi]))
    non_goi_patterns <- all_unique_pattern[
        !all_unique_pattern %in% goi_patterns]

    for (ngp in non_goi_patterns) {
        nucleotides <- unlist(strsplit(ngp, split = ""))
        if (length(nucleotides[!nucleotides %in% c("A", "C", "T", "G")])
        == length(nucleotides)) {
            cat("FAILED: ", ngp, "\n")
        }
    }
    if (! pattern[[isolate]] %in% target_seqs) {
        target_seqs <- c(target_seqs, pattern[[isolate]])
    }
    isolate_w_goi_pattern <- length(which(pattern %in% target_seqs))
    failed_to_discriminate <- isolate_w_goi_pattern - length(goi)
    result <- 1 - (failed_to_discriminate / (length(pattern) - length(goi)))
    return(result)
}

resolve_IUPAC_missing <- function(seqc, N_is_any_base = FALSE, # nolint
    log_operation = TRUE, log_file = "replace.log", output_progress = TRUE) {
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
    if (output_progress) {
        pb = txtProgressBar(min = 0, max = length(seqc[[1]]),
            initial = 0, style = 3)
    }

    for (p in seq_len(length(seqc[[1]]))) {
        # All the nucleotides at this position
        nucleotides <- lapply(seqc, `[[`, p)
        # These need to be modified
        to_modify <- which(!nucleotides %in% accepted_char)

        # All other valid uncleotides
        valid_replacements <- unlist(
            unique(
                nucleotides[which(nucleotides %in% accepted_char)]
            )
        )

        # Iterate through isolate that needs to be modified at position p
        for (t in to_modify) {

            if (seqc[[t]][p] == "N") {
                # If N is any bases, it can be any of the accepted characters
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
    return(seqc)
}


### REMOVE BELOW
resolve_IUPAC_missing <- function(seqc, N_is_any_base = FALSE) { # nolint
    accepted_char <- c("A", "C", "T", "G")

    if (N_is_any_base) {
        translation[["N"]] <- accepted_char
    }

    for (p in seq_len(length(seqc[[1]]))) {

        nucleotides <- lapply(seqc, `[[`, p)
        to_modify <- which(!nucleotides %in% accepted_char)
        valid_replacements <- unlist(
            unique(nucleotides[which(nucleotides %in% accepted_char)]))
        cat("Position:", p, "\n", file = "replace.log", append = TRUE)
        for (t in to_modify) {
            print(seqc[[t]][p])
            print(translation[[seqc[[t]][p]]])
            if (seqc[[t]][p] == "N") {
                if (N_is_any_base) {
                    candidates <- valid_replacements
                } else {
                    candidates <- unique(c(translation[[seqc[[t]][p]]],
                    valid_replacements))
                }
            } else {
                candidates <- unique(translation[[seqc[[t]][p]]])
            }
            replacement <-
                sample(candidates, 1, replace = TRUE)
            cat("Replaced for,", names(seqc)[t], ":",
                seqc[[t]][p], "with", replacement, "\n",
                    file = "replace.log", append = TRUE)
            seqc[[t]][p] <- replacement
        }
    }
    return(seqc)
}

r_vivax <- resolve_IUPAC_missing(vivax)
write_fasta(r_vivax, "n_vivax.fasta")


resolve_missing_IUPAC <- function(seqc, bp=BiocParallel::SerialParam(),  #nolint
    dash_ignore=TRUE, accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE) {

    nuc_seq <- bplapply(seq_len(length(seqc[[1]])), function(pos) {
        nucleotides <- lapply(seqc, `[[`, pos)
        which(nucleotides )
        return(seqc, `[[`, pos)
    }, BPPARAM = bp)
}