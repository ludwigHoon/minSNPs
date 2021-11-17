#' \code{merge_fasta}
#'
#' @description
#' \code{merge_fasta} is used to combine 2 fasta.
#' @param fasta_1 fasta read into memory to join
#' @param fasta_2 fasta read into memory to join
#' @param meta_1 meta file for `fasta_1` denoting all positions of SNPs
#' and position in reference genome
#' @param meta_2 meta file for `fasta_2` denoting all positions of SNPs
#' and position in reference genome
#' @param ref name of the reference genome (needs to be in both fasta files)
#' @param method how to join the 2 fasta, currently supported methods are:
#' inner, full
#' @return Will return a list containing a merged FASTA and a meta.
#' @export
merge_fasta <- function(fasta_1, fasta_2, meta_1, meta_2,
    ref, method = "full", bp = BiocParallel::SerialParam(), ...) {
    if (!ref %in% names(fasta_1)) {
        stop("fasta_1 must have the reference sequences")
    }
    if (!ref %in% names(fasta_2)) {
        stop("fasta_2 must have the reference sequences")
    }
    # Check fastas and metas are compatible
    if (!check_fasta_meta_mapping(fasta_1, meta_1)) {
        stop("There is problem with fasta_1 and meta_1")
    }
    if (!check_fasta_meta_mapping(fasta_2, meta_2)) {
        stop("There is problem with fasta_2 and meta_2")
    }

    repeated_isolate_name <- intersect(names(fasta_1), names(fasta_2))
    if (length(repeated_isolate_name[repeated_isolate_name != ref]) >= 1) {
        stop("Isolate name repeated: ",
            paste(repeated_isolate_name[repeated_isolate_name != ref],
                collapse = ", "))
    }
    result <- list()

    if (method == "inner") {
        result$merged_fasta <- list()
        result$merged_meta <- list()
        cat("Unimplemented at the moment\n")
    } else if (method == "full") {
       result <- full_merge_1(fasta_1, fasta_2, meta_1, meta_2, ref, bp,
        ... = ...)
    } else {
        stop("Method is not supported")
    }

    return(result)
}

#' \code{iterate_merge}
#'
#' @description
#' \code{iterate_merge} is used to combine > 2 fastas iteratively.
#' @param fastas list of fastas read into memory to join
#' @param metas list of metas read into memory to join
#' @param ref name of the reference genome (needs to be in both fasta files)
#' @param method how to join the 2 fasta, currently supported methods are:
#' inner, full
#' @return Will return a list containing a merged FASTA and a meta.
#' @export
iterate_merge <- function(fastas, metas, ref, method = "full",
                          bp = BiocParallel::SerialParam(), ...) {
    if (length(fastas) != length(metas)) {
        stop("fastas and metas must have the same length")
    }
    for (round in seq_len(length(fastas) - 1)) {
        cat("Running", round, "/", (length(fastas) - 1), "\n")
        if (round == 1) {
            result <- merge_fasta(fastas[[1]], fastas[[2]],
                metas[[1]], metas[[2]], ref, method, bp, ... = ...)
        } else {
            result <- merge_fasta(result$merged_fasta,
                fastas[[round + 1]], result$merged_meta,
                metas[[round + 1]], ref, method, bp, ... = ...)
        }
    }
    return(result)
}

#' \code{output_to_files}
#'
#' @description
#' \code{output_to_files} is write the result to files.
#' @param merged_result a list containing the merged fasta and meta.
#' @param filename filename to write to, will output
#' <filename>.fasta and <filename>.csv.
#' @return NULL, files written to filesystem
#' @export
output_to_files <- function(merged_result, filename = "merged") {
    write_fasta(merged_result$merged_fasta,
        paste(filename, ".fasta", sep = ""))
    write.csv(merged_result$merged_meta, paste(filename, ".csv", sep = ""),
        row.names = FALSE)
}

#' \code{full_merge_1}
#'
#' @description
#' \code{full_merge_1} is used to merge 2 fasta,
#' where a position exist only in 1 of the fasta, the fasta without allele
#' in that positions are given reference genome's allele at that position.
#' @inheritParams merge_fasta
#' @return merged fasta and meta
#' @export
full_merge_1 <- function(fasta_1, fasta_2, meta_1, meta_2, ref,
                       bp = BiocParallel::SerialParam(), ...) {
    args <- list(...)
    all_genome_positions <- sort(unique(c(meta_1$genome_position,
        meta_2$genome_position)))
    common_positions <- intersect(meta_1$genome_position,
        meta_2$genome_position)
    merged_meta <- data.frame(genome_position = all_genome_positions)
    merged_meta <- cbind(merged_meta,
        fasta_position = seq_len(nrow(merged_meta)))

    if (! is.null(args[["v"]])) {
        cat("Finished Merging Meta\n")
    }

    new_fasta_1 <- bplapply(merged_meta[, "genome_position"],
        function(g_pos) {
            final_seqs <- character()
            if (g_pos %in% meta_1$genome_position) {
                final_seqs <- lapply(fasta_1, `[`,
                    meta_1[meta_1$genome_position == g_pos,
                        "fasta_position"])
            } else {
                final_seqs <- as.list(rep(
                    fasta_2[[ref]][
                        meta_2[meta_2$genome_position == g_pos,
                        "fasta_position"]], length(fasta_1)))
                names(final_seqs) <- names(fasta_1)
            }
            return(final_seqs)
        }, BPPARAM = bp)
    new_fasta_1 <- bplapply(seq_len(length(names(fasta_1))),
        function(i) {
            return(
                unname(unlist(lapply(new_fasta_1, `[`, i)))
                )
            }, BPPARAM = bp
        )
    if (! is.null(args[["v"]])) {
        cat("Finished Generating New Fasta_1\n")
    }
    new_fasta_2 <- bplapply(merged_meta[, "genome_position"],
        function(g_pos) {
            final_seqs <- character()
            if (g_pos %in% meta_2$genome_position) {
                final_seqs <- lapply(fasta_2, `[`,
                    meta_2[meta_2$genome_position == g_pos,
                        "fasta_position"])
            } else {
                final_seqs <- as.list(rep(
                    fasta_1[[ref]][
                        meta_1[meta_1$genome_position == g_pos,
                        "fasta_position"]], length(fasta_2)))
                names(final_seqs) <- names(fasta_2)
            }
            return(final_seqs)
        }, BPPARAM = bp)
    new_fasta_2 <- bplapply(seq_len(length(names(fasta_2))),
        function(i) {
            return(
                unname(unlist(lapply(new_fasta_2, `[`, i)))
                )
            }, BPPARAM = bp
        )
    if (! is.null(args[["v"]])) {
        cat("Finished Generating New Fasta_2\n")
    }
    merged_fasta <- c(new_fasta_1, new_fasta_2)
    names(merged_fasta) <- c(names(fasta_1), names(fasta_2))
    return(list(merged_fasta = merged_fasta, merged_meta = merged_meta))
}


#' \code{full_merge}
#'
#' @description
#' \code{full_merge} is used to merge 2 fasta,
#' where a position exist only in 1 of the fasta, the fasta without allele
#' in that positions are given reference genome's allele at that position.
#' **Doesn't work for large dataset, hence the need for \code{full_merge_1}**
#' @inheritParams merge_fasta
#' @return merged fasta and meta
full_merge <- function(fasta_1, fasta_2, meta_1, meta_2, ref,
                       bp = BiocParallel::MulticoreParam(), ...) {
    args <- list(...)
    all_genome_positions <- sort(unique(c(meta_1$genome_position,
        meta_2$genome_position)))
    common_positions <- intersect(meta_1$genome_position,
        meta_2$genome_position)
    merged_meta <- data.frame(genome_position = all_genome_positions)
    merged_meta <- cbind(merged_meta,
        fasta_position = seq_len(nrow(merged_meta)))

    if (! is.null(args[["v"]])) {
        cat("Finished Merging Meta\n")
    }

    new_fasta_1 <- bplapply(names(fasta_1), function(isolate) {
        final_seqs <- character()
        for (genome_pos in merged_meta[, "genome_position"]) {
            if (genome_pos %in% meta_1$genome_position) {
                final_seqs <- c(final_seqs,
                    fasta_1[[isolate]][
                    meta_1[meta_1$genome_position == genome_pos,
                    "fasta_position"]])
            } else {
                final_seqs <- c(final_seqs,
                    fasta_2[[ref]][
                        meta_2[meta_2$genome_position == genome_pos,
                        "fasta_position"]])
            }
        }
        return(final_seqs)
    }, BPPARAM = bp)
    if (! is.null(args[["v"]])) {
        cat("Finished Generating New Fasta_1\n")
    }
    new_fasta_2 <- bplapply(names(fasta_2), function(isolate) {
        final_seqs <- character()
        for (genome_pos in merged_meta[, "genome_position"]) {
            if (genome_pos %in% meta_2$genome_position) {
                final_seqs <- c(final_seqs,
                    fasta_2[[isolate]][
                    meta_2[meta_2$genome_position == genome_pos,
                    "fasta_position"]])
            } else {
                final_seqs <- c(final_seqs,
                    fasta_1[[ref]][
                        meta_1[meta_1$genome_position == genome_pos,
                        "fasta_position"]])
            }
        }
        return(final_seqs)
    }, BPPARAM = bp)
    if (! is.null(args[["v"]])) {
        cat("Finished Generating New Fasta_2\n")
    }
    merged_fasta <- c(new_fasta_1, new_fasta_2)
    names(merged_fasta) <- c(names(fasta_1), names(fasta_2))
    return(list(merged_fasta = merged_fasta, merged_meta = merged_meta))
}



#' \code{check_fasta_meta_mapping}
#'
#' @description
#' \code{check_fasta_meta_mapping} is used to check if
#' fastas and metas are compatible.
#' @param fasta the fasta read into memory to join
#' @param meta the meta read into memory to join
#' @return TRUE/FALSE if the fasta and meta are compatible
check_fasta_meta_mapping <- function(fasta, meta) {
    is_valid <- FALSE
    if (length(unique(sapply(fasta, length))) > 1) {
        stop("fasta is not valid, found sequences of unequal length")
    }
    # 1. i.e. all positions in fasta are mapped
    # 2. all positions are mapped to a unique genome position
    if (all(c("genome_position", "fasta_position") %in% colnames(meta))) {
        is_valid <- TRUE && (length(fasta[[1]]) <= max(meta$fasta_position))
        is_valid <- is_valid && all(diff(meta$fasta_position) == 1)
        is_valid <- is_valid &&
            nrow(unique(meta[, c("genome_position", "fasta_position")])) ==
            nrow(meta)
    } else {
        stop("Meta must consist of `genome_position` and `fasta_position`")
    }
    return(is_valid)
}