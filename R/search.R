#' \code{generate_pattern}
#'
#' @description
#' \code{generate_pattern} is used to generate pattern for calculation.
#' @param seqc list of sequences
#' @param ordered_index list of indexes for the pattern in the order
#' @param append_to existing patterns to append to
#' @return Will return concatenated list of string for searching.
generate_pattern <- function(seqc, ordered_index=c(), append_to=list()) {
    pattern <- append_to
    for (isolate in names(seqc)) {
        current <-
            paste(seqc[[isolate]][ordered_index], collapse = "")
        pattern[[isolate]] <- paste(append_to[[isolate]],
            current, sep = "")
    }
    return(pattern)
}

#' \code{calculate_percent}
#'
#' @description
#' \code{calculate_percent} is used to calculate dissimilarity index,
#' proportion of isolates not in goi that have been discriminated against.
#' 1 being all and 0 being none.
#' @param pattern list of sequences
#' @param goi group of interest
#' @return Will return the dissimilarity index of the list of patterns.
calculate_percent <- function(pattern, goi) {
    target_seqs <- character()
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
        if (pattern[[isolate]] %in% target_seqs) {
            target_seqs <- c(target_seqs, pattern[[isolate]])
        }
    }
    isolate_w_goi_pattern <- length(which(pattern %in% target_seqs))
    failed_to_discriminate <- isolate_w_goi_pattern - length(goi)
    result <- 1 - (failed_to_discriminate / (length(pattern) - length(goi)))
    return(result)
}

#' \code{calculate_simpson}
#'
#' @description
#' \code{calculate_simpson} is used to calculate Simpson's index.
#' Which is in the range of 0-1, where the greater the value,
#' the more diverse the population.
#' @inheritParams calculate.percent
#' @return Will return the Simpson's index of the list of patterns.
calculate_simpson <- function(pattern) {
    unique_sequences <- unlist(unique(pattern))
    no_isolate_in_group <- list()

    for (seqs in unique_sequences) {
        no_isolate_in_group[[seqs]] <-
            length(which(pattern == seqs))
    }

    d_numerator <- Reduce(f = function(prevs, current) {
        return(prevs + (current * (current - 1)))
    }, x = no_isolate_in_group, init = 0)

    d_denominator <- length(pattern) *
        (length(pattern) - 1)

    if (d_numerator == 0 || d_denominator == 0) {
        return(0)
    }

    simpson_index <- 1 - (d_numerator / d_denominator)
    return(simpson_index)
}

#' \code{find_optimised_snps}
#'
#' @description
#' \code{find_optimised_snps} is used to find optimised SNPs set.
#' @param seqc list of sequences
#' @param bp BiocParallel backend
#' @param criteria either `simpson` or `percent`
#' @param goi group of interest, if creteria is percent,
#' must be specified, ignored otherwise
#' @param accept_multiallelic whether include positions with > 2 states
#' @param number_of_result number of results to return
#' @param max_depth maximum depth to go before terminating
#' @param included_positions included positions
#' @param excluded_positions excluded positions
#' @param iterate_included whether to calculate index
#' at each level of the included SNPs
#' @importFrom BiocParallel SerialParam
#' @return Will return the resolution-optimised SNPs set, based on the criteria.
find_optimised_snps <- function(seqc, criteria = "simpson", goi = c(),
    accept_multiallelic = FALSE, number_of_result = 1, max_depth = 1,
    included_positions = c(), excluded_positions = c(),
    iterate_included = FALSE, bp = SerialParam()) {

    print("...")
}