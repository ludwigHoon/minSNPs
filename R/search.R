#' \code{generate_pattern} is used to generate pattern for calculation.
#' @param seqc list of sequences
#' @param ordered_index list of indexes for the pattern in the order
#' @param append_to existing patterns to append to
#' @return Will return concatenated list of string for searching.
generate_pattern <- function(seqc, ordered_index=c(), append_to=list()) {

}

#' \code{calculate.percent} is used to calculate dissimilarity index.
#' @param pattern list of sequences
#' @param goi group of interest
#' @return Will return the dissimilarity index of the list of patterns.
calculate.percent <- function(pattern, goi) { # nolint

}

#' \code{calculate.simpson} is used to calculate Simpson's index.
#' @inheritParams calculate.percent
#' @return Will return the Simpson's index of the list of patterns.
calculate.simpson <- function(pattern) { # nolint

}

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
#' @import BiocParallel
#' @return Will return the resolution-optimised SNPs set, based on the criteria.
find_optimised_snps <- function(seqc, criteria = "simpson", goi = c(),
    accept_multiallelic = FALSE, number_of_result = 1, max_depth = 1,
    included_positions = c(), excluded_positions = c(),
    iterate_included = FALSE, bp = BiocParallel::SerialParam()) {

    print("...")
}