#' \code{cal_fn}
#'
#' @description
#' \code{cal_fn} is used to check if the proportion of false negative
#' fastas and metas are compatible.
#' @param pattern the pattern from \code{generate_pattern}
#' @param goi the group of interest (names of isolates)
#' @param target the target sequence(s)
#' @return proportion: no. false negative/number of isolates
#' @export
cal_fn <- function(pattern, goi, target) {
    target_seqs <- c(target)

    if (length(target) < 1) {
        stop("Target(s) need to be defined")
    }
    snp_length <- length(strsplit(pattern[[1]], split = "")[[1]])
    for (t in target) {
        t_length <- length(strsplit(t, split = "")[[1]])
        if (t_length < snp_length) {
            stop("Target sequence needs to be as long as the number of SNPs")
        }
    }

    false_negatives <- c()
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
        if (! pattern[[isolate]] %in% target_seqs) {
            false_negatives <- c(false_negatives, isolate)
        }
    }
    result <- length(false_negatives) / length(pattern)

    return(list(result = result, additional_data = false_negatives))
}

#' \code{cal_fp}
#'
#' @description
#' \code{cal_fp} is used to check if the proportion of false positive
#' fastas and metas are compatible.
#' @param pattern the pattern from \code{generate_pattern}
#' @param goi the group of interest (names of isolates)
#' @param target the target sequence(s)
#' @return proportion: no. false positive/number of isolates
#' @export
cal_fp <- function(pattern, goi, target) {
    target_seqs <- c(target)

    if (length(target) < 1) {
        stop("Target(s) need to be defined")
    }
    snp_length <- length(strsplit(pattern[[1]], split = "")[[1]])
    for (t in target) {
        t_length <- length(strsplit(t, split = "")[[1]])
        if (t_length < snp_length) {
            stop("Target sequence needs to be as long as the number of SNPs")
        }
    }

    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
    }
    isolate_w_goi_pattern <- names(pattern[pattern %in% target_seqs])
    false_positives <- isolate_w_goi_pattern[! isolate_w_goi_pattern %in% goi]
    result <- length(false_positives) / length(pattern)

    return(list(result = result, additional_data = false_positives))
}