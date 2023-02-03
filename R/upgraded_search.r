#### UPGRADED SEARCH ####

#' \code{num_of_states}
#'
#' @description
#' \code{num_of_states} is used to identify the number of states for a position.
#' @param position position to check
#' @param sequences sequences from group of interest
#' @param bp BiocParallel backend, default to SerialParam()
#' @importFrom BiocParallel bplapply SerialParam
#' @return return a vector of states of length(sequences)
num_of_states <- function(position, sequences, bp = SerialParam()) {
    n_snps = seq_len(sequences[[1]])
    states <- bplapply(n_snps, function(x, sequences) {
        states <- lapply(sequences, `[`, x)
        return(length(unique(states)))
    }, sequences = sequences, BPPARAM = bp)
    return(unlist(states))
}

#' \code{cal_met_snp1}
#'
#' @description
#' \code{cal_met_snp1} is used to calculate the metric at each position
#' @inheritParams find_optimised_snps
#' @inheritParams num_of_states
#' @return return the value at that position,
#' as well as base pattern for next iteration.
#' @keywords internal
#' @export
cal_met_snp1 <- function(position, metric, seqc, ...) {
    additional_args <- list(...)[[1]]
    if ("existing_pattern" %in% names(additional_args)) {
        append_to <- additional_args[["existing_pattern"]]
    } else {
        append_to <- list()
    }
    pattern <- generate_pattern(seqc, c(position), append_to)
    metric_function <- match.fun(metric)
    metric_expected_args <- formals(metric_function)
    all_expected_args <- names(additional_args)[
        names(additional_args) %in% names(metric_expected_args)]
    if (!all(names(metric_expected_args)
        [!names(metric_expected_args) == "pattern"] %in%
            all_expected_args)) {

            stop("Metric function requires additional arguments: ",
              paste(names(metric_expected_args)[(!names(metric_expected_args)
              %in% all_expected_args) &
              (!names(metric_expected_args) == "pattern")], collapse = ", "))
    }
    sub_args <- additional_args[all_expected_args]
    sub_args[["pattern"]] <- pattern
    res <- do.call(metric_function, args = sub_args)
    if ("get_addn" %in% names(additional_args)) {
        if (additional_args[["get_addn"]]) {
            return(list(result = res$result,
                additional_data = res$additional_data))
        }
    }
    return(list(result = res$result))#, patterns = pattern)) ## Updated this
}