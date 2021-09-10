#' \code{generate_pattern}
#'
#' @description
#' \code{generate_pattern} is used to generate pattern for calculation.
#' @param seqc list of sequences
#' @param ordered_index list of indexes for the pattern in the order
#' @param append_to existing patterns to append to
#' @keywords internal
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
        if (! pattern[[isolate]] %in% target_seqs) {
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
#' @inheritParams calculate_percent
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

#' \code{check_percent}
#'
#' @description
#' \code{check_percent} is used to check if parameters needed by
#' \code{calculate_percent} are all present.
#' @param list_of_parameters is a list of parameter passed
#' to functions that will perform the calculation
check_percent <- function(list_of_parameters) {
    return(length(list_of_parameters[["goi"]]) > 0)
}


view_simpson <- function(result, ...) {
    additional_args <- list(...)[[1]]
    if (exists(result$seqc_name)) {
        seqc <- get(result$seqc_name)
    } else if (!is.null(additional_args[["seqc"]])) {
        seqc <- additional_args[["seqc"]]
    } else {
        stop("Unable to find the sequences used to generate result")
    }
    if (class(seqc) == "processed_seqs") {
        seqc <- seqc$seqc
    }

    results <- result$results
    number_of_result <- length(results)
    isolates_in_groups <- list()
    for (n in seq_len(number_of_result)) {
        isolates_in_groups[[n]] <- list()
        result_levels <- names(results[[n]])
        result_max_depth <- result_levels[length(result_levels)]
        ordered_index <- as.numeric(
            strsplit(result_max_depth, split = ",")[[1]])
        patterns <- generate_pattern(seqc, ordered_index = c(ordered_index))
        unique_sequences <- unlist(unique(patterns))

        isolates_in_groups[[n]][["groups"]] <- list()
        for (seqs in unique_sequences) {
            isolates_in_groups[[n]][["groups"]][[seqs]] <-
                paste(names(which(patterns == seqs)), collapse = ", ")
        }
        isolates_in_groups[[n]][["result"]] <- results[[n]]
    }
    return(isolates_in_groups)
}

view_percent <- function(result, ...) { # nolint
    additional_args <- list(...)[[1]]
    if (exists(result$seqc_name) && !is.null(result$goi)) {
        seqc <- get(result$seqc_name)
    } else if (!is.null(additional_args[["seqc"]]) && !is.null(result$goi)) {
        seqc <- additional_args[["seqc"]]
    } else {
        stop("Unable to find the sequences used to generate result OR GOI")
    }

    if (class(seqc) == "processed_seqs") {
        seqc <- seqc$seqc
    }
    results <- result$results
    number_of_result <- length(results)
    goi <- result$goi
    isolates_in_groups <- list()

    for (n in seq_len(number_of_result)) {
        target_seqs <- character()
        isolates_in_groups[[n]] <- list()
        result_levels <- names(results[[n]])
        result_max_depth <- result_levels[length(result_levels)]
        ordered_index <- as.numeric(
            strsplit(result_max_depth, split = ",")[[1]])
        patterns <- generate_pattern(seqc, ordered_index = c(ordered_index))
        unique_sequences <- unlist(unique(patterns))

        isolates_in_groups[[n]][["groups"]] <- list()
        for (isolate in goi) {
            if (! patterns[[isolate]] %in% target_seqs) {
                target_seqs <- c(target_seqs, patterns[[isolate]])
            }
        }
        isolates_w_goi_pattern <- names(which(patterns %in% target_seqs))
        failed_to_discriminate <- which(!isolates_w_goi_pattern %in% goi)

        for (seqs in unique_sequences) {
            if (seqs %in% target_seqs) {
                seqs_name <- paste("*target* -", seqs)
            } else {
                seqs_name <- seqs
            }
            isolates_in_groups[[n]][["groups"]][[seqs_name]] <-
                paste(names(which(patterns == seqs)), collapse = ", ")
        }
        isolates_in_groups[[n]][["residual"]] <- paste(failed_to_discriminate,
            collapse = ", ")
        isolates_in_groups[[n]][["result"]] <- results[[n]]
    }
    return(isolates_in_groups)
}

#' \code{get_metric_fun}
#'
#' @description
#' \code{get_metric_fun} is used to get the metrics function
#' and required parameters. Additional metric may set by
#' assigning to `MinSNPs_metrics` variable.
#' @param metric_name name of the metric, by default percent/simpson
#' @return a list, including the function to calculate the
#' metric based on a position (`calc`), and function to check for
#' additional parameters the function need (`args`)
get_metric_fun <- function(metric_name) {
    if (! exists("MinSNPs_metrics")) {
        MinSNPs_metrics <- list( #nolint
            "percent" = list("calc" = calculate_percent,
                "args" = check_percent, "view" = view_percent),
            "simpson" = list("calc" = calculate_simpson,
                "view" = view_simpson)
        )
    }
    return(MinSNPs_metrics[[metric_name]])
}

#' \code{find_optimised_snps}
#'
#' @description
#' \code{find_optimised_snps} is used to find optimised SNPs set.
#' @param seqc list of sequences, either passed directly from
#' \code{process_allele} or \code{read_fasta} or equivalence
#' @param bp BiocParallel backend. Rule of thumbs: use MulticoreParam()
#' @param metric either `simpson` or `percent`
#' @param goi group of interest, if creteria is percent,
#' must be specified, ignored otherwise
#' @param accept_multiallelic whether include positions with > 1 state in goi
#' @param number_of_result number of results to return, 0 will be coerced to 1
#' @param max_depth maximum depth to go before terminating,
#' 0 means it will only calculate the metric for included position
#' @param included_positions included positions
#' @param excluded_positions excluded positions
#' @param iterate_included whether to calculate index
#' at each level of the included SNPs
#' @param ... other parameters as needed
#' @importFrom BiocParallel SerialParam bplapply
#' @return Will return the resolution-optimised SNPs set, based on the metric.
#' @export
find_optimised_snps <- function(seqc, metric = "simpson", goi = c(),
    accept_multiallelic = TRUE, number_of_result = 1, max_depth = 1,
    included_positions = c(), excluded_positions = c(),
    iterate_included = FALSE, bp = SerialParam(), ...) {

    # Define parameters
    positions <- seq_len(length(seqc[[1]]))
    result <- list()
    original_excluded <- excluded_positions
    included <- ifelse(length(included_positions) > 0, TRUE, FALSE)
    included_reached_1 <- FALSE
    additional_exclude <- numeric()
    # Automatically exclude ignored positions from processed sequences
    if (class(seqc) == "processed_seqs") {
        excluded_positions <- c(excluded_positions, seqc[["ignored_position"]])
        sequences <- seqc[["seqc"]]
    } else {
        sequences <- seqc
    }
    positions <- positions[! positions %in% excluded_positions]

    # Add multiallelic position to exclusion if percent mode
    if (!accept_multiallelic) {
        additional_exclude <- bplapply(positions, check_multistate,
            seqc[goi], BPPARAM = bp)
        names(additional_exclude) <- positions
        additional_exclude <- additional_exclude[additional_exclude == TRUE]
        additional_exclude <- as.numeric(names(additional_exclude))
        positions <- positions[! positions %in%
            additional_exclude]
    }

    # Check if any of the included positions are in excluded positions
    if (any(included_positions %in% excluded_positions)) {
        errors <- included_positions[included_positions %in% excluded_positions]
        stop(paste(errors, collapse = ","),
            " found in both included & excluded positions")
    }

    # Check if all required parameters are provided
    all_parameters <- list(...)
    all_parameters[["goi"]] <- goi

    metric_fun <- get_metric_fun(metric)[["calc"]]
    check_args <- get_metric_fun(metric)[["args"]]
    if (! is.null(check_args)) {
        if (! check_args(all_parameters)) {
            stop("GOI must be defined in percent mode")
        }
    }

    # Perform analysis on the included positions first
    existing_pattern <- list()
    if (included) {
        if (iterate_included) {
            iterations <- list()
            for (n in seq_along(included_positions)) {
                iterations[[n]] <- included_positions[1:n]
            }
            scores <- bplapply(iterations, cal_met_snp, metric = metric_fun,
                seqc = sequences, all_parameters, BPPARAM = bp)

            existing_pattern <- scores[[
                length(included_positions)]][["patterns"]]
            scores <- bplapply(scores, function(score) {
                return(score[["result"]])
            }, BPPARAM = bp)
            names(scores) <- bplapply(iterations, paste, collapse = ",",
                BPPARAM = bp)
            result <- scores
        } else {
            score <- cal_met_snp(included_positions, metric_fun,
                    sequences, all_parameters)
            result[[paste(included_positions, collapse = ",")]] <-
                score[["result"]]
            existing_pattern <- score[["patterns"]]
        }
    excluded_positions <- c(excluded_positions, included_positions)
    if (result[[paste(included_positions, collapse = ",")]] >= 1) {
            included_reached_1 <- TRUE
        }
    }

    if (max_depth == 0 ||
        included_reached_1
        ) {
        all_result <- list(result = result)
    } else {
        all_parameters[["existing_pattern"]] <- existing_pattern
        branch_result <- branch_and_search(included_positions,
            c(additional_exclude, excluded_positions),  sequences, metric_fun,
            number_of_result, max_depth,
            bp, all_parameters)
        all_result <- bplapply(branch_result, function(branch_r, inc_r) {
            return(c(inc_r, branch_r))
        }, inc_r = result)
    }

    return(list(results = all_result,
        excluded_positions = sort(unique(
            c(additional_exclude, original_excluded))),
        included_positions = included_positions,
        seqc_name = deparse(substitute(seqc)), goi = goi,
        all_sequences = names(sequences),
        max_depth = max_depth, metric = metric))
}

#' \code{branch_and_search}
#'
#' @description
#' \code{branch_and_search} is the actual function used to
#' find optimised SNPs set. This function is called by
#' \code{find_optimised_snps}, after preprocessing.
#' @inheritParams find_optimised_snps
#' @param starting_positions the starting positions that is already in
#' the SNP set.
#' @keywords internal
#' @importFrom BiocParallel SerialParam bplapply
#' @return Will return the resolution-optimised SNPs set, based on the metric.
branch_and_search <- function(starting_positions = c(),
    excluded_positions = c(), seqc, metric,
    number_of_result = 1, max_depth = 1,
    bp = SerialParam(), ...) {

    # Derived parameters
    additional_args <- list(...)[[1]]
    current_level <- length(starting_positions)
    traversed <- list()
    existing_pattern <- list()

    # Function to generate list of position to iterate through
    get_positions_to_search <- function(seqc_len, excluded_pos, traversed) {
        positions <- seq_len(seqc_len)
        positions <- positions[! positions %in% c(excluded_pos, traversed)]
        return(positions)
    }

    # Calculate metric
    positions <- get_positions_to_search(length(seqc[[1]]),
        excluded_positions, starting_positions)
    scores <- bplapply(positions,
        cal_met_snp, metric = metric, seqc = seqc,
        additional_args,
            BPPARAM = bp)
    depth_1 <- bplapply(scores, function(score) {
            return(score[["result"]])
    }, BPPARAM = bp)

    # Sorting for selection at this depth
    names(depth_1) <- positions
    position_order <- order(unlist(depth_1), decreasing = TRUE)

    result_d1 <- list()
    selected_positions <- numeric()

    # For each result, create the 1st selected SNP
    for (n in seq_len(number_of_result)) {
        traversed[[n]] <- c(starting_positions,
            positions[position_order[n]][[1]])
        result_d1[[
                paste(traversed[[n]], collapse = ",")]] <-
                    depth_1[position_order][[n]]
        existing_pattern[[n]] <- scores[[position_order[n]
            ]][["patterns"]]
        selected_positions <- c(selected_positions,
            positions[position_order[n]][[1]])
    }
    current_level <- current_level + 1

    if (number_of_result > 1) {
        multi_result <- list()
        # Loop to fulfil request result
        for (n in seq_len(number_of_result)) {
            # Only continue selecting another SNP if
            # max_depth has not been reached AND
            # if the index is not already 1
            if (max_depth - 1 > 0 &&
                result_d1[[paste(traversed[[n]], collapse = ",")]] < 1) {
                additional_args[["existing_pattern"]] <- existing_pattern[[n]]
                multi_result[[paste("result", n)]] <- branch_and_search(
                    starting_positions = c(starting_positions,
                        selected_positions[n]),
                    excluded_positions = c(excluded_positions,
                        selected_positions[1:n]),
                    seqc, metric, number_of_result = 1,
                    max_depth = (max_depth - 1), bp = bp, additional_args
                )
            }
            # Combine the multiple result with the result at 1st level
            multi_result[[paste("result", n)]] <- c(
                result_d1[n],
                multi_result[[paste("result", n)]][["result"]])
        }
        return(multi_result)
    } else {
        # Loop through multiple depth if max_depth has not been reached AND
        # if index is not already 1
        while (((current_level - length(starting_positions)) < max_depth) &&
               (depth_1[position_order][[1]] < 1))  {

            additional_args[["existing_pattern"]] <- existing_pattern[[1]]
            positions <- get_positions_to_search(length(seqc[[1]]),
                c(excluded_positions, selected_positions[1]),
                traversed[[1]])

            scores <- bplapply(positions,
                cal_met_snp, metric = metric, seqc = seqc,
                additional_args,
                    BPPARAM = bp)
            depth_1 <- bplapply(scores, function(score) {
                    return(score[["result"]])
            }, BPPARAM = bp)

            names(depth_1) <- positions
            position_order <- order(unlist(depth_1), decreasing = TRUE)


            traversed[[1]] <- c(traversed[[1]],
                positions[position_order[1]][[1]])
            result_d1[[
                    paste(traversed[[1]], collapse = ",")]] <-
                        depth_1[position_order][[1]]
            existing_pattern[[1]] <- scores[[position_order[1]
                ]][["patterns"]]

            current_level <- current_level + 1
        }
        return(list(result = result_d1))
    }
}

#' \code{check_multistate}
#'
#' @description
#' \code{check_multistate} is used to remove positions where there are
#' more than 1 state within the group of interest.
#' @param position position to check
#' @param sequences sequences from group of interest
#' @return return `TRUE` if the position contains multistate otherwise `FALSE`
#' @keywords internal
check_multistate <- function(position, sequences) {
    pattern <- generate_pattern(sequences, c(position))
    if (length(unique(pattern)) > 1) {
        return(TRUE)
    }
    return(FALSE)
}

#' \code{cal_met_snp}
#'
#' @description
#' \code{cal_met_snp} is used to calculate the metric at each position
#' @inheritParams find_optimised_snps
#' @inheritParams check_multistate
#' @return return the value at that position,
#' as well as base pattern for next iteration.
#' @keywords internal
cal_met_snp <- function(position, metric, seqc, ...) {
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
    sub_args <- additional_args[all_expected_args]
    sub_args[["pattern"]] <- pattern
    res <- do.call(metric_function, args = sub_args)
    return(list(result = res, patterns = pattern))
}