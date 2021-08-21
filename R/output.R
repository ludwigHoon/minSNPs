#' \code{output_result}
#'
#' @description
#' \code{output_result} is used to present the result
#' and save the result as csv.
#' @param result is the result from \code{find_optimised_snps}
#' @param view how to present the output, "csv" will be saved as a file.
#' Otherwise, printed to console.
#' @param ... if view is "csv", file name can be passed, e.g.,
#' file_name = "result.csv", otherwise, file is saved as <timestamp>.csv.
#' @return NULL, result either printed or saved as csv.
output_result <- function(result, view, ...) {
    additional_args <- list(...)
    if (is.null(minSNPs::get_metric_fun(result$metric)[["view"]])) {
        stop("The view function is not defined")
    }

    parsed_view <- minSNPs::get_metric_fun(result$metric)[["view"]](result,
        additional_args)
    if (view == "csv") {
        if ("file_name" %in% names(additional_args)) {
            file_name <- additional_args[["file_name"]]
        } else {
            file_name <- paste(gsub(" ", "_", as.character(Sys.time())),
                ".csv", sep = "")
            cat("No file name supplied, saved as: ", file_name, "\n")
        }
    } else {
        file_name <- ""
    }

    write_output <- function(text, new_file = FALSE) {
        cat(paste(text, "\n", sep = ""),
            sep = "", file = file_name, append = !new_file)
    }

    for (n in seq_along(parsed_view)) {
        write_output(paste("Result - ", n, sep = ""),
            ifelse(n == 1, TRUE, FALSE))
        write_output("Position(s)\tScore")
        for (level in seq_along(parsed_view[[n]]$result)) {
            pos <- names(parsed_view[[n]]$result[level])
            val <- parsed_view[[n]]$result[[level]]
            write_output(paste(pos, val, sep = "\t"))
        }
        write_output("\nGroups")
        for (group in seq_along(parsed_view[[n]]$groups)) {
            group_seq <- names(parsed_view[[n]]$groups[group])
            isolates <- parsed_view[[n]]$groups[[group]]
            write_output(paste(group_seq, isolates, sep = "\t"))
        }
        if ("residual" %in% parsed_view[[n]]) {
            write_output(paste(group_seq, isolates, sep = "\t"))
        }
        write_output("\n")
    }
    write_output("Additional details")
    if (exists(result$seqc_name)) {
        seqc_obj <- get(result$seqc_name)
    } else if (is.null(additional_args[["seqc"]])) {
        seqc_obj <- additional_args[["seqc"]]
    } else {
        seqc_obj <- list()
        seqc_obj$ignored_position <- "Some error occured, unable to retrieve"
    }
    write_output(paste("Metric:\t", result$metric))
    write_output(paste("Excluded Positions:\t",
        paste(result$excluded, collapse = ", ")))
    write_output(paste("Excluded Positions From process_allele:\t",
        paste(seqc_obj$ignored_position, collapse = ", ")))
    write_output(paste("Included Positions:\t",
        paste(result$included, collapse = ", ")))
    write_output(paste("Group of interest:\t",
        paste(result$goi, collapse = ", ")))
}