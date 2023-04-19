#' \code{output_result}
#'
#' @description
#' \code{output_result} is used to present the result
#' and save the result as tsv.
#' @param result is the result from \code{find_optimised_snps}
#' @param view how to present the output,
#' "csv" or "tsv" will be saved as a file.
#' Otherwise, printed to console.
#' @param ... if view is "tsv" or "csv", file name can be passed, e.g.,
#' file_name = "result.tsv", otherwise, file is saved as <timestamp>.tsv.
#' @return NULL, result either printed or saved as tsv.
#' @export
output_result <- function(result, view = "", ...) {
    additional_args <- list(...)

    get_seq_obj <- function(seqc_name, additional_args) {
        if (exists(seqc_name)) {
            seqc_obj <- get(seqc_name)
        } else if (is.null(additional_args[["seqc"]])) {
            seqc_obj <- additional_args[["seqc"]]
        } else {
            seqc_obj <- list()
            seqc_obj$ignored_position <-
                "Some error occured, unable to retrieve"
        }
        return(seqc_obj)
    }

    if (is.null(get_metric_fun(result$metric)[["view"]])) {
        stop("The view function is not defined")
    }

    if ("file_name" %in% names(additional_args)) {
        file_name <- additional_args[["file_name"]]
    } else {
        if (view == "csv") {
            o_view <- "tsv"
        } else{
            o_view <- view
        }
        file_name <- paste(gsub(" ", "_", as.character(Sys.time())),
            ".", o_view, sep = "")
        if (view == "rds" || view == "csv" || view == "tsv")
        cat("No file name supplied, saved as: ", file_name, "\n")
    }

    if (view == "rds") {
        seqc_obj <- get_seq_obj(result$seqc_name, additional_args)
        anlysis_result <- list(seqc = seqc_obj, result = result)
        saveRDS(anlysis_result, file = file_name)
        return(paste("Wrote to:", file_name))
    }
    if (view == "" || view == "cmd") {
        file_name <- ""
    }

    parsed_view <- get_metric_fun(result$metric)[["view"]](result,
        additional_args)
    seqc_obj <- get_seq_obj(result$seqc_name, additional_args)
    output_parsed(parsed_view, seqc_obj, result, file_name)
}

#' \code{output_parsed}
#'
#' @description
#' \code{output_parsed} is subfunction used by \code{output_result}
#' to present the result and save the result as tsv.
#' @param result is the result from \code{find_optimised_snps}
#' @param file_name file name, either tsv or "" for terminal output
#' @param parsed_view the result from \code{get_metric_fun}, passed from
#' \code{output_result}
#' @keywords internal
#' @param seqc_obj the sequences object.
#' @return NULL, result either printed or saved as tsv.
output_parsed <- function(parsed_view, seqc_obj, result, file_name) {
    write_output <- function(text, new_file = FALSE) {
        cat(paste(text, "\n", sep = ""),
            sep = "", file = file_name, append = !new_file)
    }
    for (n in seq_along(parsed_view)) {
        write_output(paste("Result - ", n, sep = ""),
            ifelse(n == 1, TRUE, FALSE))
        write_output("Position(s)\tScore")
        for (level in seq_along(parsed_view[[n]]$result)) {
            pos <- paste("\"", names(parsed_view[[n]]$result[level]),
                "\"", sep = "")
            val <- parsed_view[[n]]$result[[level]]
            write_output(paste(pos, val, sep = "\t"))
        }
        write_output("\nGroups")
        for (group in seq_along(parsed_view[[n]]$groups)) {
            group_seq <- names(parsed_view[[n]]$groups[group])
            isolates <- paste("\"", parsed_view[[n]]$groups[[group]],
                "\"", sep = "")
            write_output(paste(group_seq, isolates, sep = "\t"))
        }
        if ("residual" %in% names(parsed_view[[n]])) {
            write_output(paste("Residuals:",
                paste("\"", parsed_view[[n]][["residual"]], "\"", sep = ""),
                sep = "\t"))
        }
        other_details <- names(parsed_view[[n]])[
            startsWith(names(parsed_view[[n]]), "N.B.-")]
        for (other_val in other_details) {
            write_output(paste(paste(other_val, ":", sep = ""),
                paste(parsed_view[[n]][[other_val]], collapse = ", "),
                sep = "\t"))
        }
        write_output("\n")
    }

    write_output("Additional details")
    write_output(paste("Metric:\t", result$metric))
    write_output(paste("Excluded Positions:\t",
        paste("\"", paste(result$excluded, collapse = ", "), "\"", sep = "")
        ))
    write_output(paste("Excluded Positions From process_allele:\t",
        paste("\"", paste(seqc_obj$ignored_position, collapse = ", "),
            "\"", sep = "")
        ))
    write_output(paste("Included Positions:\t",
        paste("\"", paste(result$included, collapse = ", "), "\"", sep = "")))
    write_output(paste("Group of interest:\t",
        paste("\"", paste(result$goi, collapse = ", "), "\"", sep = "")))
    write_output(paste("All analysed sequences:\t",
        paste("\"", paste(result$all_sequences, collapse = ", "), "\"",
        sep = "")))
}