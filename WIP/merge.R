#' \code{merge_fasta}
#'
#' @description
#' \code{merge_fasta} is used to combine 2 fasta.
#' @param fasta_1 fasta read into memory to join
#' @param fasta_2 fasta read into memory to join
#' @param meta_1 meta file for `fasta_1`
#' @param meta_2 meta file for `fasta_2`
#' @param method how to join the 2 fasta, currently supported methods are:
#' inner, full
#' @return Will return a list containing a merged FASTA and a meta.
#' @export
merge_fasta <- function(fasta_1, fasta_2, meta_1, meta_2, method = "inner") {

    # Check fastas and metas are compatible
    if (!check_fasta_meta_mapping(fasta_1, meta_1)) {
        stop("There is problem with fasta_1 and meta_1")
    }
    if (!check_fasta_meta_mapping(fasta_2, meta_2)) {
        stop("There is problem with fasta_2 and meta_2")
    }

    merged_fasta <- NULL

    if (method == "inner") {
        merged_fasta <- "something"
    } else if (method == "full") {
        merged_fasta <- "Something else"
    } else {
        stop("Method is not supported")
    }

    return(merged_fasta)
}


### Same as the reference genome if not found.



#' \code{check_fasta_meta_mapping}
#'
#' @description
#' \code{check_fasta_meta_mapping} is used to check if
#' fastas and metas are compatible.
#' @param fasta the fasta read into memory to join
check_fasta_meta_mapping <- function(fasta, meta) {
    is_valid <- FALSE
    # 1. i.e. all positions in fasta_1 are mapped
    # 2. all positions are mapped to a unique genome position

    return(is_valid)
}