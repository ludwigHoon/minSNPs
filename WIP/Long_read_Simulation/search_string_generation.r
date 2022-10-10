process_result_file <- function(result_filepath) {
    f <- file(result_filepath, "r")
    preparse <- readLines(f)
    result <- list()
    cur <- 1
    is_result <- FALSE
    for (i in seq_len(length(preparse))) {
        if (length(grep("^Result", preparse[[i]])) >= 1) {
            is_result <- TRUE
        }
        if (is_result) {
            if (preparse[[i + 1]] == "") {
                result[[cur]] <- as.numeric(
                    strsplit(
                        gsub(
                            "\"", "",
                            strsplit(preparse[[i]], split = "\t")[[1]][1]
                        ),
                        split = ", "
                    )[[1]]
                )
                cur <- cur + 1
                is_result <- FALSE
            }
        }
    }
    close(f)
    final_selected <- (result)
    return(final_selected)
}

process_variant_file <- function(variant_sites) {
    result <- c()
    read_cc <- c()
    con <- file(variant_sites, "r")
    while (TRUE) {
        data <- readLines(con, n = 3)
        if (length(data) == 0) {
            break
        }
        positions <- strsplit(strsplit(data, split = "\n")[[2]], split = ", ")
        read_cc <- c(read_cc, strsplit(data, split = "\n")[[1]])
        result <- c(result, unlist(positions))
    }
    close(con)
    print(read_cc)
    result <- sort(as.numeric(unique(result)))
    return(result)
}

#' \code{match_count}
#'
#' @description
#' \code{match_count} return the number of matches of the target string in the
#' given sequence
#' @param target the search target
#' @param search_from the sequence to search from
#' @return number of matches
#' @export
match_count <- function(target, search_from) {
    matches <- gregexpr(paste(target, collapse = ""), search_from)
    found <- 0
    if (length(matches[[1]]) == 1) {
        if (matches[[1]][1] != -1) {
            found <- 1
        }
    } else {
        found <- length(matches[[1]])
    }
    return(found)
}

#' \code{reverse_complement}
#'
#' @description
#' \code{reverse_complement} returns the reverse complement of the
#' given sequence
#' @param seq the sequence to reverse complement
#' @return reverse complemented sequence
#' @export
reverse_complement <- function(seq) {
    seq <- rev(strsplit(seq, split = "")[[1]])
    result <- lapply(seq, function(x) {
        if (toupper(x) == "A") {
            return("T")
        }
        if (toupper(x) == "C") {
            return("G")
        }
        if (toupper(x) == "T") {
            return("A")
        }
        if (toupper(x) == "G") {
            return("C")
        }
        return(x)
    })
    return(paste(result, collapse = ""))
}

#' \code{generate_kmers}
#'
#' @description
#' \code{generate_kmers} generate the kmer sequences of the given length
#' @param final_string the string to generate kmers
#' @param k the length of the kmer
#' @return a vector of kmers
#' @export
generate_kmers <- function(final_string, k) {
    if (is.character(final_string) && length(final_string) == 1) {
        final_string <- strsplit(final_string, split = "")[[1]]
    }
    kmer <- list()
    for (i in seq_len(length(final_string) - k + 1)) {
        kmer[[i]] <- paste(final_string[i:(i + k - 1)], collapse = "")
    }
    return(unlist(kmer))
}

#' \code{identify_overlaps}
#'
#' @description
#' \code{identify_overlaps} identify the SNPs that will overlap
#' the search strings generated from the targeted SNPs
#' @param selected_snps list of targeted SNPs
#' @param position_reference the mapping between
#' reference genome positions and orthologous SNP matrix positions
#' @param prev number of characters before the SNP
#' @param after number of characters after the SNP
#' @param position_type type of SNPs input, "fasta"
#' (orthologous SNP matrix based) or "genome"
#' (reference genome based); Default to "fasta"
#' @param extend_length whether to extend the search string
#' before and after the SNP and ignore overlapping SNPs
#' @return a list containing 2 dataframes,
#' (1) SNP tables - 
#' (2) overlap tables
#' @export
identify_overlaps <- function(selected_snps, position_reference, prev, after,
    position_type = "fasta", extend_length = TRUE) {
    snp_table <- data.frame(
        snp_id = c(), fasta_position = c(), genome_position = c(),
        string_start = c(), string_end = c(),
        stringsAsFactors = FALSE
    )
    overlap_table <- data.frame(
        snp_id = c(), overlaps_genome_pos = c(), overlaps_fasta_pos = c(),
        overlaps_string_pos = c(), stringsAsFactors = FALSE
    )

    #<TODO>
    selected_fasta <- c(selected_snps)
    #<TODO>

    for (snp in selected_fasta) {
        genome_pos <- position_reference[
            position_reference$fasta_position == snp, "genome_position"]
        string_start <- genome_pos - prev
        string_end <- genome_pos + after

        n_test <- 0
        snp_string_pos <- prev + 1
        # Check overlaps before
        overlaps_before <- position_reference[
            position_reference$genome_position
                >= (genome_pos - (n_test + prev)) &
            position_reference$genome_position
                < genome_pos,
        c("fasta_position", "genome_position")]

        if (extend_length) {
            while ((n_test + prev) - nrow(overlaps_before) < prev) {
                n_test <- n_test + (nrow(overlaps_before) - n_test)
                overlaps_before <- position_reference[
                    position_reference$genome_position
                        >= (genome_pos - (n_test + prev)) &
                    position_reference$genome_position
                        < genome_pos, c("fasta_position", "genome_position")
                ]
            }
            string_start <- string_start - n_test
            snp_string_pos <- snp_string_pos + n_test
            overlaps_before$s_pos <-
                (overlaps_before$genome_position - string_start)
        }

        n_test <- 0
        # Check overlaps after
        overlaps_after <- position_reference[
            position_reference$genome_position
                > genome_pos &
            position_reference$genome_position
                <= (genome_pos + (n_test + after)),
        c("fasta_position", "genome_position")]

        if (extend_length) {
            while ((n_test + after) - nrow(overlaps_after) < after) {
                n_test <- n_test + (nrow(overlaps_after) - n_test)
                overlaps_after <- position_reference[
                    position_reference$genome_position
                        > genome_pos &
                    position_reference$genome_position
                        <= (genome_pos + (n_test + after)),
                c("fasta_position", "genome_position")]
            }
            string_end <- string_end + n_test
            overlaps_after$s_pos <-
                (snp_string_pos + (overlaps_after$genome_position - genome_pos))
        }

        temp_overlap_table <- rbind(overlaps_before, overlaps_after)
        if (ncol(temp_overlap_table) >= 3) {
            colnames(temp_overlap_table) <-
                c("overlaps_fasta_pos", "overlaps_genome_pos",
                    "overlaps_string_pos")
        }
        temp_overlap_table$snp_id <- rep(snp, nrow(temp_overlap_table))

        snp_table <- rbind(snp_table, data.frame(
            snp_id = snp,
            fasta_position = as.character(paste(snp, collapse = ", ")),
            genome_position = genome_pos,
            string_start = string_start,
            string_end = string_end,
            snp_string_pos = as.character(
                paste(snp_string_pos, collapse = ", ")),
            stringsAsFactors = FALSE
        ))
        overlap_table <- rbind(overlap_table, temp_overlap_table)
    }
    return(list(snp_table = snp_table, overlap_table = overlap_table))
}
