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
            if (gsub("\\s+", "", preparse[[i + 1]]) == "") {
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
#' (1) SNP table - snp_id, fasta_position, genome_position,
#' string_start, string_end, snp_string_pos
#' (2) overlap table - snp_id, overlaps_genome_pos,
#' overlaps_fasta_pos, overlaps_string_pos
#' @export
identify_overlaps <- function(selected_snps, position_reference, prev, after,
    position_type = "fasta", extend_length = TRUE) {
    snp_table <- data.frame(
        snp_id = c(), fasta_position = c(), genome_position = c(),
        string_start = c(), string_end = c(), snp_string_pos = c(),
        stringsAsFactors = FALSE
    )
    overlap_table <- data.frame(
        snp_id = c(), overlaps_genome_pos = c(), overlaps_fasta_pos = c(),
        overlaps_string_pos = c(), stringsAsFactors = FALSE
    )
    
    selected_snps <- unique(selected_snps)
    if (position_type == "fasta") {
        snp_genome_pos <- position_reference[
            match(selected_snps, position_reference$fasta_position),
                c("fasta_position", "genome_position")]
    } else{
        snp_genome_pos <- position_reference[
            match(selected_snps, position_reference$genome_position),
                c("fasta_position", "genome_position")]
    }

    for (pos in seq_len(nrow(snp_genome_pos))) {
        genome_pos <- snp_genome_pos[pos, "genome_position"]
        snp <- snp_genome_pos[pos, "fasta_position"]
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
            snp_id = as.character(snp),
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

#' \code{generate_search_string}
#'
#' @description
#' \code{generate_search_string} generate the search string based on the 2 tables
#' generated from \code{identify_overlaps}
#' @param snp_table SNP table generated from \code{identify_overlaps}
#' @param overlap_table overlap table generated from \code{identify_overlaps}
#' @param orth_matrix the relevant orthologous SNP matrix containing the SNPs
#' @param ref_seq the sequence of reference genome
#' @param include_neighbour whether to include the neighbouring SNPs in the search string
#' @param bp BiocParallel backend to use for parallelization
#' @return a list containing 2 dataframes,
#' (1) string table - string_id (snp_id), search_string, strand, type (snp or gene);
#' (2) SNP table - snp_id, snp_string, n_match_genome, n_match_genome_rev_com
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
generate_search_string <- function(snp_table, overlap_table, orth_matrix, ref_seq, include_neighbour = FALSE,
    bp = MulticoreParam()){

    string_table <- data.frame(search_string = c(), string_id = c(), strand = c(),
        type = c(), stringsAsFactors = FALSE)

    snp_table <- data.frame(
        snp_id = c(), snp_string = c(), n_match_genome = c(),
        n_match_genome_rev_com = c(), stringsAsFactors = FALSE
    )

    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }

    l_string_table <- BiocParallel::bplapply(seq_len(nrow(snp_table)), function(i){
        snp_id <- snp_table[i, "snp_id"]

        ref_string <- ref_seq[
            c(snp_table[i, "string_start"]:snp_table[i, "string_end"])]

        fasta_positions <- as.numeric(
            strsplit(snp_table[i, "fasta_position"], split = ", ")[[1]])
        overlaps <- overlap_table[overlap_table$snp_id == snp_id,]

        if (include_neighbour){
            ##### TO-DO
            next
        } else {
            variants <- unlist(unique(minSNPs:::generate_pattern(orth_matrix, fasta_positions)))
            f_string <- lapply(variants, function(var){
                ref_string[as.numeric(strsplit(snp_table[i, "snp_string_pos"], split = ", ")[[1]])] <- var
                return(paste(ref_string, collapse = ""))
            })
            f_string <- lapply(f_string, function(string){
                string <- strsplit(string, split = "")[[1]]
                string[overlaps$overlaps_string_pos] <- "."
                return(paste(string, collapse = ""))
            })
        }
        
        rc_string <- lapply(f_string, reverse_complement)
        temp_string_table <- data.frame(search_string = c(unlist(f_string), unlist(rc_string)),
            snp_id = rep(snp_id, length(f_string)*2), snp_string = c(variants, variants), stringsAsFactors = F)
        temp_string_table$n_match_genome <- unlist(lapply(temp_string_table$search_string, match_count,
            search_from = paste(ref_seq, collapse = "")))
        temp_string_table$n_match_genome_rev <- unlist(lapply(temp_string_table$search_string, match_count,
            search_from = reverse_complement(paste(ref_seq, collapse = ""))))
        temp_string_table$fasta_position <- rep(snp_table[i, "fasta_position"], nrow(temp_string_table))
        #toc()
        return(temp_string_table)
    }, BPPARAM = bp)
    string_table <- do.call(rbind, l_string_table)

    return(string_table)
}