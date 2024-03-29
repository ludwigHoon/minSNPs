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

#' \code{generate_kmer_search_string}
#'
#' @description
#' \code{generate_kmer_search_string} generate the search strings
#' to detect genes' presence
#' @param gene_seq sequences to generate k_mers from
#' @param k kmer length, default to 15
#' @param id_prefix prefix for the gene id
#' @param bp BiocParallel backend to use
#' @return a dataframe containing the search strings
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
generate_kmer_search_string <- function(gene_seq, k = 15,
    id_prefix = NULL, bp = MulticoreParam()) {

    kmers <- generate_kmers(final_string = gene_seq, k = k)
    result <- BiocParallel::bplapply(seq_len(length(kmers)),
        function(ik, kmers, id_prefix) {
            kmer <- kmers[ik]
            res <- list()
            for (strand in c("+", "-")) {
                if (strand == "-") {
                    kmer <- reverse_complement(kmer)
                }
                res[[paste(ik, strand)]] <-
                    data.frame(
                    "id" = paste(id_prefix, ik, strand, sep = "_"),
                    "type" = "KMER",
                    "sequence" = kmer,
                    "strand" = strand,
                    "result" = paste(id_prefix, collapse = ","),
                    "extra" = ""
                )
            }
            return(do.call(rbind, res))
        }, kmers = kmers, id_prefix = id_prefix, BPPARAM = bp
    )
    result <- do.call(rbind, result)
    return(result)
}

#' \code{generate_snp_search_string}
#'
#' @description
#' \code{generate_snp_search_string} identify the SNPs that will overlap
#' the search strings generated from the targeted SNPs
#' @param selected_snps list of targeted SNPs
#' @param position_reference the mapping between
#' reference genome positions and orthologous SNP matrix positions
#' @param ref_seq the reference genome sequence
#' @param prev number of characters before the SNP, default to 7
#' @param after number of characters after the SNP, default to 7
#' @param position_type type of SNPs input, "fasta"
#' (orthologous SNP matrix based) or "genome"
#' (reference genome based); Default to "fasta"
#' @param extend_length whether to extend the search string
#' before and after the SNP and ignore overlapping SNPs
#' @param bp BiocParallel backend to use
#' @return a dataframe containing the search strings
#' @export
generate_snp_search_string <-
    function(selected_snps, position_reference, ref_seq,
            snp_matrix, prev = 7, after = 7, position_type = "fasta",
            extend_length = TRUE, bp = MulticoreParam()) {

    # Take the first sequence of the reference genome
    # if there are multiple sequences
    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }
    if (position_type == "fasta") {
        snp_genome_pos <- position_reference[
            match(selected_snps, position_reference$fasta_position),
            c("fasta_position", "genome_position")
        ]
    } else {
        snp_genome_pos <- position_reference[
            match(selected_snps, position_reference$genome_position),
            c("fasta_position", "genome_position")
        ]
    }
    genome_max <- length(ref_seq)

    result <- BiocParallel::bplapply(seq_len(nrow(snp_genome_pos)),
        function(pos, snp_genome_pos, prev, after,
        position_reference, extend_length, genome_max) {

        snp_pos <- snp_genome_pos[pos, "fasta_position"]
        genome_pos <- snp_genome_pos[pos, "genome_position"]

        # Search string starts from (prev) bases
        # from the SNP in the reference genome
        string_start <- genome_pos - prev
        # Search string ends at (after) bases
        # from the SNP in the reference genome
        string_end <- genome_pos + after
        # The interested SNP is here in the search string
        snp_string_pos <- prev + 1
        # If search string starts at anywhere less than 1, bring to 1
        if (string_start < 1) {
            snp_string_pos <- snp_string_pos + string_start - 1
            string_start <- 1
        }
        # If search string ends longer than reference genome,
        # end at reference genome's end
        if (string_end > genome_max) {
            string_end <- genome_max
        }

        # Identify overlaps
        overlaps <- identify_overlaps(position_reference,
            genome_pos, prev, after)
        snps_in_string <- do.call(rbind, overlaps)[,"genome_position"] - string_start + 1 #nolint
        # Extend search string
        if (extend_length) {
            t_result <- extend_length_fun(overlaps, position_reference, genome_pos,
                prev, after, string_start, string_end, snp_string_pos,
                genome_max)
            string_start <- t_result$string_start
            string_end <- t_result$string_end
            snp_string_pos <- t_result$snp_pos
            snps_in_string <- t_result$snps_in_string
        }
        alleles <- unlist(unique(generate_pattern(snp_matrix, snp_pos)))
        # Iterate for each allele/combination of allele
        single_result <- list()
        ext_info <- sprintf("start: %s|end: %s|snp: %s|bystanders: %s",
            string_start, string_end, snp_string_pos,
                paste(snps_in_string, collapse = ","))
        for (i in seq_len(length(alleles))) {
            allele <- alleles[i]
            genotype_result <- generate_pattern(snp_matrix, snp_pos)
            genotype_result <- names(which(genotype_result == allele))
            for (strand in c("+", "-")) {
                search_seq <- ref_seq[string_start:string_end] #nolint
                # Replace allele
                search_seq[snp_string_pos] <- allele
                # Replace all the bystander SNPs to .
                search_seq[snps_in_string] <- "."
                # Make search string into a string
                search_seq <- paste(search_seq, collapse = "")
                if (strand == "-") {
                    search_seq <- reverse_complement(search_seq)
                }
                single_result[[paste(i, strand)]] <-
                    data.frame(
                        "id" = paste(snp_pos, allele, strand, sep = "_"),
                        "type" = "SNP",
                        "sequence" = search_seq,
                        "strand" = strand,
                        "result" = paste(genotype_result, collapse = ","),
                        "extra" = paste(ext_info, sep = ",", collapse = "|")
                    )
            }
        }
        single_result <- do.call(rbind, single_result)
        return(single_result)
    }, snp_genome_pos = snp_genome_pos, prev = prev, after = after,
    position_reference = position_reference,
    extend_length = extend_length, genome_max = genome_max, BPPARAM = bp)
    result <- do.call(rbind, result)
    return(result)
}

#' \code{identify_overlaps}
#'
#' @description
#' \code{identify_overlaps} identify the bystander SNPs in the search sequences
#' @param position_reference the mapping between
#' reference genome positions and orthologous SNP matrix positions
#' @param genome_position the SNPs in reference genome positions
#' @param prev number of characters before the SNP
#' @param after number of characters after the SNP
#' @return a list of overlaps before and after the SNPs
#' @export
identify_overlaps <- function(
    position_reference, genome_position, prev, after) {

    overlaps_before <- position_reference[
        position_reference$genome_position
        >= (genome_position - prev) &
            position_reference$genome_position
            < genome_position,
        c("fasta_position", "genome_position")
    ]
    overlaps_after <- position_reference[
        position_reference$genome_position
        > genome_position &
            position_reference$genome_position
            <= (genome_position + after),
        c("fasta_position", "genome_position")
    ]
    return(list(before = overlaps_before, after = overlaps_after))
}

#' \code{extend_length_fun}
#'
#' @description
#' \code{extend_length_fun} identify the bystander SNPs in the search sequences
#' and return the new starting and stopping positions of the search sequences
#' @param overlaps the result from \code{identify_overlaps}
#' @param position_reference the mapping between
#' reference genome positions and orthologous SNP matrix positions
#' @param genome_position the SNPs in reference genome positions
#' @param prev number of characters before the SNP
#' @param after number of characters after the SNP
#' @param ori_string_start the original starting position of the search string
#' @param ori_string_end the original ending position of the search string
#' @param genome_max the maximum position in the reference genome
#' @return a dataframe containing the search strings
#' @export
extend_length_fun <- function(overlaps,
    position_reference, genome_position, prev, after,
    ori_string_start, ori_string_end, ori_snp_pos, genome_max) {

        # How many additional bases are added to before/after the search string
        n_before <- 0
        n_after <- 0
        overlaps_before <- overlaps$before
        overlaps_after <- overlaps$after
        new_string_start <- ori_string_start
        new_string_end <- ori_string_end
        new_snp_pos <- ori_snp_pos
        oth_snps_pos <- c()
        repeat_backward <- TRUE

        # While (1) number of bases before the SNPs (n_before + prev)
        # excluding other overlapping SNPs is less than the required
        # number of bases; and
        # (2) the search string can still be expanded backward
        while (((n_before + prev) - nrow(overlaps_before) < prev) &&
            (ori_string_start - n_before > 1)) {
            # Expand backward by the number of overlapping SNPs
            n_before <- n_before + (nrow(overlaps_before) - n_before)
            if (ori_string_start - n_before < 1) {
                n_before <- n_before - (ori_string_start - n_before)
            }
            # Checking for new overlapping SNPs in the expanded string
            overlaps_before <- position_reference[
                position_reference$genome_position
                >= (genome_position - (n_before + prev)) &
                    position_reference$genome_position
                    < genome_position, c("fasta_position", "genome_position")
            ]
            if (ori_string_start - n_before <= 1) {
                repeat_backward <- FALSE
            }
        }
        # The expanded string has a new starting position
        new_string_start <- ori_string_start - n_before
        # The SNP position in the expanded string is now
        new_snp_pos <- ori_snp_pos + n_before
        # Expand search string forward
        # While (1) number of bases after the SNPs (n_after + prev)
        # excluding other overlapping SNPs is less than the required
        # number of bases; (2) the search string could not be expanded backward
        # in the previous step; and
        # (3) the search string can still be expanded forward
        while (
            (
                ((n_before + prev) - nrow(overlaps_before)) +
                ((n_after + after) - nrow(overlaps_after))
                < prev + after
            ) &&
            (ori_string_end + n_after < genome_max)
        ) {
            n_after <- n_after + (nrow(overlaps_after) - n_after) +
                (nrow(overlaps_before) - n_before)
            if (ori_string_end + n_after > genome_max){
                n_after <- n_after - (ori_string_end + n_after - genome_max)
            }
            overlaps_after <- position_reference[
                position_reference$genome_position
                > genome_position &
                    position_reference$genome_position
                    <= (genome_position + (n_after + after)),
                c("fasta_position", "genome_position")
            ]
        }

        # The expanded string has a new ending position
        new_string_end <- new_string_end + n_after

        if (repeat_backward) {
            # While (1) number of bases before the SNPs (n_before + prev)
            # excluding other overlapping SNPs is less than the required
            # number of bases; and
            # (2) the search string can still be expanded backward
            while ((
                ((n_before + prev) - nrow(overlaps_before)) +
                ((n_after + after) - nrow(overlaps_after))
                < prev + after) &&
                (ori_string_start - n_before > 1)) {
                # Expand backward by the number of overlapping SNPs
                n_before <- n_before + (nrow(overlaps_before) - n_before) +
                    (nrow(overlaps_after) - n_after)
                # Checking for new overlapping SNPs in the expanded string
                overlaps_before <- position_reference[
                    position_reference$genome_position
                    >= (genome_position - (n_before + prev)) &
                        position_reference$genome_position
                        < genome_position, c("fasta_position", "genome_position") #nolint
                ]
            }

            # The expanded string has a new starting position
            new_string_start <- ori_string_start - n_before
            # The SNP position in the expanded string is now
            new_snp_pos <- ori_snp_pos + n_before

        }
        # Identify the positions of all the overlapping SNPs
        oth_snps_pos <- do.call(rbind, list(overlaps_before, overlaps_after))[,"genome_position"] - new_string_start + 1 #nolint
        return(list(
            string_start = new_string_start,
            string_end = new_string_end,
            snp_pos = new_snp_pos,
            snps_in_string = oth_snps_pos
        ))
}

#' \code{identify_fasta_match}
#'
#' @description
#' \code{identify_fasta_match} identify the number of matches of
#' the search string in the reference genome
#' @param search_table the result from \code{generate_snp_search_string}
#' or \code{generate_kmer_search_string} or a combined table from the 2
#' @param reference the reference genome
#' @param bp the BiocParallel backend
#' @return a dataframe containing the search strings with the number of match in reference genome
#' @export
identify_fasta_match <- function(search_table, reference, bp) {
    search_table$match_ref_seq <- unlist(
        bplapply(result$sequence, match_count,
            search_from = paste(reference[[1]], collapse = ""),
            BPPARAM = bp)
    )
    return(search_table)
}



############
# Examples #
############
