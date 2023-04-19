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
#' @param ref_seq the reference genome sequence
#' @param prev number of characters before the SNP
#' @param after number of characters after the SNP
#' @param position_type type of SNPs input, "fasta"
#' (orthologous SNP matrix based) or "genome"
#' (reference genome based); Default to "fasta"
#' @param extend_length whether to extend the search string
#' before and after the SNP and ignore overlapping SNPs
#' @param bp BiocParallel backend to use
#' @return a list containing 2 dataframes,
#' (1) SNP table - snp_id, fasta_position, genome_position,
#' string_start, string_end, snp_string_pos
#' (2) overlap table - snp_id, overlaps_genome_pos,
#' overlaps_fasta_pos, overlaps_string_pos
#' @export
identify_overlaps <- function(selected_snps, position_reference, ref_seq, 
                              prev, after, position_type = "fasta",
                              extend_length = TRUE, bp = MulticoreParam()) {
    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }
    genome_max <- length(ref_seq)
    selected_snps <- unique(selected_snps)
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

    temp_result <- bplapply(seq_len(nrow(snp_genome_pos)),
        function(pos, snp_genome_pos, prev, after,
            position_reference, extend_length){

        snp_table <- data.frame(
            snp_id = c(), fasta_position = c(), genome_position = c(),
            string_start = c(), string_end = c(), snp_string_pos = c(),
            stringsAsFactors = FALSE
        )
        overlap_table <- data.frame(
            snp_id = c(), overlaps_genome_pos = c(), overlaps_fasta_pos = c(),
            overlaps_string_pos = c(), stringsAsFactors = FALSE
        )

        genome_pos <- snp_genome_pos[pos, "genome_position"]
        snp <- snp_genome_pos[pos, "fasta_position"]
        string_start <- genome_pos - prev
        string_end <- genome_pos + after
        snp_string_pos <- prev + 1
        if (string_start < 1) {
            snp_string_pos <- snp_string_pos + string_start - 1
            string_start <- 1
        }
        if (string_end > genome_max) {
            string_end <- genome_max
        }

        n_test <- 0
        # Check overlaps before
        overlaps_before <- position_reference[
            position_reference$genome_position
            >= (genome_pos - (n_test + prev)) &
                position_reference$genome_position
                < genome_pos,
            c("fasta_position", "genome_position")
        ]

        if (extend_length) {
            while (((n_test + prev) - nrow(overlaps_before) < prev) &&
             (string_start - n_test > 1)) {
                n_test <- n_test + (nrow(overlaps_before) - n_test)
                overlaps_before <- position_reference[
                    position_reference$genome_position
                    >= (genome_pos - (n_test + prev)) &
                        position_reference$genome_position
                        < genome_pos, c("fasta_position", "genome_position")
                ]
            }
        }
        string_start <- string_start - n_test
        snp_string_pos <- snp_string_pos + n_test
        overlaps_before$s_pos <-
            (overlaps_before$genome_position - string_start) + 1

        n_test <- 0
        # Check overlaps after
        overlaps_after <- position_reference[
            position_reference$genome_position
            > genome_pos &
                position_reference$genome_position
                <= (genome_pos + (n_test + after)),
            c("fasta_position", "genome_position")
        ]

        if (extend_length) {
            while (((n_test + after) - nrow(overlaps_after) < after) &&
             (string_end + n_test < genome_max)) {
                n_test <- n_test + (nrow(overlaps_after) - n_test)
                overlaps_after <- position_reference[
                    position_reference$genome_position
                    > genome_pos &
                        position_reference$genome_position
                        <= (genome_pos + (n_test + after)),
                    c("fasta_position", "genome_position")
                ]
            }
        }
        string_end <- string_end + n_test
        overlaps_after$s_pos <-
            (snp_string_pos + (overlaps_after$genome_position - genome_pos))

        temp_overlap_table <- rbind(overlaps_before, overlaps_after)
        if (ncol(temp_overlap_table) >= 3) {
            colnames(temp_overlap_table) <-
                c(
                    "overlaps_fasta_pos", "overlaps_genome_pos",
                    "overlaps_string_pos"
                )
        }
        temp_overlap_table$snp_id <- rep(snp, nrow(temp_overlap_table))

        snp_table <- rbind(snp_table, data.frame(
            snp_id = as.character(snp),
            fasta_position = as.character(paste(snp, collapse = ", ")),
            genome_position = as.character(paste(genome_pos, collapse = ", ")),
            string_start = string_start,
            string_end = string_end,
            snp_string_pos = as.character(
                paste(snp_string_pos, collapse = ", ")
            ),
            stringsAsFactors = FALSE
        ))
        overlap_table <- rbind(overlap_table, temp_overlap_table)
        return(list(overlap_table = overlap_table, snp_table = snp_table))
    }, snp_genome_pos = snp_genome_pos, prev = prev, after = after,
            position_reference = position_reference,
            extend_length = extend_length, BPPARAM = bp)

    comb_snp_table <- do.call(rbind,
        bplapply(temp_result, "[[", "snp_table", BPPARAM = bp))
    comb_overlap_table <- do.call(rbind,
        bplapply(temp_result, "[[", "overlap_table", BPPARAM = bp))
    result <- list(snp_table = comb_snp_table,
        overlap_table = comb_overlap_table)
    attr(result, "tables") <- "snp_overlap_tables"
    return(result)
}

#' \code{generate_snp_search_string}
#'
#' @description
#' \code{generate_snp_search_string} generate the search strings
#' for SNPs based on the 2 tables
#' generated from \code{identify_overlaps}
#' @param snp_overlap_tables a list containing the snp_table and overlap_table
#' (can be used instead of snp_table and overlap_table argument)
#' @param snp_table SNP table generated from \code{identify_overlaps}
#' @param overlap_table overlap table generated from \code{identify_overlaps}
#' @param orth_matrix the relevant orthologous SNP matrix containing the SNPs
#' @param ref_seq the sequence of reference genome
#' @param include_neighbour whether to include the neighbouring SNPs
#' in the search string
#' @param unique_only returns only unique search strings
#' @param bp BiocParallel backend to use for parallelization
#' @return a list containing 2 dataframes,
#' (1) string table - string_id (snp_id), search_string, strand,
#' type (snp or kmer);
#' (2) SNP table - snp_id, snp_string, fasta_position,
#' n_match_genome, n_match_genome_rev_com
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
generate_snp_search_string <- function(snp_overlap_tables = NULL,
    snp_table = NULL, overlap_table = NULL, orth_matrix,
    ref_seq, include_neighbour = FALSE, unique_only = TRUE,
    bp = MulticoreParam()) {

    if (is.null(snp_table) && is.null(overlap_table)
        && is.null(snp_overlap_tables)) {
        stop(paste0("Please provide either [snp_table and overlap_table]",
            " or snp_overlap_tables"))
    }

    if (is.null(snp_table) || is.null(overlap_table)) {
        if (is.null(snp_overlap_tables)) {
            stop("snp_table and overlap_table cannot be NULL")
        }
        snp_table <- snp_overlap_tables$snp_table
        overlap_table <- snp_overlap_tables$overlap_table
    }

    string_table <- data.frame(
        search_string = c(), string_id = c(), strand = c(),
        type = c(), stringsAsFactors = FALSE
    )

    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }

    l_string_table <- BiocParallel::bplapply(seq_len(nrow(snp_table)),
        function(i, snp_table, ref_seq, overlap_table, orth_matrix,
            match_count) {
        snp_id <- snp_table[i, "snp_id"]

        ref_string <- ref_seq[
            c(snp_table[i, "string_start"]:snp_table[i, "string_end"])
        ]

        overlaps <- overlap_table[overlap_table$snp_id == snp_id, ]

        fasta_positions <- as.numeric(
            strsplit(snp_table[i, "fasta_position"], split = ", ")[[1]]
        )
        genome_positions <- as.numeric(
            strsplit(snp_table[i, "genome_position"], split = ", ")[[1]]
        )

        strings_pos <- as.numeric(
                    strsplit(snp_table[i, "snp_string_pos"],
                        split = ", ")[[1]]
                )

        if (include_neighbour) {
            fasta_positions <- sort(c(fasta_positions,
                as.numeric(
                    overlaps[overlaps$snp_id == snp_id, "overlaps_fasta_pos"])),
                decreasing = FALSE)
            genome_positions <- sort(c(genome_positions,
                as.numeric(
                    overlaps[overlaps$snp_id == snp_id, "overlaps_genome_pos"])),
                decreasing = FALSE)
            strings_pos <- sort(c(strings_pos,
                as.numeric(
                    overlaps[overlaps$snp_id == snp_id, "overlaps_string_pos"])),
                decreasing = FALSE)
        }

        stopifnot(
            "FASTA positions and string positions are not the same length"
            =
            length(fasta_positions) == length(strings_pos),
            "FASTA positions and genome positions are not the same length"
            =
            length(genome_positions) == length(fasta_positions))

        variants <-
            minSNPs:::generate_pattern(orth_matrix, fasta_positions)
        snp_ids <- names(variants)
        variants <- unlist(variants)
        if (unique_only) {
            variants <- unique(variants)
            snp_ids <- paste0("snp_", snp_id, "_", seq_along(variants))
        } else {
            snp_ids <- paste0("snp_", snp_id, "_", snp_ids)
        }

        f_string <- lapply(variants, function(var, ref_string) {
            ref_string[strings_pos] <- strsplit(var, split = "")[[1]]
            return(paste(ref_string, collapse = ""))
        }, ref_string = ref_string)

        if (! include_neighbour) {
            f_string <- lapply(f_string, function(string) {
                string <- strsplit(string, split = "")[[1]]
                string[overlaps$overlaps_string_pos] <- "."
                return(paste(string, collapse = ""))
            })
        }

        rc_string <- lapply(f_string, reverse_complement)
        temp_string_table <- data.frame(
            search_string = c(unlist(f_string), unlist(rc_string)),
            snp_id = rep(snp_ids, 2),
            strand = c(rep("+", length(f_string)),
                rep("-", length(rc_string))),
            snp_string = rep(variants, 2), stringsAsFactors = FALSE
        )
        temp_string_table$n_match_reference <- unlist(
            lapply(temp_string_table$search_string, match_count,
                search_from = paste(ref_seq, collapse = "")
        ))
        temp_string_table$n_match_reference_rev_com <- unlist(
            lapply(temp_string_table$search_string, match_count,
                search_from = reverse_complement(paste(ref_seq, collapse = ""))
        ))
        temp_string_table$fasta_position <- rep(
            paste(fasta_positions, collapse = ", "), nrow(temp_string_table))
        temp_string_table$genome_position <- rep(
            paste(genome_positions, collapse = ", "), nrow(temp_string_table))

        return(temp_string_table)
    },
    snp_table = snp_table, ref_seq = ref_seq, overlap_table = overlap_table,
    orth_matrix = orth_matrix, match_count = match_count, BPPARAM = bp)

    temp_table <- do.call(rbind, l_string_table)

    string_table <- temp_table[, c("snp_id", "search_string", "strand")]
    string_table$type <- "snp"
    names(string_table)[which(names(string_table) == "snp_id")] <- "string_id"
    snp_table <- temp_table[temp_table$strand == "+",
        c("snp_id", "snp_string", "fasta_position", "genome_position",
        "n_match_reference", "n_match_reference_rev_com")]
    return(list(string_table = string_table, snp_table = snp_table))
}

#' \code{string_table_to_fasta}
#'
#' @description
#' \code{string_table_to_fasta} convert the string table into format that can
#' be written to fasta with \code{write_fasta}
#' @param string_table string table from either:
#'  (1) \code{generate_snp_search_string} or
#'  (2) \code{generate_kmer_search_string}
#' @param strand filter to only include specific strand, default to only "+"
#' @return Will return a list 
string_table_to_fasta <- function(string_table, strand = "+") {
    selected <- string_table[which(string_table$strand %in% strand), ]
    seqs <- as.list(selected$search_string)
    names(seqs) <- selected$string_id
    return(seqs)
}

#' \code{generate_kmer_search_string}
#'
#' @description
#' \code{generate_kmer_search_string} generate the search strings
#' to detect genes' presence
#' @param gene_seq sequences to generate k_mers from
#' @param ref_seq reference genome sequence
#' @param k kmer length
#' @param id_prefix prefix for the gene id
#' @param bp BiocParallel backend to use for parallelization
#' @return a list containing 2 dataframes,
#' (1) string table - string_id (kmer_id), search_string, strand,
#' type (snp or kmer);
#' (2) kmer table - kmer_id, match_gene,
#' n_match_genome, n_match_genome_rev_com
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
generate_kmer_search_string <- function(gene_seq, ref_seq, k,
    id_prefix = "gene", bp = MulticoreParam()) {

    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }

    search_strings <- bplapply(gene_seq, function(gene, generate_kmers, k) {
        kmers <- generate_kmers(final_string = gene, k = k)
        return(kmers)
    }, k = k, generate_kmers = generate_kmers, BPPARAM = bp)

    search_strings <- unique(unlist(search_strings))

    string_ids <- bplapply(seq_along(search_strings),
        function(id, search_strings, id_prefix) {
            return(paste0(id_prefix, "_", id))
        },
    search_strings = search_strings, id_prefix = id_prefix,
    BPPARAM = bp)

    genes_matches <- rep(id_prefix, length(search_strings))

    f_string <- unlist(search_strings)
    rc_string <- unlist(bplapply(f_string, reverse_complement, BPPARAM = bp))

    string_table <- data.frame(
        string_id = rep(unlist(string_ids), 2),
        search_string = c(f_string, rc_string),
        strand = c(rep("+", length(f_string)),
            rep("-", length(rc_string))),
        stringsAsFactors = FALSE
    )
    string_table$type <- "kmer"

    n_match_ref <- unlist(bplapply(f_string, match_count,
        search_from = paste(ref_seq, collapse = "")))

    n_match_ref_rc <- unlist(bplapply(f_string, match_count,
        search_from = reverse_complement(paste(ref_seq, collapse = ""))))

    kmer_table <- data.frame(
        string_id = unlist(string_ids),
        match_gene = genes_matches,
        n_match_reference = n_match_ref,
        n_match_reference_rev_com = n_match_ref_rc,
        stringsAsFactors = FALSE
    )

    return(list(string_table = string_table, kmer_table = kmer_table))
}

#' \code{read_sequences_from_fastq}
#'
#' @description
#' \code{read_sequences_from_fastq} get the sequences
#' from a fastq file, it completely ignores the quality scores
#' @param fastq_file location of the fastq file
#' @param force_to_upper whether to transform sequences
#' to upper case, default to TRUE
#' @param bp BiocParallel backend to use for parallelization
#' @return will return a list of sequences
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
read_sequences_from_fastq <- function(fastq_file, force_to_upper = TRUE,
    bp = MulticoreParam()) {
    lines <- readLines(fastq_file)
    seqs_id <- lines[seq(1, length(lines), 4)]
    seqs <- lines[seq(2, length(lines), 4)]

    seqs_id <- unlist(bplapply(seqs_id, function(id) {
        id <- strsplit(id, split = "")[[1]]
        return(paste0(id[2:length(id)], collapse = ""))
    }, BPPARAM = bp))

    if (force_to_upper) {
        seqs <- unlist(bplapply(seqs, function(seq) {
            return(toupper(seq))
        }, BPPARAM = bp))
    }

    seqs <- as.list(seqs)
    names(seqs) <- seqs_id
    return(seqs)
}

#' \code{search_from_fastq_reads}
#'
#' @description
#' \code{search_from_fastq_reads} identify the matches
#' from a list of search strings
#' @param fastq_file fastq file containing the runs to search from
#' @param search_tables either a single dataframe or a list of dataframes
#' of the string table from either:
#'  (1) \code{generate_snp_search_string} or
#'  (2) \code{generate_kmer_search_string}
#' @param bp BiocParallel backend to use for parallelization (per file)
#' @param bp_per_read BiocParallel backend to use for parallelization (per read)
#' @param output_temp_result whether to output the temporary results
#' @param output_temp_result_dir directory to output the temporary results
#' @return will return a dataframe containing: -
#' file_path, read_id, matched_string_id, strand, match_count
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
search_from_fastq_reads <- function(fastq_file, search_tables,
    output_temp_result = TRUE, temp_result_folder = "./temp_results",
    bp = MulticoreParam(), bp_per_read = SerialParam()) {

    reads <- read_sequences_from_fastq(fastq_file, bp = bp)
    read_ids <- names(reads)

    if (class(search_tables) == "list") {
        search_tables <- do.call(rbind, search_tables)
    }

    stopifnot(
        "search_tables is not data frame" =
            class(search_tables) == "data.frame",
        "search_tables does not contain all the needed columns" =
            all(
                c("string_id", "search_string", "strand", "type")
                %in% colnames(search_tables)
            )
    )

    temp_result <- bplapply(read_ids, function(id, reads, search_tables,
        bp_per_read, output_temp_result, temp_result_folder, fastq_file) {

        u_string <- unique(search_tables$search_string)
        read <- paste(reads[[id]], collapse = "")

        found_string <- bplapply(u_string, function(string, read) {
            return(match_count(string, read))
        }, read = read, BPPARAM = bp_per_read)

        names(found_string) <- u_string
        found_string <- found_string[found_string > 0]
        found_strings <- names(found_string)
        found_strings_count <- unlist(found_string)

        stopifnot("Strings count and found strings are not the same" =
            length(found_strings) == length(found_strings_count))

        if (length(found_strings) == 0) {
            if (output_temp_result) {
                writeLines(paste0("\"search_string\",\"string_id\",\"strand\"",
                    ",\"type\",\"found_string_count\",\"read_id\",",
                    "\"file_path\""),
                    con =
                    paste0(temp_result_folder, "/",
                        gsub(".gz", "", gsub(".fastq", "", fastq_file)),
                        "_", id, ".csv")
                    )
            }
            return(NULL)
        } else {
            temp_result <- search_tables[
                search_tables$search_string %in% found_strings, ]
            merged_result <- merge(temp_result, data.frame(
                found_string = found_strings,
                found_string_count = found_strings_count,
                stringsAsFactors = FALSE
            ), by.x = "search_string", by.y = "found_string")
            merged_result$read_id <- id
            merged_result$file_path <- fastq_file

            if (output_temp_result) {
                write.csv(merged_result, file =
                    paste0(temp_result_folder, "/",
                        gsub(".gz", "", gsub(".fastq", "", fastq_file)),
                            "_", id, ".csv"),
                    row.names = FALSE)
            }
            return(merged_result)
        }
    }, reads = reads, search_tables = search_tables, bp_per_read = bp_per_read,
    output_temp_result = output_temp_result,
    temp_result_folder = temp_result_folder,
    fastq_file = fastq_file, BPPARAM = bp)

    result <- do.call(rbind, temp_result)
    return(result)
}

#' \code{fastq_reads_length}
#'
#' @description
#' \code{fastq_reads_length} identifies the read length in the fastq file
#' @param fastq_file fastq file containing the reads
#' @param bp BiocParallel backend to use for parallelization
#' @importFrom BiocParallel bplapply SerialParam
#' @return a dataframe of the reads and read_length
fastq_reads_length <- function(fastq_file, bp = SerialParam()) {
    reads <- read_sequences_from_fastq(fastq_file, bp = bp)
    result <- bplapply(names(reads), function(read_id, reads) {
        read <- strsplit(reads[[read_id]], split = "")[[1]]
        read_length <- length(read)
        return(data.frame(
            read_id = read_id,
            read_length = read_length,
            stringsAsFactors = FALSE
        ))
    }, reads = reads, BPPARAM = bp)
    result <- do.call(rbind, result)
    return(result)
}

search <- function(ids) {
    results <- lapply(list.files(pattern = "*.fastq")[ids],
        function(f) {
            print(paste0("DOING ", f))
            return(search_from_fastq_reads(f,
                search_tables = list(snp_string, kmer_string),
                output_temp_result = TRUE,
                temp_result_folder = "./temp_results",
                bp = MulticoreParam(workers = 64, progress = TRUE),
                bp_per_read = MulticoreParam(workers = 1, progress = FALSE))
        )}
    )
    results <- do.call(rbind, results)
    return(results)
}

#' \code{gather_temp_result}
#'
#' @description
#' \code{gather_temp_result} will go to the temporary result folder
#' to obtained all the current results
#' @param output_pattern pattern of the files to concatenate into a dataframe
#' @param temp_result_folder the folder containing the temporary results
#' @return Will return a dataframe of matched strings from the temporary results
gather_temp_result <- function(output_pattern,
    temp_result_folder = "./temp_results") {
    temp_result_files <- list.files(temp_result_folder,
        pattern = output_pattern)
    temp_result <- lapply(temp_result_files, function(f){
        read.csv(paste0(temp_result_folder, "/", f))
    })
    temp_result <- do.call(rbind, temp_result)
    return(temp_result)
}

#' \code{collapse_result}
#'
#' @description
#' \code{collapse_result} will send the gathered matched strings
#' to relevant functions for generate result
#' @param results result from \code{gather_temp_result} or
#' \code{search_from_fastq_reads}
#' @param matched_type_functions a list of functions to apply
#' to the matched strings for different type
#' @param arguments list of arguments to pass to different functions
#' for different type of matched strings
#' @param output_match_strand whether to output the strand statistics
#' for each reads
#' @param bp_per_type BiocParallel backend to use for
#' parallelization in each function
#' @param bp_overall BiocParallel backend to use for
#' parallelizing the different functions
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam
#' @importFrom dplyr %>% group_by summarise
#' @return Will return result from the relevant functions
#' for different type of matches
collapse_result <- function(results, matched_type_functions =
    list(snp = analyse_snps_matches, kmer = analyse_kmer_matches),
    arguments = list(snp = list(), kmer = list()),
    output_match_strand = TRUE,
    bp_overall = SerialParam(),
    bp_per_type = MulticoreParam()) {

    matched_types <- unique(results$type)
    if (! all(matched_types %in% names(matched_type_functions))) {
        not_in_matched_type <- matched_types[
            which(!matched_types %in% names(matched_type_functions))]
        warning(paste(not_in_matched_type, "is not in matched_type_functions"))
    }
    if (! all(matched_types %in% names(arguments))) {
        not_in_arguments <- matched_types[
            which(!matched_types %in% names(arguments))]
        warning(paste(not_in_arguments, "is not in arguments"))
    }

    isolate_result <- bplapply(matched_types,
        function(type, results, matched_type_functions,
            bp_per_type, arguments) {

            if (type %in% names(matched_type_functions)) {
                return(matched_type_functions[[type]](
                    results[results$type == type, ],
                    arguments = arguments[[type]],
                    bp = bp_per_type))
            } else {
                return(NA)
            }
        }, results = results, matched_type_functions = matched_type_functions,
    arguments = arguments, bp_per_type = bp_per_type, BPPARAM = bp_overall)
    names(isolate_result) <- matched_types
    if (output_match_strand) {
        strand_stat <- as.data.frame(results %>%
            dplyr::group_by(read_id, #nolint
                strand) %>% #nolint
            dplyr::summarise(count = dplyr::n(), .groups = "rowwise"))
        colnames(strand_stat) <- c("read_id", "strand", "count")
        isolate_result$strand_stat <- strand_stat
    }
    return(isolate_result)
}

#' \code{analyse_kmer_matches}
#'
#' @description
#' \code{analyse_kmer_matches} will analyse the matched kmers
#' @param kmer_matches result from \code{gather_temp_result}
#' restricted to type "kmer"
#' @param arguments a list containing (1) Kmer table
#' @param bp BiocParallel backend to use for parallelization
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% group_by summarise
#' @return Will return the resulting ST
analyse_kmer_matches <- function(kmer_matches, arguments,
    bp = MulticoreParam()) {
    kmer_table <- arguments[["kmer_table"]]
    min_match_per_read <- arguments[["min_match"]]

    kmer_matches_with_ID <- merge( #nolint
        kmer_matches,
        kmer_table[, c("string_id", "match_gene")],
    by.x = "string_id", by.y = "string_id")

    read_kmer_count <- kmer_matches_with_ID %>%
        dplyr::group_by(read_id, match_gene) %>% #nolint
        dplyr::summarise(
            count = dplyr::n_distinct(search_string), .groups = "rowwise") #nolint

    checked_genes <- unique(kmer_table$match_gene)

    gene_results <- bplapply(checked_genes,
        function(gene, kmer_matches, min_match_per_read,
            read_kmer_count) {

        relevant_reads <- read_kmer_count[
            read_kmer_count$match_gene == gene &
            read_kmer_count$count >= min_match_per_read[[gene]], ] %>%
            pull(read_id)

        kmer_id_count <- kmer_matches %>%
            dplyr::filter(match_gene == gene) %>% #nolint
            dplyr::filter(read_id %in% relevant_reads) %>% #nolint
            dplyr::group_by(search_string) %>% #nolint
            dplyr::summarise(
                count = sum(found_string_count), .groups = "keep" #nolint
            )

        if (nrow(kmer_id_count) < 1) {
            return(
                data.frame(gene = gene, total_count = 0,
                    unique_kmer_count = 0,
                    unique_kmer = "",
                    kmer_count = "")
            )
        }

        total_count <- sum(kmer_id_count[, "count"])

        unique_kmer <- kmer_id_count[sort(kmer_id_count$count,
                decreasing = TRUE),
            "search_string"]

        kmer_count <- kmer_id_count[sort(kmer_id_count$count,
                decreasing = TRUE),
            "count"]

        unique_kmer_count <- length(unlist(unique_kmer))
        return(data.frame(gene = gene, total_count = total_count,
            unique_kmer_count = unique_kmer_count,
            unique_kmer = paste0(unlist(unique_kmer), collapse = "_"),
            kmer_count = paste0(unlist(kmer_count), collapse = "_")))
    }, read_kmer_count = read_kmer_count,
        kmer_matches = kmer_matches_with_ID,
        min_match_per_read = min_match_per_read, BPPARAM = bp)

    result <- do.call(rbind, gene_results)
    return(result)
}

#' \code{analyse_kmer_matches_2}
#'
#' @description
#' \code{analyse_kmer_matches_2} will analyse the matched kmers
#' @param kmer_matches result from \code{gather_temp_result}
#' restricted to type "kmer"
#' @param arguments a list containing (1) Kmer table
#' @param bp BiocParallel backend to use for parallelization
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% group_by summarise filter
#' @return Will return the resulting ST
analyse_kmer_matches_2 <- function(kmer_matches, arguments,
    bp = MulticoreParam()) {
    kmer_table <- arguments[["kmer_table"]]
    min_match_per_read <- arguments[["coeffs"]]

    kmer_matches_with_ID <- merge( #nolint
        kmer_matches,
        kmer_table[, c("string_id", "match_gene")],
    by.x = "string_id", by.y = "string_id")

    read_kmer_count <- kmer_matches_with_ID %>%
        dplyr::group_by(read_id, match_gene) %>% #nolint
        dplyr::summarise(
            count = dplyr::n_distinct(search_string), .groups = "rowwise") #nolint

    read_kmer_count <- merge(read_kmer_count,
        unique(kmer_matches[, c("read_id", "read_length")])
    )
    checked_genes <- unique(kmer_table$match_gene)

    gene_results <- bplapply(checked_genes,
        function(gene, kmer_matches, min_match_per_read,
            read_kmer_count) {

        relevant_reads_check <- read_kmer_count[
            read_kmer_count$match_gene == gene, ]

        n_searched_kmer <- 0
        if (gene == "mecA") {
            n_searched_kmer <- 2477
        }
        if (gene == "pvl_p") {
            n_searched_kmer <- 1259
        }
        if (gene == "lukS") {
            n_searched_kmer <- 943
        }
        if (gene == "lukF") {
            n_searched_kmer <- 970
        }
        p_acceptance <- -2.911843e+00 +
            relevant_reads_check$count * 2.283203e-02 +
            relevant_reads_check$read_length * 5.455275e-07 +
            n_searched_kmer * 1.016444e-03
        p_acceptance <- exp(p_acceptance) / (1 + exp(p_acceptance))

        relevant_reads <- relevant_reads_check[
            which(p_acceptance >= .5), "read_id"]

        kmer_id_count <- kmer_matches %>%
            dplyr::filter(match_gene == gene) %>% #nolint
            dplyr::filter(read_id %in% relevant_reads) %>% #nolint
            dplyr::group_by(search_string) %>% #nolint
            dplyr::summarise(
                count = sum(found_string_count), .groups = "keep" #nolint
            )

        if (nrow(kmer_id_count) < 1) {
            return(
                data.frame(gene = gene, total_count = 0,
                    unique_kmer_count = 0,
                    unique_kmer = "",
                    kmer_count = "")
            )
        }

        total_count <- sum(kmer_id_count[, "count"])

        unique_kmer <- kmer_id_count[sort(kmer_id_count$count,
                decreasing = TRUE),
            "search_string"]

        kmer_count <- kmer_id_count[sort(kmer_id_count$count,
                decreasing = TRUE),
            "count"]

        unique_kmer_count <- length(unlist(unique_kmer))
        return(data.frame(gene = gene, total_count = total_count,
            unique_kmer_count = unique_kmer_count,
            unique_kmer = paste0(unlist(unique_kmer), collapse = "_"),
            kmer_count = paste0(unlist(kmer_count), collapse = "_")))
    }, read_kmer_count = read_kmer_count,
        kmer_matches = kmer_matches_with_ID,
        min_match_per_read = min_match_per_read, BPPARAM = bp)

    result <- do.call(rbind, gene_results)
    return(result)
}

#' \code{analyse_snps_matches}
#'
#' @description
#' \code{analyse_snps_matches} will analyse the matched snps
#' @param snp_matches result from \code{gather_temp_result}
#' restricted to type "snp"
#' @param arguments a list containing (1) SNP table, (2) SNP matrix
#' @param bp BiocParallel backend to use for parallelization
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% group_by summarise
#' @return Will return the resulting ST
analyse_snps_matches <- function(snp_matches,
    arguments, bp = MulticoreParam()) {

    snp_table <- arguments[["snp_table"]]
    snp_matrix <- arguments[["snp_matrix"]]

    all_matches <- merge(snp_matches[, c("string_id", "strand",
        "found_string_count", "read_id")],
    snp_table, by.x = "string_id", by.y = "snp_id", all.x = TRUE)

    temp <- all_matches %>% dplyr::group_by(fasta_position, #nolint
        snp_string) %>% #nolint
        dplyr::summarise(count = sum(found_string_count), .groups = "keep") #nolint
    temp_max <- temp %>% dplyr::group_by(fasta_position) %>% #nolint
        dplyr::summarise(max_count = max(count), .groups = "keep") #nolint
    temp_result <- merge(temp, temp_max, by.x = c("fasta_position", "count"),
        by.y = c("fasta_position", "max_count"))

    result <- bplapply(seq_len(nrow(temp_result)),
        function(i, temp_result, snp_matrix) {
            fasta_position <- as.numeric(temp_result[i, "fasta_position"])
            patterns <- minSNPs:::generate_pattern(snp_matrix, fasta_position)
            matching_pat <- names(patterns)[
                which(patterns == temp_result[i, "snp_string"])]
            return(data.frame(
                fasta_position = rep(fasta_position, length(matching_pat)),
                matched = matching_pat)
            )
        },
    temp_result = temp_result, snp_matrix = snp_matrix, BPPARAM = bp)

    result <- do.call(rbind, result)
    return(result)
}