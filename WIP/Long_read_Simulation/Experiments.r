# This is used to extract the sequences for mecA, lukS, lukF gene
# from the assembled genomes of JKD6159 and JRA307.
# The positions are pre-extraced with:
# `grep -A 5 -B 5 "lukS" *.gb >> extracted_relevant_gene_features.txt` #nolint
# `grep -A 5 -B 5 "lukF" *.gb >> extracted_relevant_gene_features.txt` #nolint
# `grep -A 5 -B 5 "mecA" *.gb >> extracted_relevant_gene_features.txt` #nolint
# stored in extracted_relevant_gene_features.txt

library(minSNPs)

extract_sequences_from_genome <-
    function(start, end, genome, rev, name) {
        sequence <- paste(genome[start:end], collapse = "")
        if (rev) {
            sequence <- reverse_complement(sequence)
        }
        return(data.frame(name = name, sequence = sequence))
    }

jkd6159 <- read_fasta("JKD6159.fasta")[[1]]
jra307 <- read_fasta("JRA307.fasta")[[1]]

start <- c(1520057, 1519078, 39177, 955935, 39295)
end <- c(1520995, 1520061, 41186, 956654, 41301)
iso <- list(jkd6159, jkd6159, jkd6159, jkd6159, jra307)
rev <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
name <- c("lukS", "lukF", "mecA_l", "mecA_s", "mecA")

all_sequences <- mapply(extract_relevant_sequences_from_genome, 
    start, end, iso, rev, name
)
result <- data.frame(
    name = unlist(all_sequences[1, ]),
    sequences = unlist(all_sequences[2, ]))

write.csv(result, "gene_sequences.csv", row.names = FALSE)

## T
##
#ID, Type, Sequence, Strand, Result, Additional Information
#"FASTA_POS_ALLELE_STRAND", SNP, +, CC1|CC2|, "source: | SNP FASTA Position: | SNP Genome Position: | SNP String position: | Allele: ""
#"GENE_SEQ_STRAND", "KMER", -	mecA	source: 


### These utilities may not be needed in the future if minSNPs output is changed in the upcoming version
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

#### Generate the search strings for SNPs
snps_loc <- process_result_file("sup_5_6_mega_hd_2_steps.txt")
position_reference <- read.csv("ref_OUTF1.csv")
snps <- unique(unlist(snps_loc))
ref_seq <- read_fasta("Mu50.fasta")
snp_matrix <- read_fasta("balanced_mat_without_SRR798725_single_final.fasta")
source("search_string.R")
library(BiocParallel)
bp <- MulticoreParam(workers = 64, progress = TRUE)
result <- generate_snp_search_string(snps, position_reference, ref_seq,
            snp_matrix, 7, 7, position_type = "fasta",
            extend_length = TRUE, bp)
row.names(result) <- seq_len(nrow(result))
result$match_ref_seq <- unlist(
    bplapply(result$sequence, match_count,
    search_from = paste(ref_seq[[1]], collapse = ""), BPPARAM = bp))

#### Generate the search strings for KMER
result2 <- list()
genes <- read.csv("gene_sequences.csv")
for (i in seq_len(nrow(genes))) {
    result2[[i]] <- generate_kmer_search_string(
        genes$sequences[i], 15, genes$name[i], bp)
}
result2 <- do.call(rbind, result2)
row.names(result2) <- seq_len(nrow(result2))
result2$match_ref_seq <- unlist(
    bplapply(result2$sequence, match_count,
    search_from = paste(ref_seq[[1]], collapse = ""), BPPARAM = bp))

### Combine everything together
fin_result <- do.call(rbind, list(result, result2))
write.csv(fin_result, "all_search_string.csv", row.names = FALSE)
