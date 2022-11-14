library(dplyr)
library(minSNPs)
library(BiocParallel)

source("search_string_generation.r")

kmer_table <- read.csv("kmer_table.csv")
snp_table <- read.csv("snp_table.csv")
read_length <- read.csv("number_bases_reads.csv")

snp_matrix <- read_fasta("balanced_mat_without_SRR798725_single_final.fasta")
names(snp_matrix) <- gsub("_1$", "", names(snp_matrix))

kmer_table[grepl("pvl", kmer_table$match_gene), "match_gene"] <- "pvl_p"
kmer_table[grepl("mecA", kmer_table$match_gene), "match_gene"] <- "mecA"
kmer_table[grepl("lukS", kmer_table$string_id), "match_gene"] <- "lukS"

arguments <- list(kmer =
    list(kmer_table = kmer_table,
        min_match = list(
            mecA = 155,
            pvl_p = 315,
            lukS = 30,
            lukF = 70
        )),
    snp = list(snp_table = snp_table, snp_matrix = snp_matrix))

all_sim_out <- paste0(
    gsub(".fastq", "", list.files(pattern = "*fastq")),
    ".*.csv")

n_bases_all <- read.csv("n_bases_all.csv", stringsAsFactors = FALSE)

sim <- function(sim_iso, seed) {
    fname <- gsub("\\.\\*\\.csv", "", sim_iso)
    set.seed(seed)
    sim <- list.files(path = "./temp_results", pattern = sim_iso)
    reordered <- sample(sim, length(sim))
    if (length(reordered) < 500) {
        return(NA)
    }
    sampling <- c(seq(10, 190, 10),
        seq(200, 490, 50),
        seq(500, length(reordered), 100))
    if (tail(sampling, n = 1) != length(reordered)) {
        sampling <- c(sampling, length(reordered))
    }
    sim_result <- bplapply(sampling, function(i){
        total_read_bases <- sum(as.numeric(
            unlist(
                lapply(
                strsplit(
                    gsub(".csv", "", reordered[1:i]), split = "_"),
                `[`, 4)
            )
        ))
        t_results <- lapply(reordered[1:i], function(f){
            read.csv(paste0("./temp_results", "/", f))
        })
        t_results <- do.call(rbind, t_results)
        t_results <- merge(t_results, n_bases_all)
        result <- collapse_result(t_results, arguments = arguments,
            bp_per_type = SerialParam())

        sim_result <- data.frame(file = fname,
            n_read = i,
            n_bases = total_read_bases)
        CC_count <- table(result[["snp"]][["matched"]])
        T5_CC <- head(CC_count[order(CC_count, decreasing = TRUE)], n = 5)
        #sim_result$T5_ST <- paste(names(T5_ST), collapse = "_")
        #sim_result$T5_ST_match_count <- paste(T5_ST, collapse = "_")
        for (i in 1:5) {
            sim_result[paste0("T", i, "_CC")] <- ifelse(
                is.null(names(T5_CC)[i]),
                NA,
                names(T5_CC)[i])
            sim_result[paste0("T", i, "_CC_match_count")] <- T5_CC[i]
        }
        kmer_res <- result[["kmer"]]
        for (gene in unique(kmer_table$match_gene)){
            if (is.null(kmer_res[kmer_res$gene == gene, ])) {
                sim_result[[paste(gene, "total", sep = "_")]] <- 0
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <- 0
            } else {
                sim_result[[paste(gene, "total", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene, "total_count"],
                            integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "total_count"])
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene,
                            "unique_kmer_count"], integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "unique_kmer_count"])
            }
        }
        return(sim_result)
    }, BPPARAM = MulticoreParam(workers = 62, progress = TRUE))
    sim_result <- do.call(rbind, sim_result)
    return(sim_result)
}

sim_2 <- function(sim_iso, seed){
    fname <- gsub("\\.\\*\\.csv", "", sim_iso)
    set.seed(seed)
    sim <- list.files(path = "./temp_results", pattern = sim_iso)
    reordered <- sample(sim, length(sim))
    if (length(reordered) < 500) {
        return(NA)
    }
    sampling <- c(seq(10, 190, 10),
        seq(200, 490, 50),
        seq(500, length(reordered), 100))
    if (tail(sampling, n = 1) != length(reordered)) {
        sampling <- c(sampling, length(reordered))
    }
    sim_result <- bplapply(sampling, function(i){
        total_read_bases <- sum(as.numeric(
            unlist(
                lapply(
                strsplit(
                    gsub(".csv", "", reordered[1:i]), split = "_"),
                `[`, 4)
            )
        ))
        t_results <- lapply(reordered[1:i], function(f){
            read.csv(paste0("./temp_results", "/", f))
        })
        t_results <- do.call(rbind, t_results)
        t_results <- merge(t_results, n_bases_all)
        result <- collapse_result(t_results, matched_type_functions =
            list(snp = analyse_snps_matches, kmer = analyse_kmer_matches_2),
            arguments = arguments,
            bp_per_type = SerialParam())

        sim_result <- data.frame(file = fname,
            n_read = i,
            n_bases = total_read_bases)
        CC_count <- table(result[["snp"]][["matched"]])
        T5_CC <- head(CC_count[order(CC_count, decreasing = TRUE)], n = 5)
        #sim_result$T5_ST <- paste(names(T5_ST), collapse = "_")
        #sim_result$T5_ST_match_count <- paste(T5_ST, collapse = "_")
        for (i in 1:5) {
            sim_result[paste0("T", i, "_CC")] <- ifelse(
                is.null(names(T5_CC)[i]),
                NA,
                names(T5_CC)[i])
            sim_result[paste0("T", i, "_CC_match_count")] <- T5_CC[i]
        }
        kmer_res <- result[["kmer"]]
        for (gene in unique(kmer_table$match_gene)){
            if (is.null(kmer_res[kmer_res$gene == gene, ])) {
                sim_result[[paste(gene, "total", sep = "_")]] <- 0
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <- 0
            } else {
                sim_result[[paste(gene, "total", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene, "total_count"],
                            integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "total_count"])
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene,
                            "unique_kmer_count"], integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "unique_kmer_count"])
            }
        }
        return(sim_result)
    }, BPPARAM = MulticoreParam(workers = 62, progress = TRUE))
    sim_result <- do.call(rbind, sim_result)
    return(sim_result)
}

sim_kmer_only <- function(sim_iso, seed, arguments_list) {
    fname <- gsub("\\.\\*\\.csv", "", sim_iso)
    set.seed(seed)
    sim <- list.files(path = "./temp_results", pattern = sim_iso)
    reordered <- sample(sim, length(sim))
    if (length(reordered) < 500) {
        return(NA)
    }
    sampling <- c(seq(10, 190, 10),
        seq(200, 490, 50),
        seq(500, length(reordered), 100))
    if (tail(sampling, n = 1) != length(reordered)) {
        sampling <- c(sampling, length(reordered))
    }

    sim_result <- bplapply(sampling, function(i) {
        t_results <- lapply(reordered[1:i], function(f){
            read.csv(paste0("./temp_results", "/", f))
        })
        t_results <- do.call(rbind, t_results)
        result <- collapse_result(t_results,
            matched_type_functions = list(kmer = analyse_kmer_matches),
            arguments = arguments_list,
            bp_per_type = SerialParam())

        sim_result <- data.frame(file = fname,
            n_read = i)
        kmer_res <- result[["kmer"]]
        for (gene in names(arguments_list$kmer$min_match)){
            sim_result[[paste(gene, "min_match", sep = "_")]] <-
                arguments_list$kmer$min_match
            if (is.null(kmer_res[kmer_res$gene == gene, ])) {
                sim_result[[paste(gene, "total", sep = "_")]] <- 0
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <- 0
            } else {
                sim_result[[paste(gene, "total", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene, "total_count"],
                            integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "total_count"])
                sim_result[[paste(gene, "n_unique_kmer", sep = "_")]] <-
                    ifelse(
                    identical(kmer_res[kmer_res$gene == gene,
                            "unique_kmer_count"], integer(0)),
                        0, kmer_res[kmer_res$gene == gene, "unique_kmer_count"])
            }
        }
        return(sim_result)
    }, BPPARAM = MulticoreParam(workers = 62, progress = TRUE))
    sim_result <- do.call(rbind, sim_result)
    return(sim_result)
}