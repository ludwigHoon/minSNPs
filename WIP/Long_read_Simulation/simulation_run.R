source("inference_common.R")
library(minSNPs)
library(dplyr)
library(BiocParallel)
library(tictoc)
args <- commandArgs(trailingOnly = TRUE)

isolate_file <- read.csv("isolate_sim.tsv", sep = "\t",
    stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
isolate <- isolate_file[as.numeric(args[[1]]), ]
search_type <- c("hd_200", "all")[as.numeric(args[[2]])]
seed <- as.numeric(args[[3]])
search_table <- c("search_hd_200.csv", "search_all.csv")[as.numeric(args[[2]])]

isolate_name <- isolate$isolate
isolate$CC
file_path <- isolate$filename

matrix <- read_fasta("balanced_mat_without_SRR798725_single_final.fasta")
names(matrix) <- unlist(
    lapply(strsplit(names(matrix), split = "_"), function(x) {
        return(x[1])
    }))

print(search_table)
searches <- read.csv(search_table, stringsAsFactors = FALSE)

bp_backend <- BiocParallel::SerialParam()

## echo "sim_id,step,length,elapsed" > time_taken.csv
sim_id <- paste(isolate_name, search_type, seed, sep = "_")


fq_data <- get_data(file_path)
all_reads <- length(fq_data)
sampled <- sample(seq_len(all_reads))
res_file <- paste("./result/", sim_id, ".csv", sep = "")

i <- 1
sampled_set <- c()
temp_result_list <- list()
while (TRUE) {
    sampled_set <- c(sampled_set, sampled[i])
    temp_data <- paste(fq_data[sampled[i]], collapse = "")
    tic()
    re <- matched_snp(temp_data, searches)
    re$isolate <- rep(isolate_name, nrow(re))
    re$read_id <- rep(i, nrow(re))
    if (nrow(re) == 0) {
        re$count <- rep(i, nrow(re))
    }
    temp_result <- re
    temp_result_list[[i]] <- temp_result
    s1_rt <- toc()
    rt_1 <- unname(s1_rt$toc - s1_rt$tic)
    cat(paste(
        paste(sim_id, 1, 1, rt_1, sep = ","),
        "\n", sep = ""), file =
            paste("./time_taken/time_taken_", sim_id, ".csv", sep = ""),
        append = TRUE)

    tic()
    all_result <- do.call(rbind, temp_result_list)
    all_result$snp_id <- as.character(all_result$snp_id)
    all_result$count <- as.numeric(all_result$count)
    transformed_table <- combine_snps(all_result, matrix)
    temp <- infer_most_likely(transformed_table)
    most_likely <- temp$most_likely
    snp_count <- temp$snp_count
    n_candidate <- temp$n_candidate
    min_diff <- temp$min_diff
    n_groups <- temp$group_count
    group_diffs <- temp$group_diff
    cat(paste(paste(isolate_name, length(sampled_set),
            paste("\"", paste(most_likely, collapse = ", "), "\"", sep = ""),
            snp_count, n_candidate, min_diff, n_groups, group_diffs,
        sep = ","), "\n"), file = res_file, append = TRUE)
    s2_rt <- toc()
    rt_2 <- unname(s2_rt$toc - s2_rt$tic)
    cat(paste(
        paste(sim_id, 2, length(sampled_set), rt_2, sep = ","),
        "\n", sep = ""), 
        file = paste("./time_taken/time_taken_", sim_id, ".csv", sep = ""),
        append = TRUE)

    i <- i + 1
    if (i > length(sampled)) {
        break
    }
    if (search_type == "all") {
        if (snp_count >= 300 && min_diff >= 100 && length(most_likely) == 1) {
            break
        }
    }
    if (search_type == "hd_200") {
        if (snp_count >= 30 && min_diff >= 10 && length(most_likely) == 1) {
            break
        }
    }
}
