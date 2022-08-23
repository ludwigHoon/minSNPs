library(BiocParallel)

hd_200_result_files <- Sys.glob("result/*_hd_200_*.csv")
all_snp_result_files <- Sys.glob("result/*_all_*.csv")

isolate_file <- read.csv("isolate_sim.tsv", sep = "\t",
    stringsAsFactors = FALSE)

get_if_complete <- function(result_file, min_snp, min_diff) {
    c_pipe <- pipe(
        paste("tail -n 1", result_file)
    )
    data <- read.table(text = readLines(
        c_pipe
    ), sep = ",", header = FALSE, col.names = c("isolate", "sampled",
        "most_likely", "snp_conf", "candidate", "min_diff", "group_snp",
        "group_diff"))
    close(c_pipe)
    if (data$snp_conf >= min_snp && data$min_diff >= min_diff) {
        sim_id <- gsub(".csv", "", gsub("result/", "", result_file))
        t_taken <- read.csv(
            paste("time_taken/time_taken_", sim_id, ".csv", sep = ""),
        header = FALSE)
        total_time <- sum(t_taken$V4)
        is_right <- (isolate_file[isolate_file$isolate == data$isolate, "CC"]
            == data$most_likely)
        result <- data.frame(sim_id = sim_id, most_likely = data$most_likely,
            snp_conf = data$snp_conf, min_diff = data$min_diff,
            is_right = is_right, time_taken = total_time)
        return(result)
    } else {
        return(NULL)
    }
}

hd_200_res <- bplapply(hd_200_result_files, get_if_complete,
    min_snp = 30, min_diff = 10, BPPARAM = MulticoreParam(workers = 4))
hd_200_comb <- do.call(rbind, hd_200_res)
all_snp_res <- bplapply(all_snp_result_files, get_if_complete,
    min_snp = 300, min_diff = 100, BPPARAM = MulticoreParam(workers = 4))
all_snp_comb <- do.call(rbind, all_snp_res)

all_res <- rbind(all_snp_comb, hd_200_comb)
write.csv(all_res, "simulated_run.csv", row.names = FALSE)