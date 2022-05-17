devtools::load_all("/home/lhoon/minSNPs")
library(BiocParallel)

bp <- MulticoreParam(workers = 16)


processFile <- function(variant_sites){
        result <- c()
        read_cc <- c()
        con <- file(variant_sites, "r")
        while(TRUE){
                data <- readLines(con, n = 3)
                if (length(data) == 0){
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

all_exc_pos <- processFile("variant_t10.csv")

f1 <- read_fasta("./balanced_mat_without_SRR798725_single_final.fasta")
ref_f1 <- read.csv("./ref_OUTF1.csv", stringsAsFactors=F)

windows <- list()
prev_results <- list()
start <- 1
end <- max(ref_f1$genome_position) + 1
while (start <= max(ref_f1$genome_position)){
  inc <- na.omit(ref_f1[(ref_f1$genome_position >= start & ref_f1$genome_position < start + 100), "fasta_position"])
  inc <- inc[! inc %in% all_exc_pos]
  if (! is.null(prev_results[[paste(inc, collapse = ", ")]])){
    windows[[as.character(start)]] <- prev_results[[paste(inc, collapse = ", ")]]
    print(paste("Completed:", start/end))
    start <- start + 1
    next
  }
  if (length(inc) > 0){
    tryCatch({
      hd_res <- find_optimised_snps(f1, included_positions = inc, max_depth = 0, bp = bp)
      prev_results[[paste(inc, collapse = ", ")]] <- tail(hd_res$results[[1]])[[1]]
      windows[[as.character(start)]] <- tail(hd_res$results[[1]])[[1]]
    }, error = function(err){
      print(err)
    }, warning = function(warning){
      print(warning)
    })
  }else{
    windows[[as.character(start)]] <- 0
  }
  print(paste("Completed:", start/end))
  start <- start + 1
}

df <- data.frame(matrix(unlist(windows), nrow=length(windows), byrow=TRUE), stringsAsFactors=FALSE)
colnames(df) <- c("d_val")
df[["window_id"]] <- row.names(df)
row.names(df) <- NULL
write.csv(df[c("window_id", "d_val")], "all_high_d_sets.csv", row.names = F)