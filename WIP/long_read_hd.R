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
run <- (ref_f1$genome_position %/% 10000)+ 1
ref_f1["run"] <- run

hd <- list()
for (i in unique(ref_f1$run)){
  exc <- ref_f1[ref_f1$run != i, "fasta_position"]
  tryCatch({
    hd_res <- find_optimised_snps(f1, excluded_positions = unique(c(exc, all_exc_pos)), max_depth = 5, bp = bp)
    output_result(hd_res, view = "csv", file_name = paste("hd_run", i, ".tsv", sep = ""))
    hd[[as.character(i)]] <- names(tail(hd_res$results[[1]], n = 1))
  }, error = function(err){
    print(err)
    if (is.null(hd[[as.character(i)]])){
      hd[[as.character(i)]] <- -1
    }
  }, warning = function(warning){
    print(warning)
    if (is.null(hd[[as.character(i)]])){
      hd[[as.character(i)]] <- -1
    }
  }, finally={
    print(paste("Completed:", i))
  })
}

df <- data.frame(matrix(unlist(hd), nrow=length(hd), byrow=TRUE), stringsAsFactors=FALSE)
colnames(df) <- c("high_d_sets")
write.csv(df, "all_high_d_sets.csv", row.names = F)


###
#
#  All the SNPs 
#  awk 'FNR == 5' hd_run*.tsv | awk -v FS='\t' '{print $1}' | more
#
###