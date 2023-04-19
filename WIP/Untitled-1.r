devtools::load_all("~/minSNPs")
f1 <- read_fasta("./OUTF1.fasta")
library(BiocParallel)
bp <- MulticoreParam(workers = 32)
res_d1 <- bplapply(1:length(f1[[1]]), cal_met_snp, seqc = f1, metric = calculate_simpson, list(), BPPARAM = bp)


calculate_minor_allele <- function(pattern) {
    
    total <- length(pattern)

    sorted_freq <- sort(unname(table(unlist(pattern))), decreasing=T)
    maf <- 0
    if (length(sorted_freq) > 1) {
        maf <- sorted_freq[2]
    }
    return(list(result = maf/total))
}

res_minor <- bplapply(1:length(f1[[1]]), cal_met_snp, seqc = f1, metric = calculate_minor_allele, list(), BPPARAM = bp)

ref_f1 <- read.csv("ref_OUTF1.csv")
res_df <- data.frame(fasta_pos = 1:length(f1[[1]]), genome_pos = ref_f1$genome_position)
res_df[["simpson"]] <- unlist(res_d1)
res_df[["minor_allele"]] <- unlist(res_minor)

write.csv(res_df, "d1_simpson.csv", row.names = FALSE)

a[["new_minor"]] <- unlist(res_minor)
a[["new_dnorm"]] <- a$simpson/a$new_minor
write.csv(a, "d1_simpson2.csv", row.names = FALSE)

 ln -s ../*_DeepSimu .; ln -s ../*.fasta .; ln -s ../*.csv .; cp ../ARG_Header_build.sh .;cp ../long_read_inf_ARG.R .