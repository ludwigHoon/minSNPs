library(minSNPs)
library(BiocParallel)
library(data.table)

meta <- read.csv("meta.csv")
regions_exc_list <- list()
countries_exc_list <- list()

bp <- MulticoreParam(workers = 64)
pv40 <- read_fasta("pv40-799-snps.fasta", bp = bp)
g_pos <- c(1:length(pv40[[1]]))


u_region <- unique(meta$region)
i <- 1
for (reg in u_region) {
    print(paste0(i, "/", length(u_region)))
    goi <- meta[meta$region == reg, "sample"]
    temp_result <- bplapply(g_pos, function(pos) {
        temp_result <- minSNPs:::check_multistate(pos, pv40[goi])
        return(temp_result)
    }, BPPARAM = MulticoreParam(workers = 64, progress = TRUE))
    names(temp_result) <- as.character(g_pos)
    exc_snp <- names(temp_result[which(temp_result == TRUE)])
    regions_exc_list[[reg]] <- data.frame(region = reg, exclusion = paste(exc_snp, collapse = ", "))
    i <- i + 1
}
result_region <- rbindlist(regions_exc_list)
write.csv(result_region, "result_region.csv", row.names = FALSE)

u_country <- unique(meta$country)
i<- 1
for (reg in u_country) {
    print(paste0(i, "/", length(u_country)))
    goi <- meta[meta$country == reg, "sample"]
    temp_result <- bplapply(g_pos, function(pos) {
        temp_result <- minSNPs:::check_multistate(pos, pv40[goi])
        return(temp_result)
    }, BPPARAM = MulticoreParam(workers = 64, progress = TRUE))
    names(temp_result) <- as.character(g_pos)
    exc_snp <- names(temp_result[which(temp_result == TRUE)])
    countries_exc_list[[reg]] <- data.frame(country = reg, exclusion = paste(exc_snp, collapse = ", "))
    i <- i + 1
}
result_country <- rbindlist(countries_exc_list)
write.csv(result_country, "result_country.csv", row.names = FALSE)

############

# iterative percent-mode
library(minSNPs)
library(BiocParallel)
library(data.table)

meta <- read.csv("meta.csv")
regions_result <- list()
countries_result <- list()

bp <- MulticoreParam(workers = 64)
pv40 <- read_fasta("pv40-799-snps.fasta", bp = bp)
ppv40 <- process_allele(pv40, check_base = FALSE, bp = bp)
u_region <- unique(meta$region)
i <- 1
for (reg in u_region) {
    print(paste0(i, "/", length(u_region)))
    goi <- meta[meta$region == reg, "sample"]
    temp_result <- find_optimised_snps(seqc = ppv40, metric = "percent", number_of_result = 20, max_depth = 5,
        goi = goi, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
    regions_result[[reg]] <- temp_result
    i <- i + 1
}
saveRDS(regions_result, "regions_result.rds")

u_country <- unique(meta$country)
i <- 1
for (coun in u_country) {
    print(paste0(i, "/", length(u_country)))
    goi <- meta[meta$country == coun, "sample"]
    temp_result <- find_optimised_snps(seqc = ppv40, metric = "percent", number_of_result = 20, max_depth = 5,
        goi = goi, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
    countries_result[[coun]] <- temp_result
    i <- i + 1
}
saveRDS(countries_result, "countries_result.rds")


####
library(minSNPs)
library(BiocParallel)
library(data.table)

meta <- read.csv("meta.csv")
regions_result <- list()
countries_result <- list()

bp <- MulticoreParam(workers = 64)
pv40 <- read_fasta("pv40-799-snps.fasta", bp = bp)
ppv40 <- process_allele(pv40, check_base = FALSE, bp = bp)
u_region <- unique(meta$region)
i <- 1
for (reg in u_region) {
    print(paste0(i, "/", length(u_region)))
    goi <- meta[meta$region == reg, "sample"]
    temp_result <- find_optimised_snps(seqc = ppv40, metric = "percent", number_of_result = 10, max_depth = 1,
        goi = goi, accept_multiallelic = FALSE, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
    regions_result[[reg]] <- temp_result
    i <- i + 1
}
saveRDS(regions_result, "regions_result_wo_multi.rds")

u_country <- unique(meta$country)
i <- 1
for (coun in u_country) {
    print(paste0(i, "/", length(u_country)))
    goi <- meta[meta$country == coun, "sample"]
    temp_result <- find_optimised_snps(seqc = ppv40, metric = "percent", number_of_result = 10, max_depth = 1,
        goi = goi, , accept_multiallelic = FALSE, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
    countries_result[[coun]] <- temp_result
    i <- i + 1
}
saveRDS(countries_result, "countries_result_wo_multi.rds")


#############
library(BiocParallel)
library(minSNPs)

bp <- MulticoreParam(workers = 64)

meta <- read.csv("meta.csv")
pv <- read_fasta("pv40-799-snps.fasta", bp = bp)


selected_iso <- sapply(unique(meta$region), function(x){
    return(meta[meta$region == x, "sample"][1])
})
ppv <- process_allele(pv[selected_iso], check_base = FALSE, bp = bp)
exc_region <- as.numeric(readLines("exc_regions.txt"))
high_d_snps <- find_optimised_snps(seqc = ppv, metric = "simpson", number_of_result = 20, max_depth = 5,
    included_positions = c(), excluded_positions = exc_region, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
output_result(high_d_snps, view = "tsv", file_name = "hd_20x5_exc_region.tsv")

selected_iso <- sapply(unique(meta$country), function(x){
    return(meta[meta$country == x, "sample"][1])
})
ppv <- process_allele(pv[selected_iso], check_base = FALSE, bp = bp)
exc_country <- as.numeric(readLines("exc_countries.txt"))
high_d_snps <- find_optimised_snps(seqc = ppv, metric = "simpson", number_of_result = 20, max_depth = 5, 
    included_positions = c(), excluded_positions = exc_country, output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = TRUE))
output_result(high_d_snps, view = "tsv", file_name = "hd_20x5_exc_country.tsv")