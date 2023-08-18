library(minSNPs)
library(BiocParallel)
library(data.table)

meta <- read.csv("meta.csv")
pv40 <- readRDS("ppv40.rds")
colnames(meta)[1] <- "isolate"
colnames(meta)[which(colnames(meta) == "country")] <- "target"

bp <- MulticoreParam(workers = 128)
source("mcc_custom_metric.r")
priority <- generate_prioritisation(meta)
res <- find_optimised_snps(pv40,  number_of_result = 5, max_depth = 10, metric = "mcc_multi", bp = bp, meta = meta, target = "target", priority = priority, output_progress = TRUE)
saveRDS(res, "mcc_multi_5x10.rds")

res <- find_optimised_snps(pv40,  number_of_result = 5, max_depth = 10, metric = "simpson_mcc_multi", bp = bp, meta = meta, target = "target", priority = priority, output_progress = TRUE)
saveRDS(res, "simpson_mcc_multi_5x10.rds")

res <- find_optimised_snps(pv40,  number_of_result = 5, max_depth = 10, metric = "simpson", bp = bp, output_progress = TRUE)
saveRDS(res, "simpson_5x10.rds")



srun --ntasks=1 --nodes=1 --cpus-per-task=64 -t 7-00:00:00 -J simpson_5x10 bash -c "eval \"\$(conda shell.bash hook)\"; conda activate r_minsnps; R CMD BATCH --vanilla minsnps_5x10_simpson.R" &
srun --ntasks=1 --nodes=1 --cpus-per-task=64 -t 7-00:00:00 -J mcc_multi_5x10 bash -c "eval \"\$(conda shell.bash hook)\";conda activate r_minsnps; R CMD BATCH --vanilla minsnps_5x10_multi_mcc.R" &
srun --ntasks=1 --nodes=1 --cpus-per-task=64 -t 7-00:00:00 -J mcc_simpson_5x10 bash -c "eval \"\$(conda shell.bash hook)\";conda activate r_minsnps; R CMD BATCH --vanilla minsnps_5x10_simpson_mcc.R" &


library(minSNPs)
library(data.table)
meta <- read.csv("meta.csv")
pv40 <- readRDS("ppv40.rds")
colnames(meta)[1] <- "isolate"
colnames(meta)[which(colnames(meta) == "country")] <- "target"

unsupervised_sets <- c("simpson_5x10.rds", "mcc_multi_5x10.rds", "simpson_mcc_multi_5x10.rds")
priority <- generate_prioritisation(meta)

all_result <- list()
for (mode_sets in unsupervised_sets){
    result_sets <- readRDS(mode_sets)
    mode <- gsub(".rds", "", mode_sets)
    results <- lapply(result_sets$results, function(res){
        snps <- as.numeric(strsplit(names(tail(res, n = 1)), split = ", ")[[1]])
        ori_result = as.numeric(tail(res, n = 1)[[1]])
        simpson = unlist(tail(find_optimised_snps(included_positions=snps, metric ="simpson", seqc = pv40, max_depth = 0, meta = meta, priority = priority, target = "target")$results, n = 1))
        mcc = unlist(tail(find_optimised_snps(included_positions=snps, metric ="mcc_multi", seqc = pv40, max_depth = 0, meta = meta, priority = priority, target = "target")$results, n = 1))
        simpson_mcc = unlist(tail(find_optimised_snps(included_positions=snps, metric ="simpson_mcc_multi", seqc = pv40, max_depth = 0, meta = meta, priority = priority, target = "target")$results, n = 1))
        data.frame(snps = paste(snps, collapse = ", "), ori_result = ori_result, simpson = simpson, mcc = mcc, simpson_mcc = simpson_mcc)
    })
    mode_result <- rbindlist(results)
    mode_result$mode <- mode 
    all_result[[mode]] <- mode_result
}

all_result <- rbindlist(all_result)



supervised_set <- c("mcc_countries_5x5.rds", "percent_countries_5x5.rds")
colnames(meta)[which(colnames(meta) == "target")] <- "country"
countries_result <- list()



supervised_set <- c("coun_percent.txt", "coun_single_mcc.txt")
modes <- c("percent", "mcc_single")
all_result <- list()
u_country <- unique(meta$country)
k <- 1
for (i in c(1:2)){
    mode <- modes[i]
    res <- readLines(supervised_set[i])
    coun_i <- 0
    for (j in res){
        if (startsWith(j, "[")){
            coun_i <- as.numeric(strsplit(gsub("\"", "", gsub("\\[.\\] ", "", j)), split = "/")[[1]][1])
            next
        }
        coun <- u_country[coun_i]
        goi <- meta[meta$country == coun,]$isolate
        snps <- as.numeric(strsplit(gsub("Selected SNPs are:", "", j), split = " ")[[1]])
        snps <- snps[!is.na(snps)]
        mcc_single <- unlist(tail(find_optimised_snps(seqc = pv40, metric = "mcc_single", included_positions = snps, max_depth = 0, goi = goi)$results, n = 1))
        percent <- unlist(tail(find_optimised_snps(seqc = pv40, metric = "percent", included_positions = snps, max_depth = 0, goi = goi)$results, n = 1))
        all_result[[k]] <- data.frame(mode = mode, goi = coun, snps = paste(snps, collapse = ", "), percent = percent, mcc_single = mcc_single)
        k <- k + 1
    }
}
supervised_result <- rbindlist(all_result)


all_result <- list()
for (mode_set in supervised_set){
    result_sets <- readRDS(mode_sets)
    mode <- gsub(".rds", "", mode_sets)

}
for (coun in u_country) {
    print(paste0(i, "/", length(u_country)))
    goi <- meta[meta$country == coun, "sample"]
    temp_result <- find_optimised_snps(seqc = ppv40, metric = "mcc_single",
        number_of_result = 5, max_depth = 5, goi = goi, accept_multiallelic = TRUE,
        output_progress = TRUE, bp = MulticoreParam(workers = 64, progress = FALSE))
    countries_result[[coun]] <- temp_result
    saveRDS(countries_result, "mcc_countries_5x5.rds")
}

rr <- rbindlist(all_result)


result <- list()
for (coun in names(mcc_res)) {
    goi <- meta[meta$country == coun, "sample"]
    full_sets <- rrr <- lapply(mcc_res[[coun]]$results, tail, n = 1)
    snps <- as.character(sapply(full_sets, names))
    result[[coun]] <- rbindlist(lapply(snps, function(res){
        snps <- as.numeric(strsplit(res, split = ", ")[[1]])
        #ori_result = as.numeric(tail(res, n = 1)[[1]])
        percent = unlist(tail(find_optimised_snps(included_positions=snps, metric ="percent", seqc = pv40, max_depth = 0, goi = goi)$results, n = 1))
        mcc = unlist(tail(find_optimised_snps(included_positions=snps, metric ="mcc_single", seqc = pv40, max_depth = 0, goi = goi)$results, n = 1))
        return(data.frame(
            mode = "mcc_single_fixed",
            snps = paste(snps, collapse = ", "),
            percent = percent,
            mcc_single = mcc,
            goi = coun))
    }))
}

for (i in 1:nrow(rere)){
    snp <- as.numeric(strsplit(rere[i,]$snps, split = ", ")[[1]])
    goi <- meta[meta$country == rere[i,]$goi, "sample"]
    mcc_single_2_score[[i]] = unlist(tail(find_optimised_snps(included_positions=snp, metric ="mcc_single", seqc = pv40, max_depth = 0, goi = goi)$results, n = 1))
}





library(data.table)
library(minSNPs)
pv40 <- readRDS("ppv40.rds")
options(stringsAsFactors = FALSE)
all_supervised <- fread("all_supervised_fixed_w2.csv")
all_unsupervised <- fread("all_unsupervised.csv")

#candidate_sets <- head(all_unsupervised[order(-all_unsupervised$simpson_mcc, -all_unsupervised$mcc),], n = 5)
candidate_sets <- all_unsupervised[all_unsupervised$mode == "simpson_mcc_multi_5x10"]
candidate_sets_snps <- lapply(strsplit(candidate_sets$snps, split = ", "), as.numeric)

source("mcc_custom_metric.r")
meta943 <- fread("meta-943.csv")

colnames(meta943)[1] <- "isolate"
colnames(meta943)[which(colnames(meta943) == "country")] <- "target"
meta <- meta943[meta943$isolate %in% names(pv40$seqc)]
priority <- generate_prioritisation(meta[,c("isolate", "target")])

expanded_sets <- list()
set_id <- 1
for (sets in candidate_sets_snps) {
    snps <- sets
    ii <- 0
    while(TRUE){

        mixed_target <- identify_mixed_targets(snps, meta = meta, target = "target", priority = priority, seqc = pv40$seqc)
        
        if (nrow(mixed_target$missed_target_count) < ii + 1){
            print("NO MORE")
            break
        }
        
        target_goi <- tail(mixed_target$missed_target_count[ii+1,], n = 1)$target
        print(paste0("Adding SNPs for: ", target_goi))
        mcc_single_target <- all_supervised[all_supervised$goi == target_goi 
            & all_supervised$mode == "mcc_single_2", ]$snps[(c(set_id:(set_id+4))-1)%%5+1]
        high_mcc_single_snps <- strsplit(mcc_single_target, split = ", ")
        NO_SNPS_FOUND <- FALSE
        for (snp_candidate in high_mcc_single_snps){
            candidate <- as.numeric(snp_candidate)
            FOUND_ONE <- FALSE
            for (asnps in candidate){
                print(paste0("testing ", asnps))
                if (!asnps %in% snps){
                    snps <- c(snps, asnps)
                    FOUND_ONE <- TRUE
                    break
                }
            }
            if (FOUND_ONE){
                break
            }
            NO_SNPS_FOUND <- TRUE
        }
        if (NO_SNPS_FOUND){
            print(paste0("FAILED to find SNPs for:", target_goi))
            ii <- ii + 1
        }
        print(paste0("Current SNPs: ", paste0(snps, collapse = ", ")))
    }
    expanded_sets[[set_id]] <- snps
    set_id <- set_id + 1
}

lapply(expanded_sets, length)
lapply(expanded_sets, identify_mixed_targets, meta = meta, target = "target", priority = priority, seqc = pv40$seqc)

find_optimised_snps(included_positions=snps, metric ="mcc_multi", seqc = pv40, max_depth = 0, meta = meta, target = "country")



library(naivebayes)
TRAINING <- SNP_trans[match(set3, SNP_trans$fasta_pos)][, -c(1,2)]
TRAINING_1 <- cbind(TRAINING, do.call(cbind, strsplit(unlist(pattern), split = "")))
TRAINING_1[,
    (cols) := lapply(.SD, function(x){
        return(
            as.numeric(x==TRAINING_1$REF)*2
        )
    }), .SDcols = cols]
positions <- TRAINING_1[,1][[1]]
TRAINING_1 <- TRAINING_1[,-c(1:3)]
isolates <- colnames(TRAINING_1)
country <- meta[match(colnames(TRAINING_1), meta$isolate),]$country
TRAINING_TABLE <- as.data.table(t(TRAINING_1))
TRAINING_TABLE[, colnames(TRAINING_TABLE) := lapply(.SD, function(x){
    return(factor(x, levels=c(0,1,2)))
}), .SDcols = colnames(TRAINING_TABLE)]
TRAINING_TABLE[, colnames(TRAINING_TABLE)[1:31] := lapply(.SD, as.numeric),
    .SDcols = colnames(TRAINING_TABLE)[1:31]]

TRAINING_TABLE$target <- country

excluded_positions <- 1:pv40$length
tested_set <- FIN_SET
excluded_positions <- which(!excluded_positions %in% as.numeric(tested_set))
res <- find_optimised_snps(search_from = tested_set, seqc = pv40, meta = meta,
    target = "country", max_depth = length(FIN_SET), metric = "mcc_multi",
    bp = bp, output_progress = TRUE)

selected <- as.numeric(strsplit(names(tail(res$results[[1]], n = 1)), split = ", ")[[1]])
tested_set <- tested_set[!tested_set %in% selected]
excluded_positions <- which(!excluded_positions %in% as.numeric(tested_set))



1. Selected SNPs are:  228502 213179 39967 47157 213300 136086 213309 37275 81671 152000 222749 190246 15148 113621 185010 6468 30704 222194 1478 832 41777 92 127391
2. Selected SNPs are:  228521 213069 111525 138164 58701 215503 167346 217836 39098 229040 208179 109747 55296 31294 111902 213376 4559 49585 40031 108910 221881 54946 194967 223125 142106 38689 40181 198744 26688 9595 24050 228056 16880 28295

 [1]   4559  54946 221881  40031  26688   9595 111525 229040 213376 138164
[11] 111902  40181 217836  55296 208179 228521  16880 228056  28295 223125
[21] 142106 108910  31294 194967  24050  38689 213069  58701 215503 198744
[31] 167346  49585 109747  39098

 [1] 228502 213179  39967  47157 213300 136086 213309  37275  81671 152000
[11] 222749 190246  15148 113621 185010   6468  30704 222194   1478    832
[21]  41777     92 127391

tested_set <- FIN_SET
while(length(tested))







library(minSNPs)
library(data.table)
pv40 <- readRDS("ppv40.rds")
library(BiocParallel)
snps_sets <- readLines("all_fin_snps.csv")
sim_gsk <- fread("sim+gsk_test_2.tsv")
snp_trans <- fread("SNPs_pos.csv")
meta943 <- fread("meta-943.csv")

source("mcc_custom_metric.r")
country_region <- unique(meta943[,2:3])
colnames(meta943)[1] <- "isolate"
colnames(meta943)[2] <- "target"
colnames(meta943)[3] <- "region"
snps_sets <- lapply(strsplit(snps_sets, split = ", "), as.numeric)
sim_gsk_pos <- lapply(snps_sets, function(x){
    return(snp_trans[match(x, snp_trans$fasta_pos),]$POS)
})

sets_profiles <- lapply(seq_len(length(sim_gsk_pos)), function(i){
    x <- sim_gsk_pos[[i]]
    xx <- c("SAMPLE", x)
    all_profiles <- sim_gsk[sim_gsk$SAMPLE %in% meta943$isolate,..xx]
    profiles <- lapply(seq_len(nrow(all_profiles)), function(y){paste0(all_profiles[y,-1], collapse="")})
    val_set <- data.frame(set = i, isolate = all_profiles$SAMPLE, profile = unlist(profiles), stringsAsFactors = FALSE)
    return(val_set)
})


profile_targets <- lapply(seq_len(length(snps_sets)), function(x){
    result <- map_profile_to_target(meta943,
        generate_pattern(pv40$seqc, snps_sets[[x]]),
        generate_prioritisation(meta=meta943))
    result$set <- x
    return(result)
    })

all_sets_profiles <- rbindlist(sets_profiles)
all_profile_targets <- rbindlist(profile_targets)

result <- merge(all_sets_profiles, all_profile_targets, by.x = c("set", "profile"), by.y = c("set", "profile"), all.x = TRUE)
colnames(result)[which(colnames(result) == "target")] <- "INFERRED"

result <- merge(result, meta943, by.x = c("isolate"), by.y = c("isolate"), all.x = TRUE)
colnames(result)[which(colnames(result) == "target")] <- "country"
result <- result[,c("set", "isolate", "country", "region", "profile", "INFERRED")]

#colnames(meta943)[2] <- "country"
#colnames(meta943)[3] <- "target"
#profile_targets <- lapply(seq_len(length(snps_sets)), function(x){
#    result <- map_profile_to_target(meta943,
#        generate_pattern(pv40$seqc, snps_sets[[x]]),
#        generate_prioritisation(meta=meta943))
#    result$set <- x
#    return(result)
#    })
#all_profile_targets <- rbindlist(profile_targets)
#result <- merge(result, all_profile_targets, by.x = c("set", "profile"), by.y = c("set", "profile"), all.x = TRUE)
#colnames(result)[which(colnames(result) == "target")] <- "INFERRED-region"

result$`IS_CORRECT` <- result$`INFERRED` == result$country
#result$`IS_CORRECT-region` <- result$`INFERRED-region` == result$region

#for (i in 1:nrow(snp_allele)){
#    pat <- generate_pattern(pv40$seqc, snp_allele[i,]$snp)
#    print(snp_allele[i,]$allele %in% pat)
#}


similarity_score <- function(profile, search_from, return_position = FALSE){
    compare_this <- strsplit(profile, split = "")[[1]]
    with_this <- strsplit(search_from, split = "")
    res <- lapply(with_this, function(x){
        return(length(which(x == compare_this))/length(compare_this))
    })
    names(res) <- search_from
    return(res)
}

infer_from_similarity_to_profile <- function(profile, profile_targets, return_all = FALSE){
    scores <- similarity_score(profile, profile_targets$profile)
    setDT(profile_targets)
    profile_targets_prop <- profile_targets[,.(prop=.N/profile_targets[,.N]), target]
    all_scores <- data.frame(profile = names(scores), similarity = as.vector(unlist(scores)))
    all_scores <- merge(all_scores, profile_targets, by.x = "profile", by.y = "profile", all.x = TRUE)
    setDT(all_scores)
    ordered_scores <- all_scores[, .(target, probability = similarity)][order(-probability)]
    setDF(ordered_scores)
    if (return_all){
        return(ordered_scores)
    }
    return(ordered_scores[ordered_scores$probability == max(ordered_scores$probability),])
}


rr <- list()
for (set in 1:5){
    temp <-
        bplapply(result[eval(result$set == set),]$profile, infer_from_similarity_to_profile,
            profile_targets = all_profile_targets[eval(all_profile_targets$set == set),], BPPARAM = MulticoreParam(workers = 4, progress = TRUE))
    names(temp) <- result[eval(result$set == set),]$isolate
    temp <- rbindlist(temp, idcol = "isolate")
    temp <- temp[, .(result = paste0(paste(target, probability, sep = ":"), collapse = "; ")), isolate]
    temp$set <- set
    rr[[set]] <- temp
}
rr <- rbindlist(rr)
result <- merge(result, rr, by=c("set", "isolate"))
colnames(result)[which(colnames(result) == "result")] <- "INFERRED_2"
result$IS_MULTI_RESULT <- result$IS_MULTI_RESULT <- grepl(";", result$`INFERRED_2`)

result$CORRECT_INFERRED_2 <- sapply(1:nrow(result), function(x){
    return(grepl(result[x,]$country, result[x,]$`INFERRED_2`))
})



snps_to_target_prob <- function(snp, meta, seqc){
    pattern <- generate_pattern(seqc, snp)
    target_count <- as.data.frame(table(meta[meta$isolate %in% names(pattern),]$target))
    colnames(target_count) <- c("target", "count.x")
    allele_prob <- data.frame(isolate = names(pattern), allele = unname(unlist(pattern)))
    allele_prob <- merge(allele_prob, meta, by.x = "isolate", by.y = "isolate", all.x = TRUE)[,c("target", "allele")]
    setDT(allele_prob)
    allele_prob <- allele_prob[,.(count.y=.N), c("allele", "target")]
    allele_prob <- merge(allele_prob, target_count, by = "target", all.x = TRUE)[,.(allele, target, probability = count.y/count.x)]
    allele_prob$snp <- snp
    return(allele_prob)
}

all_snps_prob <- lapply(unique(unlist(snps_sets)), snps_to_target_prob, meta = meta943, seqc = pv40$seqc)
all_snps_prob <- rbindlist(all_snps_prob)

profile_to_target_prob <- function(profile, snp_prob, snps, return_all = FALSE){
    snp_allele <- data.frame(snp = snps, allele = strsplit(profile, split = "")[[1]])
    snp_allele <- merge(snp_allele, snp_prob, by.x=c("snp", "allele"), by.y = c("snp", "allele"))
    setDT(snp_allele)
    n_snps <- length(unique(snp_allele$snp))
    prob_t <- snp_allele[,.(probability = mean(probability)), c("target")]
    if (return_all){
        return(prob_t)
    }
    max_prob <- max(prob_t$probability)
    return(prob_t[probability == max_prob,])
}

rr <- list()
for (set in 1:5){
    temp <-
        bplapply(result[eval(result$set == set),]$profile, profile_to_target_prob,
               snp_prob = all_snps_prob, snps = snps_sets[[set]], BPPARAM = MulticoreParam(workers = 4, progress = TRUE))
    names(temp) <- result[eval(result$set == set),]$isolate
    temp <- rbindlist(temp, idcol = "isolate")
    temp <- temp[, .(result = paste0(paste(target, probability, sep = ":"), collapse = "; ")), isolate]
    temp$set <- set
    rr[[set]] <- temp
}
rr <- rbindlist(rr)
result <- merge(result, rr, by=c("set", "isolate"))
colnames(result)[which(colnames(result) == "result")] <- "INFERRED_3"
result$IS_MULTI_RESULT_3 <- grepl(";", result$`INFERRED_3`)

result$CORRECT_INFERRED_3 <- sapply(1:nrow(result), function(x){
    return(grepl(result[x,]$country, result[x,]$`INFERRED_3`))
})

infer_combined <- function(profile, snp_prob, profile_targets, snps, use_sum = TRUE, return_all=FALSE) {
    res1 <- profile_to_target_prob(profile, snp_prob, snps, return_all = TRUE)
    colnames(res1)[2] <- "probability.x"
    res2 <- infer_from_similarity_to_profile(profile, profile_targets, return_all = TRUE)
    colnames(res2)[2] <- "probability.y"
    res <- merge(res2, res1, by = "target")

    if (use_sum) {
        res$probability <- res$probability.x + res$probability.y
    } else {
        res$probability <- res$probability.x * res$probability.y
    }
    if (return_all){
        return(res)
    }
    rr <- res[res$probability == max(res$probability),]
    setDT(rr)
    
    return(
        rr[, .(result = paste0(paste(target, probability, sep = ":"),
            collapse = "; "))]
    )
}

rr <- list()
for (set in 1:5){
    temp <-
        bplapply(result[eval(result$set == set),]$profile, infer_combined,
               snp_prob = all_snps_prob,
               profile_targets = all_profile_targets[
                    eval(all_profile_targets$set == set),],
                snps = snps_sets[[set]], BPPARAM = MulticoreParam(workers = 4, progress = TRUE))
    names(temp) <- result[eval(result$set == set),]$isolate
    temp <- rbindlist(temp, idcol = "isolate")
    temp$set <- set
    rr[[set]] <- temp
}
rr <- rbindlist(rr)

result <- merge(result, rr, by=c("set", "isolate"))
colnames(result)[which(colnames(result) == "result")] <- "INFERRED_combined_sum"
result$IS_MULTI_combined_sum <- grepl(";", result$`INFERRED_combined_sum`)

result$CORRECT_INFERRED_combined_sum <- sapply(1:nrow(result), function(x){
    return(grepl(result[x,]$country, result[x,]$`INFERRED_combined_sum`))
})


rr <- list()
for (set in 1:5){
    temp <-
        bplapply(result[eval(result$set == set),]$profile, infer_combined,
               snp_prob = all_snps_prob,
               profile_targets = all_profile_targets[
                    eval(all_profile_targets$set == set),],
                use_sum = FALSE,
                snps = snps_sets[[set]], BPPARAM = MulticoreParam(workers = 4, progress = TRUE))
    names(temp) <- result[eval(result$set == set),]$isolate
    temp <- rbindlist(temp, idcol = "isolate")
    temp$set <- set
    rr[[set]] <- temp
}
rr <- rbindlist(rr)

result <- merge(result, rr, by=c("set", "isolate"))
colnames(result)[which(colnames(result) == "result")] <- "INFERRED_combined_mul"
result$IS_MULTI_combined_mul <- grepl(";", result$`INFERRED_combined_mul`)

result$CORRECT_INFERRED_combined_mul <- sapply(1:nrow(result), function(x){
    return(grepl(result[x,]$country, result[x,]$`INFERRED_combined_mul`))
})

fwrite(result, "ALL_RESULT.csv", row.names = FALSE)

infer_from_similarity_to_profile(result[eval(result$set == set),]$profile[1],
            profile_targets = all_profile_targets[eval(all_profile_targets$set == set),])