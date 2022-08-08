# calculate difference between 2 strings
cal_dis <- function(s1, s2){
    temp <- strsplit(c(s1, s2), split = "")
    s_length <- unique(sapply(temp, length))
    stopifnot(length(s_length) == 1)
    diff <- 0
    for (n in seq_len(s_length[1])){
        if (temp[[1]][n] != temp[[2]][n]){
            diff = diff + 1
        }
    }
    return(diff)
}

# infer result
infer_result <- function(data, searches, seqs){
    s_strings <- unique(searches$search_string)
    temp_result <- data.frame(fragment_id = c(), snp_id = c(),
        extracted_snps = c(), match_count = c(),
        stringsAsFactors = F)

    for (i in seq_len(length(data))) {
        c_count <- BiocParallel::bplapply(s_strings, function(string) {
            match <- unlist(gregexpr(string, paste(data[i], collapse = "")))
            return(match[which(match > 0)])
        }, BPPARAM = BiocParallel::MulticoreParam(workers = 62))
        s_m <- lapply(c_count, length)
        matched_string <- s_strings[which(s_m > 0)]
        match_count <- unlist(s_m[which(s_m > 0)])
        match_indexes <- unlist(lapply(c_count[which(s_m > 0)], paste, collapse = ", "))
        snp_id <- searches[match(matched_string, searches$search_string), "fasta_position"]
        target_pos <- searches[match(matched_string, searches$search_string), "target_position"]
        ext_snps <- searches[match(matched_string, searches$search_string), "snp_string"]
        stopifnot(length(matched_string) == length(match_count))
        stopifnot(length(matched_string) == length(match_indexes))
        new_data <- data.frame(fragment_id = rep(i, length(matched_string)), snp_id = snp_id,
            extracted_snps = ext_snps, match_count = match_count,
            stringsAsFactors = F)
        temp_result <- rbind(temp_result, new_data)
    }

    result_to_use <- data.frame(fragment_id = c(), snp_id = c(),
        extracted_snps = c(), match_count = c(),
        stringsAsFactors = F)

    temp_result <- as.data.frame(temp_result %>%
        group_by(snp_id, extracted_snps) %>%
        summarise(match_count=sum(match_count), .groups="keep"))
    # Check each SNP set
    for (snp in unique(temp_result$snp_id)){

        subset <- temp_result[temp_result$snp_id == snp,]
        if (length(unique(subset$extracted_snps)) > 1){

            highest_count <- max(subset$match_count)

            if (length(unique(subset[subset$match_count == highest_count, "extracted_snps"])) > 1){
                print("Found 2 match with highest/equal counts, discarding result")
                next
            } else {

                result_to_use <- rbind(result_to_use, subset[subset$match_count == highest_count,])
            }
        } else{

            result_to_use <- rbind(result_to_use, subset)
        }
    }
    if (nrow(result_to_use) > 1){
        snps_found <- as.numeric(strsplit(paste(result_to_use[, "snp_id"], collapse = ", "), ", ")[[1]])
        out_snps <- paste(result_to_use[, "extracted_snps"], collapse = "")
        pattern_ref <- minSNPs:::generate_pattern(seqs, snps_found)
        possible_profiles <- lapply(names(pattern_ref)[which(pattern_ref == out_snps)], function(x) {
            return(strsplit(x, split = "_")[[1]][1])
        })
    } else{
        possible_profiles <- lapply(names(seqs), function(x) {
            return(strsplit(x, split = "_")[[1]][1])
        })
    }
    return(possible_profiles)
}

read_fastq <- function(seq_loc){
    c_pipe <- pipe(
            paste("sed -n '2~4p' ", seq_loc, sep = "")
        )
    data <- readLines(
        c_pipe
    )
    close(c_pipe)
    return(data)
}
