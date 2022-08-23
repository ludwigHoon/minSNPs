### COMMON FUNCTIONS

### HELPERS
get_data <- function(file_path, samples = 0){
    c_pipe <- pipe(
        paste("sed -n '2~4p' ", file_path, sep = "")
    )
    data <- readLines(
        c_pipe
    )
    close(c_pipe)
    if (samples > 0){
        return(data[samples])
    }
    return(data)
}

### ACTUAL FUNCTIONS
#### STEP 1.
matched_snp <- function(ref, searches){
    s_strings <- unique(searches$search_string)
    c_count <- lapply(s_strings, function(string) {
            match <- unlist(gregexpr(string, ref))
            return(match[which(match > 0)])
        })
    s_m <- lapply(c_count, length)
    matched_string <- s_strings[which(s_m > 0)]
    match_count <- unlist(s_m[which(s_m > 0)])
    match_indexes <- unlist(lapply(c_count[which(s_m > 0)], paste, collapse = ", "))
    snp_id <- searches[match(matched_string, searches$search_string), "fasta_position"]
    #target_pos <- searches[match(matched_string, searches$search_string), "target_position"]
    ext_snps <- searches[match(matched_string, searches$search_string), "snp_string"]
    stopifnot(length(matched_string) == length(match_count))
    stopifnot(length(matched_string) == length(match_indexes))
    result <- data.frame(snp_id = snp_id,
        extracted = ext_snps, count = match_count,
        stringsAsFactors = F)
    return(result)
}

#step_1_res <- data.frame(isolate = c(), read_id = c(), snp_id = c(), 
#    count = c(), extracted = c(),
#    stringsAsFactor = F)

#### STEP 2.
combine_snps <- function(result_list, matrix){
    result_to_use <- data.frame(isolate = c(), snp_id = c(),
        extracted = c(), matched_cc = c(),
        stringsAsFactors = F)
    
    temp_data <- as.data.frame(result_list %>%
        group_by(snp_id, extracted) %>%
        summarise(count=sum(count), .groups="keep"))

    all_snps <- unique(temp_data$snp_id)
    for (snp in all_snps){
        subset <- temp_data[temp_data$snp_id == snp,]
        if (length(unique(subset$extracted)) > 1){
            highest_count <- max(subset$count)
            if (length(unique(subset[subset$count == highest_count, "extracted"])) > 1){
                print("Found 2 match with highest/equal counts, discarding result")
                next
            } else {
                base_result <- subset[subset$count == highest_count, ]
                out_snps <- subset[subset$count == highest_count, "extracted"]
            }
        } else{
            out_snps <- subset$extracted
            base_result <- subset
        }
        snps_found <- as.numeric(strsplit(snp, split = ", ")[[1]])
        pattern_ref <- minSNPs:::generate_pattern(matrix, snps_found)
        base_result$matched_cc <- paste(unlist(names(pattern_ref)[which(pattern_ref == out_snps)]), collapse = ", ")
        
        base_result$isolate <- unique(result_list$isolate)[1]
        base_result <- base_result[,c("isolate", "snp_id", "extracted", "count", "matched_cc")]
        result_to_use <- rbind(result_to_use, base_result)
    }
    return(result_to_use)
}

#### STEP 3.
infer_most_likely <- function(processed_result){
    if (ncol(processed_result) == 0){
        return(NULL)
    }
    matched_cc <- strsplit(processed_result$matched_cc, split = ", ")
    snp_diff <- 0
    cc_m_count <- table(unlist(matched_cc))
    m_count <- max(cc_m_count)
    snp_count <- as.vector(cc_m_count)[order(as.vector(cc_m_count), decreasing=T)]
    if (length(snp_count) > 1){
        snp_diff <- snp_count[1] - snp_count[2]
    }
    cc_group_count <- unique(cc_m_count)
    cc_group_diff <- diff(cc_group_count[order(cc_group_count)])
    res <- list(most_likely = names(which(cc_m_count == m_count)),
        n_candidates = length(names(cc_m_count)),
        min_diff = snp_diff,
        snp_count = m_count,
        group_count = paste(cc_group_count[order(cc_group_count)], collapse = "_"),
        group_diff = paste(cc_group_diff, collapse ="_"))
    return(res)
}