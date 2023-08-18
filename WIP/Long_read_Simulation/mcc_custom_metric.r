# Because if all the profiles that has at least 1 goi in it are considered target profile will inevitably leads to maximising sensitivity, we only consider those profiles that have a majority of goi, i.e.,
# a profile is consider a target profile iff the number of goi >= other non-target sample in the profile
#
calculate_mcc <- function(pattern, goi, MUST_HAVE_TARGET = FALSE) {
    target_seqs <- character()
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
    }
    candidate_target_seqs <- ""
    candidate_count <- 0
    for (pat in unique(pattern)){
        iso_current_pat <- names(which(pattern == pat))
        PAT_in_GOI <- length(which(iso_current_pat %in% goi))
        PAT_not_in_GOI <- length(which(!iso_current_pat %in% goi))
        if (PAT_in_GOI >= PAT_not_in_GOI) { ### Condition for the profile to be added to diagnostic profile
            target_seqs <- c(target_seqs, pat)
        }
        if (MUST_HAVE_TARGET) {
            ## Such that there is always a profile for the targeted samples
            if (PAT_in_GOI > candidate_count) {
                candidate_target_seqs <- pat
                candidate_count <- PAT_in_GOI
            }
        }
    }
    if (MUST_HAVE_TARGET) {
        if (length(target_seqs) < 1) {
                target_seqs <- c(target_seqs, candidate_target_seqs)
            }
    }
    classed_positive <- names(pattern[pattern %in% target_seqs])
    classed_negative <- names(pattern[!pattern %in% target_seqs])
    TP <- as.numeric(length(which(classed_positive %in% goi)))
    TN <- as.numeric(length(which(!classed_negative %in% goi))) ### not in goi and in negative
    FP <- as.numeric(length(which(!classed_positive %in% goi))) ### in goi but not positive
    FN <- as.numeric(length(which(classed_negative %in% goi))) ### negative but in goi
    numerator <- (TP * TN) - (FP * FN)
    denominator <- as.numeric((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- (numerator) / (denominator ** (1/2))
    if (denominator == 0) {
        mcc <- 0
    }
    return(list(result = mcc))
}

parse_group_mcc <- function(pattern, goi, MUST_HAVE_TARGET = FALSE) {
    #data.frame(profile = , is_target =, goi = , non_goi = )
    all_result <- list()
    candidate_target_seqs <- ""
    candidate_count <- 0
    target_seqs <- c()
    u_pattern <- unlist(unique(pattern))
    pat_sample_non_goi <- c()
    pat_sample_goi <- c()
    for (pat in u_pattern){
        iso_current_pat <- names(which(pattern == pat))
        PAT_in_GOI <- length(which(iso_current_pat %in% goi))
        PAT_not_in_GOI <- length(which(!iso_current_pat %in% goi))
        if (PAT_in_GOI >= PAT_not_in_GOI) {
            target_seqs <- c(target_seqs, pat)
        }
        pat_sample_non_goi <- c(pat_sample_non_goi, paste0(iso_current_pat[which(!iso_current_pat %in% goi)], collapse = ", "))
        pat_sample_goi <- c(pat_sample_goi, paste0(iso_current_pat[which(iso_current_pat %in% goi)], collapse = ", "))
        if (MUST_HAVE_TARGET) {
            ## Such that there is always a profile for the targeted samples
            if (PAT_in_GOI > candidate_count) {
                candidate_target_seqs <- pat
                candidate_count <- PAT_in_GOI
            }
        }
    }
    if (MUST_HAVE_TARGET) {
        if (length(target_seqs) < 1) {
                target_seqs <- candidate_target_seqs
            }
    }
    fin_result <- data.frame(profile = u_pattern,
        is_target = u_pattern %in% target_seqs,
        goi = pat_sample_goi,
        non_goi = pat_sample_non_goi)
    sorted_fin_result <- fin_result[with(fin_result, order(-is_target, profile)), ]
    output <- "\nProfile\tIs_Target\tGOI\tNon_GOI\n"
    for (row in seq_len(nrow(sorted_fin_result))){
        output <- paste0(output, paste(sorted_fin_result[row, ], collapse = "\t"), "\n")
    }
    return(output)
}

view_mcc <- function(result, ...) {
    additional_args <- list(...)[[1]]
    print("***")
    print(additional_args)
    if (exists(result$seqc_name) && !is.null(result$goi)) {
        seqc <- get(result$seqc_name)
    } else if (!is.null(additional_args[["seqc"]]) && !is.null(result$goi)) {
        seqc <- additional_args[["seqc"]]
    } else {
        stop("Unable to find the sequences used to generate result OR GOI")
    }

    if (inherits(seqc, "processed_seqs")) {
        seqc <- seqc$seqc
    }
    if (is.null(additional_args[["MUST_HAVE_TARGET"]])) {
        MUST_HAVE_TARGET <- FALSE
    } else {
        MUST_HAVE_TARGET <- additional_args[["MUST_HAVE_TARGET"]]
    }
    results <- result$results
    number_of_result <- length(results)
    goi <- result$goi
    isolates_in_groups <- list()

    for (n in seq_len(number_of_result)) {
        target_seqs <- character()
        isolates_in_groups[[n]] <- list()
        result_levels <- names(results[[n]])
        result_max_depth <- result_levels[length(result_levels)]
        ordered_index <- as.numeric(
            strsplit(result_max_depth, split = ", ")[[1]])
        patterns <- generate_pattern(seqc, ordered_index = c(ordered_index))
        parsed_group <- parse_group_mcc(patterns, goi, MUST_HAVE_TARGET)
        isolates_in_groups[[n]][["groups"]] <- list()
        isolates_in_groups[[n]][["N.B.- groups"]] <- parsed_group
        isolates_in_groups[[n]][["result"]] <- results[[n]]
    }
    return(isolates_in_groups)
}

MinSNPs_metrics <- get_metric_fun()
MinSNPs_metrics[["mcc_single"]] <-
    list(
        "calc" = calculate_mcc,
        "args" = check_percent,
        "view" = view_mcc
    )
MinSNPs_metrics[["mcc_single_2"]] <-
    list(
        "calc" = function(pattern, goi){
            return(calculate_mcc(pattern, goi, MUST_HAVE_TARGET = TRUE))
        },
        "args" = check_percent,
        "view" = function(result, ...){
            args <- list(...)
            args[["MUST_HAVE_TARGET"]] <- TRUE
            return(view_mcc(result, args))
        }
    )

#' \code{calculate_simpson}
#'
#' @description
#' \code{calculate_simpson} is used to calculate Simpson's index.
#' Which is in the range of 0-1, where the greater the value,
#' the more diverse the population.
#' @inheritParams calculate_percent
#' @return Will return the Simpson's index of the list of patterns.
#' @export
calculate_simpson_by_group <- function(pattern, meta, target) {
    setDT(meta)
    colnames(meta)[colnames(meta) == target] <- "target"
    npattern <- data.frame(pattern = unname(unlist(pattern)),
        target = meta[match(names(pattern), meta$isolate), ]$target)
    npattern <- unique(npattern)
    pattern2 <- as.list(npattern$pattern)
    names(pattern2) <- npattern$target
    return(calculate_simpson(pattern2))
}

calculate_simpson_mcc_multi_2 <- function(pattern, meta, target, priority = NULL) {
    mcc <- calculate_mcc_multi(pattern, meta, target, priority)$result
    simpson <- calculate_simpson_by_group(pattern, meta, target)$result
    res <- (simpson + mcc) /2
    return(list(result = res))
}

calculate_simpson_mcc_multi <- function(pattern, meta, target, priority = NULL) {
    mcc <- calculate_mcc_multi(pattern, meta, target, priority)$result
    simpson <- calculate_simpson(pattern)$result
    res <- (simpson + mcc) /2
    return(list(result = res))
}

#' \code{generate_prioritisation}
#' 
#' @description 
#' \code{generate_prioritisation} create a vector of the targets in order of priority.
#' Targets with the less samples are prioritised first, breaking ties by alphabetical order.
#' @param meta A data.table containing the meta data, expect the 
#' @return A data.table of the targets in order of priority
#' @importFrom data.table setDT .N
#' @export
calculate_mcc_multi <- function(pattern, meta, target, priority = NULL) {
    colnames(meta)[which(colnames(meta) == target)] <- "target"
    setDT(meta)
    if (is.null(priority)) {
        priority <- generate_prioritisation(meta[,c("isolate", "target")])
    }
    profile_target <- map_profile_to_target(meta, pattern, priority)
    
    # Translate pattern to target
    isolate <- names(pattern)
    isolate_target <- profile_target[match(unlist(pattern[isolate]), profile_target$profile), "target"][[1]]
    result <- data.frame(isolate = isolate, target = isolate_target, actual = meta[match(isolate, meta$isolate), "target"][[1]])
    setDT(result)
    # Calculate MCC
    predicted <- correct <- NULL
    predicted_table <- meta[, .(actual = .N), by = "target"]
    predicted_table <- merge(predicted_table, result[, .(predicted = .N), "target"], all.x = TRUE )
    predicted_table <- merge(predicted_table, result[, .(correct = sum(target == actual)), "target"], all.x = TRUE )
    predicted_table$predicted[is.na(predicted_table$predicted)] <- 0
    predicted_table$correct[is.na(predicted_table$correct)] <- 0

    numerator <- (sum(predicted_table$actual) * sum(predicted_table$correct)) - crossprod(predicted_table$actual, predicted_table$predicted)[[1]]
    denominator <- ((((sum(predicted_table$actual) ** 2) - crossprod(predicted_table$predicted, predicted_table$predicted) )**(1/2)) * (((sum(predicted_table$actual) ** 2)- crossprod(predicted_table$actual, predicted_table$actual) )**(1/2))) [[1]]
    mcc <- numerator/denominator
    return(list(result = mcc))
}

#' \code{generate_prioritisation}
#' 
#' @description 
#' \code{generate_prioritisation} create a vector of the targets in order of priority.
#' Targets with the less samples are prioritised first, breaking ties by alphabetical order.
#' @param meta A data.table containing the meta data, expect the 
#' @return A data.table of the targets in order of priority
#' @importFrom data.table setDT .N
#' @export
generate_prioritisation <- function(meta) { 
    setDT(meta)
    priority <- target <- NULL
    return(
        data.frame(
            target = as.character(meta[, .(priority = .N), by = target][order(priority, target),target]),
            priority = seq_len(length(unique(meta$target)))
        )
    )
}

#' \code{map_profile_to_target}
#' 
#' @description 
#' \code{map_profile_to_target} creates a mapping of the profile to the target, breaking the ties by the priority.
#' @param patterns A list of the patterns from \code{generate_pattern}
#' @param meta A data.table containing the meta data
#' @param priority A data.table of the targets and priority, either generated by \code{generate_prioritisation} or supplied by user
#' @return A vector of the targets in order of priority
#' @importFrom data.table rbindlist
#' @export
map_profile_to_target <- function(meta, patterns, priority) {
    result_list <- list()
    unique_patterns <- unique(unlist(patterns))
    for (pat in unique_patterns){
        pattern_target_breakdown <- table(meta[meta$isolate %in% names(which(patterns == pat)), "target"][[1]])
        profile_groups <- data.frame(
            target = names(pattern_target_breakdown),
            count = as.numeric(pattern_target_breakdown)
        )
        profile_groups <- merge(profile_groups, priority, by = "target", all.x = TRUE, all.y = FALSE)
        # sort profile_groups
        selected_target <- head(profile_groups[order(-profile_groups$count, profile_groups$priority),], 1)[, "target"][[1]]
        # get the head
        result_list[[pat]] <- data.frame(profile = pat, target = selected_target)
    }
    return(rbindlist(result_list))
}

check_meta_target <- function(list_of_parameters) {
    if (all(c("meta", "target") %in% names(list_of_parameters))){
        if (all(c(list_of_parameters[["target"]], "isolate") %in% colnames(list_of_parameters[["meta"]]))) {
            return(TRUE)
        }
    }
    return(FALSE)
}

parse_group_mcc_multi <- function(result, as_string = TRUE) {
    #data.frame(isolate = , profile = , actual = , group = )
    all_result <- list()
    output <- list()
    for (pat in unique(unlist(result$profile))){
        target <- unique(unlist(result[which(result$profile == pat), "group"]))
        if (length(target)>1) {
            warning("More than one target found for profile ", pat, ". Using the first target.")
            target <- target[1]
        }
        correct_iso <- result[result$profile == pat & result$actual == target, "isolate"]
        incorrect_iso <- result[result$profile == pat & result$actual != target, "isolate"]
        all_result[[pat]] <- data.frame(profile = pat, target = target,
            correct_isolate = paste(correct_iso, collapse = ", "), incorrect_isolate = paste(incorrect_iso, collapse = ", "))
    }
    fin_result <- rbindlist(all_result)
    sorted_fin_result <- fin_result[with(fin_result, order(profile)), ]
    if (as_string){
        output <- "\nProfile\tTarget\tCorrect_Isolate\tIncorrect_Isolate\n"
        for (row in seq_len(nrow(sorted_fin_result))){
            output <- paste0(output, paste(sorted_fin_result[row, ], collapse = "\t"), "\n")
        }
        return(output)
    }
    return(sorted_fin_result)
}

view_mcc_multi <- function(result, ...) {
    additional_args <- list(...)[[1]]
    if (exists(result$seqc_name)) {
        seqc <- get(result$seqc_name)
    } else if (!is.null(additional_args[["seqc"]])) {
        seqc <- additional_args[["seqc"]]
    } else {
        stop("Unable to find the sequences used to generate result")
    }

    if (!is.null(additional_args[["meta"]])) {
        meta <- additional_args[["meta"]]
    } else {
        stop("Unable to find the meta")
    }

    if (!is.null(additional_args[["target"]])) {
        target <- additional_args[["target"]]
    } else {
        stop("Unable to find the target")
    }

    setDT(meta)
    colnames(meta)[colnames(meta) == target] <- "target"
    if (!is.null(additional_args[["priority"]])) {
        priority <- additional_args[["priority"]]
    } else {
        warning("Unable to find the priority list, using default method to generate")
        priority <- generate_prioritisation(meta[,c("isolate", "target")])
    }

    if (inherits(seqc, "processed_seqs")) {
        seqc <- seqc$seqc
    }

    results <- result$results
    number_of_result <- length(results)
    isolates_in_groups <- list()

    for (n in seq_len(number_of_result)) {
        target_seqs <- character()
        isolates_in_groups[[n]] <- list()
        result_levels <- names(results[[n]])
        result_max_depth <- result_levels[length(result_levels)]
        ordered_index <- as.numeric(
            strsplit(result_max_depth, split = ", ")[[1]])
        patterns <- generate_pattern(seqc, ordered_index = c(ordered_index))
        
        profile_target <- map_profile_to_target(meta, patterns, priority)
        isolate <- names(patterns)
        isolate_target <- profile_target[match(unlist(patterns[isolate]), profile_target$profile), "target"][[1]]
        pred_result <- data.frame(isolate = isolate, profile = unlist(patterns[isolate]), actual = meta[match(isolate, meta$isolate), "target"][[1]], group = isolate_target)
        parsed_group <- parse_group_mcc_multi(pred_result)
        
        isolates_in_groups[[n]][["groups"]] <- list()
        isolates_in_groups[[n]][["N.B.- groups"]] <- parsed_group
        isolates_in_groups[[n]][["result"]] <- results[[n]]
    }
    return(isolates_in_groups)
}

MinSNPs_metrics <- get_metric_fun()
MinSNPs_metrics[["mcc_multi"]] <-
    list(
        "calc" = calculate_mcc_multi,
        "args" = check_meta_target,
        "view" = view_mcc_multi
    )

MinSNPs_metrics[["simpson_mcc_multi"]] <-
    list(
        "calc" = calculate_simpson_mcc_multi,
        "args" = check_meta_target,
        "view" = view_mcc_multi
    )

MinSNPs_metrics[["simpson_mcc_multi_2"]] <-
    list(
        "calc" = calculate_simpson_mcc_multi_2,
        "args" = check_meta_target,
        "view" = view_mcc_multi
    )

MinSNPs_metrics[["simpson_by_group"]] <-
    list(
        "calc" = calculate_simpson_by_group,
        "args" = check_meta_target,
        "view" = view_mcc_multi
    )


### Identify the target from mcc_multi that failed to be identified by the profile
identify_mixed_targets <- function(snps, meta, target, priority = NULL, seqc){
    setDT(meta)
    colnames(meta)[colnames(meta) == target] <- "target"
    if (is.null(priority)) {
        priority <- generate_prioritisation(meta[,c("isolate", "target")])
    }
    patterns <- generate_pattern(seqc, ordered_index = c(snps))
    profile_target <- map_profile_to_target(meta, patterns, priority)
    isolate <- names(patterns)
    isolate_target <- profile_target[match(unlist(patterns[isolate]), profile_target$profile), "target"][[1]]
    pred_result <- data.frame(
        isolate = isolate,
        profile = unlist(patterns[isolate]),
        actual = meta[match(isolate, meta$isolate), "target"][[1]],
        group = isolate_target)
    profile_isolates <- parse_group_mcc_multi(pred_result, as_string = FALSE)

    incorrect <- profile_isolates[profile_isolates$incorrect_isolate != "",]
    if (nrow(incorrect) == 0){
        return(
            list(
                missed_target_count = data.frame(target = c(), n_missed =c(), prop = c()),
                mixed_target = data.frame(mixed_group = c(), n_isolate = c(), n_target = c())
            )
        )
    }
    results <- lapply(seq_len(nrow(incorrect)), function(i){
        row <- incorrect[i,]
        mixed_iso <- strsplit(row$incorrect_isolate, split = ", ")[[1]]
        mixed_target <- table(meta[meta$isolate %in% mixed_iso,]$target)
        mt_names <- names(mixed_target)
        mt_count <- as.numeric(mixed_target)
        return(
            data.frame(mixed_group =
                paste0(sort(c(row$target, mt_names)), collapse = " -> "),
            target = row$target, 
            mixed = mt_names,
            mixed_count = mt_count))
    })
    res <- rbindlist(results)
    mixed_count = n_isolate = mixed_group = NULL
    res <- res[, .(n_isolate = sum(mixed_count)), mixed_group]
    res$n_target <- sapply(strsplit(res$mixed_group, split = " -> "), length)
    missed_target_count <- as.data.frame(table(meta[match(unlist(strsplit(incorrect$incorrect_isolate, split = ", ")), meta$isolate), ]$target))
    colnames(missed_target_count) <- c("target", "n_missed")
    missed_target_count$prop <- missed_target_count$n_missed / sapply(missed_target_count$target, function(x){return(nrow(meta[meta$target == x]))})

    fin_result <- list(
        missed_target_count = missed_target_count[order(-missed_target_count$prop, -missed_target_count$n_missed),],
        mixed_target = res[order(-res$n_target, -res$n_isolate),]
    )
    return(fin_result)
}

### data.frame
# SNPs_set | score
### List
# [[1]] SNP1, SNP2, SNP3
# [[N]] SNP1, SNP2, SNP3
# attr(SNPs, "score") <- c(score1, scoreN)
#
get_snps_set <- function(results, as = "data.frame") {
    if (!as %in% c("data.frame", "list"))
        stop("as must be either data.frame or list")

    snp_sets <- sapply(results$results, function(x){
        return(
            tail(names(x), n = 1)
        )
    })
    score <- sapply(results$results, function(x){
        return(
            tail(x, n = 1)[[1]]
        )
    })
    if (as == "data.frame"){
        result <- data.frame(SNPs_set = snp_sets, score = score)
    } else {
        result <- unname(lapply(strsplit(snp_sets, split = ", "), as.numeric))
        attr(result, "score") <- unname(score)
    }
    return(result)
}

read_vcf <- function(file, header_line = NULL) {
    read_data <- readLines(file)
    if (is.null(header_line)){
        header_line <- max(which(startsWith(read_data, "#")))
    }
    header <- gsub("#", "", read_data[header_line])
    data <- paste(c(header, read_data[(header_line + 1):length(read_data)]), collapse = "\n")
    result <- read.table(text = data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(result)
}

dual_bases <- list("AC"="M",
                   "AG"="R",
                   "AT"="W",
                   "CG"="S",
                   "CT"="Y",
                   "GT"="K")

extract_gt <- function(data, format, ref, alt, ads, min_ad = 1, min_alt = 1, het_ratio = 0.67, major = TRUE){
    data2 <- strsplit(unlist(data), split = ":")
    format_ind <- match(c("GT", "AD"), strsplit(format, split = ":")[[1]])
    gts <- lapply(data2, function(x){
        strsplit(x[format_ind[1]], split = "/")[[1]]
    })
    ads <- sapply(data2, function(x){
        strsplit(x[format_ind[2]], split = ",")[[1]]
    })
    r_depth <- as.numeric(ads[1,])
    a_depth <- as.numeric(ads[2,])
    total_depth <- r_depth + a_depth
    major <- 
    minor <- 
    if (!major){

    }else{
        
    }
}

vcf_data_to_fasta <- function(vcf_tab, min_ad = 1, min_alt = 1, het_ratio = 0.67, major = TRUE, bp = MulticoreParam()) {
    infocols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    stopifnot(all(infocols %in% colnames(vcf_tab)))
    isolates <- colnames(vcf_tab)[colnames(vcf_tab) %in% infocols == FALSE]
    FORMAT <- vcf_tab$FORMAT
    REF <- vcf_tab$REF
    ALT <- vcf_tab$ALT
    n_iso <- length(isolates)
    n_snps <- nrow(vcf_tab)
    POS <- paste0(vcf_tab$CHROM, ":", vcf_tab$POS)
    result <- bplapply(seq_len(nrow(vcf_tab)), function(r){
        data <- vcf_tab[r,isolates]
        return(extract_gt(data, FORMAT[r], REF[r], ALT[r], min_ad, min_alt, het_ratio, major))
    }, BPPARAM = bp)
    result <- do.call(cbind, result)
    result <- lapply(seq_len(nrow(result)), function(r){
        paste0(result[r,], collapse = "")
    })
    names(result) <- isolates
    return(
        list(seqc = result,
        pos = POS
        )
    )
}