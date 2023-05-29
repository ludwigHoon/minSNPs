# Because if all the profiles that has at least 1 goi in it are considered target profile will inevitably leads to maximising sensitivity, we only consider those profiles that have a majority of goi, i.e.,
# a profile is consider a target profile iff the number of goi >= other non-target sample in the profile
calculate_mcc <- function(pattern, goi) {
    target_seqs <- character()
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
    }
    for (pat in unique(pattern)){
        iso_current_pat <- names(which(pattern == pat))
        PAT_in_GOI <- length(which(iso_current_pat %in% goi))
        PAT_not_in_GOI <- length(which(!iso_current_pat %in% goi))
        if (PAT_in_GOI >= PAT_not_in_GOI) { ### Condition for the profile to be added to diagnostic profile
            target_seqs <- c(target_seqs, pat)
        }
    }
    classed_positive <- names(pattern[pattern %in% target_seqs])
    classed_negative <- names(pattern[!pattern %in% target_seqs])
    TP <- length(which(classed_positive %in% goi))
    TN <- length(which(!classed_negative %in% goi)) ### not in goi and in negative
    FP <- length(which(!classed_positive %in% goi)) ### in goi but not positive
    FN <- length(which(classed_negative %in% goi)) ### negative but in goi
    numerator <- (TP * TN) - (FP * FN)
    denominator <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
    mcc <- (numerator) / (denominator ** (1/2))
    return(list(result = mcc))
}

parse_group_mcc <- function(pattern, goi) {
    #data.frame(profile = , is_target =, goi = , non_goi = )
    all_result <- list()
    for (pat in unlist(unique(pattern))){
        is_target <- FALSE
        iso_current_pat <- names(which(pattern == pat))
        PAT_in_GOI <- length(which(iso_current_pat %in% goi))
        PAT_not_in_GOI <- length(which(!iso_current_pat %in% goi))
        if (PAT_in_GOI >= PAT_not_in_GOI) {
            is_target <- TRUE
        }
        all_result[[pat]] <- data.frame(profile = pat, is_target = is_target, goi = paste(PAT_in_GOI, collapse = ", "), non_goi = paste(PAT_not_in_GOI, collapse = ", "))
    }
    fin_result <- rbindlist(all_result)
    sorted_fin_result <- fin_result[with(fin_result, order(-is_target, profile)), ]
    output <- "Profile\tIs_Target\tGOI\tNon_GOI\n"
    for (row in seq_len(nrow(sorted_fin_result))){
        output <- paste0(output, paste(sorted_fin_result[row, ], collapse = "\t"), "\n")
    }
    return(output)
}

view_mcc <- function(result, ...) {
    additional_args <- list(...)[[1]]
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
        parsed_group <- parse_group_mcc(patterns, goi)
        isolates_in_groups[[n]][["groups"]] <- parsed_group
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


calculate_simpson_mcc_multi <- function(pattern, meta, target, priority = NULL) {
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
    predicted_table <- merge(predicted_table, result[, .(predicted = .N), by = "target"], all.x = TRUE )
    predicted_table <- merge(predicted_table, result[, .(correct = sum(target == actual)), "target"], all.x = TRUE )
    predicted_table$predicted[is.na(predicted_table$predicted)] <- 0
    predicted_table$correct[is.na(predicted_table$correct)] <- 0

    numerator <- (sum(predicted_table$actual) * sum(predicted_table$correct)) - crossprod(predicted_table$actual, predicted_table$predicted)[[1]]
    denominator <- ((((sum(predicted_table$actual) ** 2) - crossprod(predicted_table$predicted, predicted_table$predicted) )**(1/2)) * (((sum(predicted_table$actual) ** 2)- crossprod(predicted_table$actual, predicted_table$actual) )**(1/2))) [[1]]
    mcc <- numerator/denominator

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
    print(length(list_of_parameters))
    print(names(list_of_parameters))

    if (all(c("meta", "target") %in% names(list_of_parameters))){
        if (all(c(list_of_parameters[["target"]], "isolate") %in% colnames(list_of_parameters[["meta"]]))) {
            return(TRUE)
        }
    }
    print(all(c("meta", "target") %in% names(list_of_parameters)))
    print(all(c(list_of_parameters[["target"]], "isolate") %in% colnames(list_of_parameters[["meta"]])))
    return(FALSE)
}

view_mcc_multi <- function(result, ...) {
    additional_args <- list(...)[[1]]
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

    if (is.null(additional_args[["priority"]])) {
        priority <- generate_prioritisation(seqc$meta[,c("sample", "target")])
    } else {
        priority <- additional_args[["priority"]]
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
        parsed_group <- parse_group_mcc(patterns, goi)
        isolates_in_groups[[n]][["groups"]] <- parsed_group
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


### RESULT view:
# result - 
#    [[snps]] score
#    [[snps 2]] score
# groups - 
#    [[group ID]] isolates
### paste(group ID, isolate, sep="\t")
# other_details: N.B.-
#    N.B.- ......: paste(value, collapse=", ")

#   ./gsk-test/njtree/sim+gsk-229317.recode.vcf.gz

#source /data/Public/htrimarsanto/opt/etc/bashrc_py3.11
/data/Public/htrimarsanto/env/pys/wgs (env)
spcli

/data/malaria/work/anto/pv-geobarcode/vivaxgen-geo/training-sets/Pv40-943