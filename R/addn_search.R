#' \code{calculate_mcc}
#'
#' @description
#' \code{calculate_mcc} is used to calculate the MCC score given the SNP profile.
#' @param pattern the SNP profile for each samples
#' @param goi the samples belonging to the group of interest
#' @param MUST_HAVE_TARGET whether to force the profile to have at least 1 target profile
#' (the profile containing the most goi)
#' @return the MCC score
#' @export
calculate_mcc <- function(pattern, goi, MUST_HAVE_TARGET = TRUE) {
    # Because if all the profiles that has at least 1 goi in it are considered target profile
    # it will inevitably leads to maximising sensitivity,
    # we only consider those profiles that have the number of goi >= other non-target sample 
    # to be a target profile and additional constraining that there is at least a target profile
    # (the profile containing the most goi)
    target_seqs <- character()
    for (isolate in goi) {
        if (is.null(pattern[[isolate]])) {
            stop(paste(isolate, " in group of interest, ",
                "but not found in list of isolates", sep = ""))
        }
    }
    candidate_target_seqs <- ""
    candidate_count <- 0
    candidate_fp <- 0
    for (pat in unique(pattern)){
        iso_current_pat <- names(which(pattern == pat))
        PAT_in_GOI <- length(which(iso_current_pat %in% goi))
        PAT_not_in_GOI <- length(which(!iso_current_pat %in% goi))
        if (PAT_in_GOI >= PAT_not_in_GOI) { ### Condition for the profile to be added to diagnostic profile
            target_seqs <- c(target_seqs, pat)
        }
        if (MUST_HAVE_TARGET) {
            ## Such that there is always a profile for the targeted samples
            if (candidate_target_seqs == "") {
                candidate_target_seqs <- pat
                candidate_count <- PAT_in_GOI
                candidate_fp <- PAT_not_in_GOI
            } else {
                if (PAT_in_GOI >= candidate_count && PAT_not_in_GOI < candidate_fp) {
                    candidate_fp <- PAT_not_in_GOI
                    candidate_target_seqs <- pat
                    candidate_count <- PAT_in_GOI
                }
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
    sme <- 0.000001
    denominator <- as.numeric((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) + sme
    mcc <- (numerator) / (denominator ** (1/2))

    return(list(result = mcc))
}

#' \code{parse_group_mcc}
#'
#' @description
#' \code{parse_group_mcc} is used to group the sample according to SNPs profile and present in a table format
#' @param pattern the SNP profile for each samples
#' @param goi the samples belonging to the group of interest
#' @param MUST_HAVE_TARGET whether to force the profile to have at least 1 target profile
#' (the profile containing the most goi)
#' @return the parsed group views
#' @export
parse_group_mcc <- function(pattern, goi, MUST_HAVE_TARGET = TRUE) {
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
            if (candidate_target_seqs == "") {
                candidate_target_seqs <- pat
                candidate_count <- PAT_in_GOI
                candidate_fp <- PAT_not_in_GOI
            } else {
                if (PAT_in_GOI >= candidate_count && PAT_not_in_GOI < candidate_fp) {
                    candidate_fp <- PAT_not_in_GOI
                    candidate_target_seqs <- pat
                    candidate_count <- PAT_in_GOI
                }
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

#' \code{view_mcc}
#'
#' @description
#' \code{view_mcc} is used to present the result of selected SNPs set
#' based on the MCC score
#' @param result is the result from \code{find_optimised_snps}
#' @param ... other optional parameters
#' @return formatted result list to be saved or presented.
#' @export
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
    if (is.null(additional_args[["MUST_HAVE_TARGET"]])) {
        MUST_HAVE_TARGET <- TRUE
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

#' \code{check_meta_target}
#'
#' @description
#' \code{check_meta_target} is used to check if parameters needed by
#' \code{calculate_mcc_multi} and \code{simpson_by_group} are all present.
#' @param list_of_parameters is a list of parameter passed
#' to functions that will perform the calculation
#' @return TRUE if the parameters exists, else FALSE
#' @export
check_meta_target <- function(list_of_parameters) {
    if (all(c("meta", "target") %in% names(list_of_parameters))){
        if (all(c(list_of_parameters[["target"]], "isolate") %in% colnames(list_of_parameters[["meta"]]))) {
            return(TRUE)
        }
    }
    return(FALSE)
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
            target = as.character(meta[, list(priority = .N), by = target][order(priority, target),target]),
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
#' @param sensitive_to_1 whether to be completely sensitive to the first target (percent default), set to TRUE for percent
#' @return A vector of the targets in order of priority
#' @importFrom data.table rbindlist setDT
#' @importFrom utils head
#' @export
map_profile_to_target <- function(meta, patterns, priority, sensitive_to_1 = FALSE) {
    result_list <- list()
    setDT(meta)
    unique_patterns <- unique(unlist(patterns))
    if (sensitive_to_1){
        stopifnot(nrow(priority) == 2)
        GOI <- priority[1, ]$target
        NON_GOI <- priority[2, ]$target
    }
    for (pat in unique_patterns){
        if (sensitive_to_1){
            if (any(meta[meta$isolate %in% names(which(patterns == pat))]$target %in% GOI)){
                result_list[[pat]] <- data.frame(profile = pat, target = GOI)
            } else{
                result_list[[pat]] <- data.frame(profile = pat, target = NON_GOI)
            }
        } else {
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
    }
    return(rbindlist(result_list))
}

#' \code{calculate_mcc_multi}
#' 
#' @description 
#' \code{calculate_mcc_multi} Calculate the multi-class MCC score for the SNPs.
#' It assigns each SNP profile to a class, based on the majority of the samples having the profile,
#' targets with the less samples are prioritised first, breaking ties by alphabetical order.
#' @param pattern the SNP profile for each samples
#' @param meta A data.table containing the meta data
#' @param target the column name of the target in the meta data, default to target
#' @param priority A data.table of the targets and priority,
#' either supplied by user, or by default generated by \code{generate_prioritisation}.
#' @return multiclass-MCC score
#' @importFrom data.table setDT .N
#' @export
calculate_mcc_multi <- function(pattern, meta, target = "target", priority = NULL) {
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
    predicted <- correct <- actual <- NULL
    predicted_table <- meta[, list(actual = .N), by = "target"]
    predicted_table <- merge(predicted_table, result[, list(predicted = .N), "target"], all.x = TRUE )
    predicted_table <- merge(predicted_table, result[, list(correct = length(which(target == actual))), "target"], all.x = TRUE )
    predicted_table$predicted[is.na(predicted_table$predicted)] <- 0
    predicted_table$correct[is.na(predicted_table$correct)] <- 0
    sme <- 0.000001
    numerator <- (sum(predicted_table$actual) * sum(predicted_table$correct)) - crossprod(predicted_table$actual, predicted_table$predicted)[[1]]
    denominator <- ((((sum(predicted_table$actual) ** 2) - crossprod(predicted_table$predicted, predicted_table$predicted) )**(1/2)) * (((sum(predicted_table$actual) ** 2)- crossprod(predicted_table$actual, predicted_table$actual) )**(1/2)) +sme) [[1]]
    mcc <- numerator/denominator
    return(list(result = mcc))
}

#' \code{parse_group_mcc_multi}
#'
#' @description
#' \code{parse_group_mcc_multi} is used to put samples according to SNP profile, and
#' put them into a table format.
#' @param result result from \code{find_optimised_snps}
#' @param as_string whether to return the result as string or data.frame
#' @return Will return the grouped samples.
#' @export
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


#' \code{view_mcc_multi}
#'
#' @description
#' \code{view_mcc_multi} is used to present the result of selected SNPs set
#' based on the multi-MCC score
#' @param result is the result from \code{find_optimised_snps}
#' @param ... other optional parameters
#' @return formatted result list to be saved or presented.
#' @export
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
        reordered <- match(isolate, meta$isolate)
        pred_result <- data.frame(isolate = isolate, profile = unlist(patterns[isolate]), actual = meta[reordered, "target"][[1]], group = isolate_target)
        parsed_group <- parse_group_mcc_multi(pred_result)
        
        isolates_in_groups[[n]][["groups"]] <- list()
        isolates_in_groups[[n]][["N.B.- groups"]] <- parsed_group
        isolates_in_groups[[n]][["result"]] <- results[[n]]
    }
    return(isolates_in_groups)
}

#' \code{calculate_simpson_by_group}
#'
#' @description
#' \code{calculate_simpson_by_group} is used to calculate Simpson's index.
#' Which is in the range of 0-1, where the greater the value,
#' the more diverse the population.
#' @inheritParams calculate_simpson
#' @param meta the metadata
#' @param target the target column name
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



#### Utilities ####

#' \code{profile_to_group_result}
#'
#' @description
#' \code{profile_to_group_result} given profile target, return the result
#' @param patterns the SNP profile for each samples
#' @param profile_target the profile target - should be from samples previously seen, generate with \code{map_profile_to_target}
#' @return Will return the result, given the SNP profile.
#' @export
profile_to_group_result <- function(patterns, profile_target){
    t_patterns <- data.frame(isolate = names(patterns), profile = unlist(patterns))
    result <- merge(t_patterns, profile_target, by = "profile", all.x = TRUE)
    colnames(result)[which(colnames(result) == "target")] <- "predicted_target"
    return(result)
}

#' \code{mcc_calculation}
#'
#' @description
#' \code{mcc_calculation} calculate the MCC score given the truth and predicted target.
#' @param result_with_truth the dataframe containing the truth and predicted target
#' @param is_multi Whether to use MCC-multi or MCC
#' @param return_all_intermediate whether to return all intermediate values, only possible for binary class
#' @importFrom data.table setDT .N
#' @return Will return the mcc score
#' @export
mcc_calculation <- function(result_with_truth, is_multi = TRUE, return_all_intermediate = FALSE) {
    predicted <- correct <- predicted_target <- target <- actual <- NULL
    untypeable <- 0
    setDT(result_with_truth)

    if (is_multi){
        predicted_table <- result_with_truth[, list(actual = .N), by = "target"]
        predicted_table <- merge(predicted_table, result_with_truth[, list(predicted = .N), "predicted_target"], 
            by.x = "target", by.y = "predicted_target", all = TRUE )
        predicted_table <- merge(predicted_table, result_with_truth[, list(correct = length(which(predicted_target == target))), "target"], all.x = TRUE )
        predicted_table$predicted[is.na(predicted_table$predicted)] <- 0
        predicted_table$correct[is.na(predicted_table$correct)] <- 0
        predicted_table$actual[is.na(predicted_table$actual)] <- 0

        if (nrow(predicted_table[is.na(target)]) > 0){
            untypeable <- sum(predicted_table[is.na(target)]$predicted)
            predicted_table <- predicted_table[!is.na(target)]
        }
        sme <- 0.000001
        numerator <- (sum(predicted_table$actual) * sum(predicted_table$correct)) - crossprod(predicted_table$actual, predicted_table$predicted)[[1]]
        denominator <- ((((sum(predicted_table$actual) ** 2) - crossprod(predicted_table$predicted, predicted_table$predicted) )**(1/2)) * (((sum(predicted_table$actual) ** 2)- crossprod(predicted_table$actual, predicted_table$actual) )**(1/2)) + sme) [[1]]
        mcc <- numerator/denominator

        return(data.frame(mcc = mcc, untypeable = untypeable))
    } else {
        # Calculate the MCC for each class
        mcc_each_target <- c()
        truth_target <- unique(result_with_truth$target)
        if (return_all_intermediate){
            stopifnot(length(truth_target) == 2)
            stopifnot("GOI" %in% truth_target)
            TPs <- TNs <- FPs <- FNs <- c()
            goi_index <- which(truth_target == "GOI") 
        }
        for (i in truth_target){
            TP <- as.numeric(nrow(result_with_truth[target == i & predicted_target == i]))
            TN <- as.numeric(nrow(result_with_truth[target != i & predicted_target != i]))
            FP <- as.numeric(nrow(result_with_truth[target != i & predicted_target == i]))
            FN <- as.numeric(nrow(result_with_truth[target == i & predicted_target != i]))
            numerator <- (TP * TN) - (FP * FN)
            sme <- 0.000001
            denominator <- as.numeric((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) + sme
            mcc <- (numerator) / (denominator ** (1/2))
            mcc_each_target <- c(mcc_each_target, mcc)
            if (return_all_intermediate){
                TPs <- c(TPs, TP)
                TNs <- c(TNs, TN)
                FPs <- c(FPs, FP)
                FNs <- c(FNs, FN)
            }
        }
        result_df <- setNames(data.frame(t(mcc_each_target)), truth_target)
        if (return_all_intermediate){
            result_df <- cbind(result_df, data.frame(mcc = mcc_each_target[goi_index], TP = TPs[goi_index], TN = TNs[goi_index], FP = FPs[goi_index], FN = FNs[goi_index]))
            result_df$sensitivity <- 
                result_df$TP / (result_df$TP + result_df$FN)
            result_df$specificity <-
                result_df$TN / (result_df$TN + result_df$FP)
        }
        result_df$untypeable <- nrow(result_with_truth[is.na(predicted_target)])
        return(result_df)
    }
}


#' \code{summarise_result}
#'
#' @description
#' \code{summarise_result} calculate the MCC score given the SNP sets, training, validation and metadata(s).
#' @param snp_sets the dataframe containing the truth and predicted target
#' @param training_seqs the training sequences
#' @param validation_seqs the validation sequences
#' @param training_metadata the training metadata
#' @param validation_metadata the validation metadata
#' @param priority the priority of the target, generated by \code{generate_prioritisation}
#' @param is_multi Whether to use MCC-multi or MCC
#' @param return_all_intermediate whether to return all intermediate values, only possible for binary class
#' @param is_percent whether to return result by considering all the profiles having a GOI as target profile
#' @importFrom data.table rbindlist
#' @return Will return the summarised result
#' @export
summarise_result <- function(snp_sets, training_seqs, validation_seqs, training_metadata, validation_metadata, priority, is_multi = TRUE, return_all_intermediate = FALSE, is_percent = FALSE){
    l_result <- list()
    for (sset in snp_sets) {
        # all the needed seqs
        training_pat <- generate_pattern(training_seqs, sset)
        validation_pat <- generate_pattern(validation_seqs, sset)
        all_pat <- c(training_pat, validation_pat)

        # use the training only to generate the profile map 
        profile_target <- map_profile_to_target(training_metadata, training_pat, priority, is_percent)
        
        # mcc score for training only
        temp_training <- profile_to_group_result(training_pat, profile_target)
        temp_training <- merge(temp_training, training_metadata, by = "isolate", all.x = TRUE)
        ## MCC score calculation here
        tr_score <- mcc_calculation(temp_training, is_multi, return_all_intermediate)

        # mcc score for validation only
        temp_validation <- profile_to_group_result(validation_pat, profile_target)
        temp_validation <- merge(temp_validation, validation_metadata, by = "isolate", all.x = TRUE)
        va_score <- mcc_calculation(temp_validation, is_multi, return_all_intermediate)

        # mcc score for training and validation
        temp_all <- profile_to_group_result(all_pat, profile_target)
        metadata <- rbindlist(list(training_metadata, validation_metadata), fill = TRUE)
        temp_all <- merge(temp_all, metadata, by = "isolate", all.x = TRUE)
        co_score <- mcc_calculation(temp_all, is_multi, return_all_intermediate)

        if (is_multi){
            l_result[[length(l_result) + 1]] <- data.frame(snp_sets = paste(sset, collapse = ", "), tr_score = tr_score$mcc, va_score = va_score$mcc, co_score = co_score$mcc, 
                   tr_untypeable = tr_score$untypeable, va_untypeable = va_score$untypeable, co_untypeable = co_score$untypeable)
        } else {
            temp_result_df <- rbindlist(list(tr_score, va_score, co_score), fill = TRUE)
            temp_result_df$snp_sets <- paste(sset, collapse = ", ")
            temp_result_df$partition <- c("training", "validation", "all")
            l_result[[length(l_result) + 1]] <- temp_result_df
        }
    }
    return(rbindlist(l_result))
}