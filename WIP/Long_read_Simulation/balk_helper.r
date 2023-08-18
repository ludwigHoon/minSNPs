std_bases <- c("A", "C", "G", "T")
dual_bases <- list("AC"="M",
                   "AG"="R",
                   "AT"="W",
                   "CG"="S",
                   "CT"="Y",
                   "GT"="K")

transform_snp <- function(pat, binomial_n, levels = c(), get=c("levels", "transformed")){
    pat <- toupper(pat)
    unique_bases <- unique(pat)
    pat_ambig_bases <- unique_bases[unique_bases %in% dual_bases]
    pat_std_bases <- sort(unique_bases[unique_bases %in% std_bases])

    ### Determine the transformation levels
    if (length(levels) == 0){
        if (length(unique_bases) > (binomial_n + 1))
            stop("Unable to transform: unique_bases > binomial_n")
        if (length(unique_bases) <= 1)
            stop("Unable to transform: unique_bases <= 1")
        if (length(pat_ambig_bases) > 1)
            stop("Unable to transform: multiple ambiguity codes found")
        if (length(pat_std_bases) == 2 && length(pat_ambig_bases) == 1)
            stopifnot(dual_bases[paste0(pat_std_bases, collapse = "")] == pat_ambig_bases)
        if (length(pat_std_bases) == 1){
            inferred_std_bases <- strsplit(names(dual_bases[which(dual_bases == pat_ambig_bases)]), split = "")[[1]]
            stopifnot(pat_std_bases %in% inferred_std_bases)
            pat_std_bases <- inferred_std_bases
        }
        if (binomial_n == 2 && length(pat_ambig_bases) == 0){
            pat_ambig_bases <- dual_bases[paste0(pat_std_bases, collapse = "")][[1]]
        }
        all_bases <- c(pat_std_bases[1], pat_ambig_bases, pat_std_bases[2])
        attr(all_bases, "order") <- seq_along(all_bases) - 1
        levels <- all_bases
    }
    if (all(get == "levels"))
        return(levels)
    ### Transformation steps
    temp_levels <- levels
    if (binomial_n == 1 && length(levels) > 2)
        temp_levels <- temp_levels[c(1,3)]
    stopifnot(length(temp_levels) == (binomial_n + 1))
    if (binomial_n == 2){
        pat <- gsub("N", levels[2], pat)
    }
    transformed_pat <- match(pat, temp_levels) - 1
    if (all(get == "transformed"))
        return(transformed_pat)
    return(list(transformed_pat = transformed_pat,
                levels = levels))
}

train_balk <- function(seqc, snps_pos, meta, binomial_n = 1, laplace = 1, snp_id = NULL, prior = NULL,
                    fit_prior = FALSE){
    if (!all(names(seqc) %in% meta$isolate)){
        stop("Not all sequences in seqc are found in the metadata")
    }
    if (!is.null(snp_id)){
        if((length(snp_id) != length(snps_pos)))
            stop("snp_id & snps_pos lengths differ")
    } else {
        snp_id <- as.character(snps_pos)
    }
    rows <- names(seqc)
    columns <- snp_id 
    levels_transformed <- lapply(snps_pos, function(x){
        pat <- minSNPs::generate_pattern(seqc, x)
        return(transform_snp(pat, binomial_n))
    })
    X <- sapply(levels_transformed, function(x){
        return(rbind(x[["transformed_pat"]]))
    })
    colnames(X) <- columns
    rownames(X) <- rows
    Y <- meta[match(rows, meta$isolate),]$target
    levels <- sapply(levels_transformed, function(x){
        return(rbind(x[["levels"]]))
    })
    colnames(levels) <- columns
    classifier <- binomial_naive_bayes(X, Y, prior = prior, laplace = laplace, fit_prior = fit_prior,
                        binomial_n = binomial_n)
    result <- list(classifier = classifier, levels = levels)
    return(result)
}

### Object <- from above
### new data <- list(isolate_n = "aaaaaattttt", isolate_n = "aaaaaattttt", ...")
predict_balk <- function(object, newdata = NULL, snp_id = NULL, type = "prob"){
    binomial_n <- object$classifier$binomial_n
    levels <- object$levels
    if (is.null(newdata)){
        X <- newdata
    } else {
        rows <- names(newdata)
        sequences <- strsplit(unlist(newdata), split = "")
        if (! all(snp_id %in% colnames(levels))){
            warning("Not all snp_id are found in the training data, using subset of supplied data")
            ind <- which(snp_id %in% colnames(levels))
            snp_id <- snp_id[ind]
            sequences <- lapply(sequences, function(x){
                return(x[ind])
            })
        }
        columns <- snp_id
        available_snps <- levels[,match(columns, colnames(levels))]
        t_subject <- lapply(sequences, function(pat_all){
            res <- sapply(seq_along(pat_all), function(i){
                pat_one <- pat_all[i]
                snp_one <- available_snps[,i]
                return(transform_snp(pat_one, binomial_n, snp_one, get = "transformed"))
            })
            return(rbind(res))
        })
        X <- do.call(rbind, t_subject)
        colnames(X) <- columns
        rownames(X) <- rows
    }
    return(
        predict.binomial_naive_bayes(object$classifier, newdata = X, type = type))
}