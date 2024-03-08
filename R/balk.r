#' \code{get_binomial_tables}
#'
#' @description
#' \code{get_binomial_tables} is an internal function that returns a table of probability for binomial naive bayes.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param prob1 a matrix of probabilities
#' @return return table of probability for binomial naive bayes
get_binomial_tables <- function(prob1) {
    if (!is.matrix(prob1))
        stop("prob1 must be a matrix and an element of the binomial_naive_bayes object")
    n_tables <- nrow(prob1)
    vars <- rownames(prob1)
    tables <- lapply(seq_len(n_tables), function(i) {
        ith_row <- prob1[i, ]
        ith_tab <- as.table(rbind(1 - ith_row, ith_row))
        rownames(ith_tab) <- c("0", "1")
        ith_tab
    })
    names(tables) <- vars
    class(tables) <- "naive_bayes_tables"
    attr(tables, "cond_dist") <- stats::setNames(rep("Binomial", n_tables), vars)
    return(tables)
}

#' \code{binomial_naive_bayes}
#' 
#' @description
#' \code{binomial_naive_bayes} is an implementation of the binomial naive bayes algorithm.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param x a matrix with numeric: 0,1, up to binomial_n columns
#' @param y a factor or character or logical vector
#' @param prior a vector of prior probabilities
#' @param laplace a numeric value for Laplace smoothing
#' @param ... additional arguments
#' @importFrom stats setNames na.omit
#' @return return a binomial_naive_bayes object
binomial_naive_bayes <- function (x, y, prior = NULL, laplace = 1, ...)  {
    ### {Changes
    argl <- list(...)
    if (is.null(argl$fit_prior))
        argl$fit_prior <- FALSE
    if (is.null(argl$binomial_n))
        binomial_n <- 1
    else
        binomial_n <- argl$binomial_n
    ### Changes}
    if (!is.factor(y) & !is.character(y) & !is.logical(y))
        stop("binomial_naive_bayes(): y must be either a factor or character or logical vector", call. = FALSE)
    if (!is.factor(y))
        y <- factor(y)
    levels <- levels(y)
    nlev <- nlevels(y)
    vars <- colnames(x)
    class_x <- class(x)[1]
    use_Matrix <- FALSE
    if (!is.matrix(x) & !use_Matrix) {
        warning("binomial_naive_bayes(): x was coerced to matrix.", call. = FALSE)
        x <- as.matrix(x)
        if (mode(x) != "numeric")
            stop("binomial_naive_bayes(): x must be a matrix with with numeric {0,1,...binomial_n} columns.", call. = FALSE)
    }
     if (any(na.omit(x) < 0) | any(na.omit(x) > binomial_n))
            stop("binomial_naive_bayes(): x must be a matrix with with numeric {0,1,...binomial_n} columns.", call. = FALSE)

    if (nlev < 2)
        stop("binomial_naive_bayes(): y must contain at least two classes. ", call. = FALSE)
    if (is.null(vars))
        stop(paste0(match.call(), " failed, binomial_naive_bayes(): x must have unique column names.\n"), call. = FALSE)
    NAy <- anyNA(y)
    NAx <- anyNA(x)
    if (NAy) {
        na_y_bool <- is.na(y)
        len_na <- sum(na_y_bool)
        warning(paste0("binomial_naive_bayes(): y contains ", len_na, " missing",
                       ifelse(len_na == 1, " value", " values"), ". ",
                       ifelse(len_na == 1, "It is", "They are"),
                       " not included (also the corresponding rows in x) ",
                       "into the estimation process."), call. = FALSE)
        y <- y[!na_y_bool]
        x <- x[!na_y_bool, ]
    }
    if (NAx) {
        na_x <- is.na(x) * 1
        len_nax <- sum(na_x)
        warning(paste0("binomial_naive_bayes(): x contains ", len_nax, " missing",
                       ifelse(len_nax == 1, " value", " values"), ". ",
                       ifelse(len_nax == 1, "It is", "They are"),
                       " not included into the estimation process."), call. = FALSE)
    }
    y_counts <- stats::setNames(tabulate(y), levels)
    y_min <- y_counts < 1
    if (any(y_min))
        stop(paste0("binomial_naive_bayes(): y variable must contain at least ",
                    "one observation per class for estimation process.",
                    " Class ", paste0(levels[y_min], collapse =  ", "),
                    " has less than 1 observation."), call. = FALSE)
    if (is.null(prior)) {
        prior <- prop.table(y_counts)
    } else {
        if (length(prior) != nlev)
            stop(paste0("binomial_naive_bayes(): vector with prior probabilities must have ",
                        nlev, " entries"))
        prior <- stats::setNames(prior / sum(prior), levels)
    }
    ### {changes
    prob1 <- {
        xtemp <- (!is.na(x))*binomial_n
        params <- rowsum.default(xtemp, y) + laplace * 2
        x2 <- x
        x2[is.na(x2)] <- 0
        fc <- rowsum.default(x2, y) + laplace
        log(fc) - log(params)
    }
    ### changes}
    
    if (any(prob1 == 0)) {
        nempty <- length(which(prob1 == 0, arr.ind = TRUE)[ ,1])
        warning(paste0("binomial_naive_bayes(): there ", ifelse(nempty == 1, "is ", "are "),
                       nempty, " empty ", ifelse(nempty == 1, "cell ", "cells "),
                       "leading to zero estimates. Consider Laplace smoothing."), call. = FALSE)

    }
    structure(list(data = list(x = x, y = y), levels = levels, binomial_n = binomial_n,
                   laplace = laplace, prob1 = prob1, prior = prior, fit_prior = argl$fit_prior,
                   call = match.call()), class = "binomial_naive_bayes")
}

#' \code{predict.binomial_naive_bayes}
#' 
#' @description
#' \code{predict.binomial_naive_bayes} is an implementation of the predict method for the binomial naive bayes algorithm.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param object a binomial_naive_bayes object
#' @param newdata a matrix with numeric: 0,1,up to binomial_n columns
#' @param type a character string specifying the type of output: "class" or "prob"
#' @param ... additional arguments
#' @return return a factor or matrix of class probabilities
predict.binomial_naive_bayes <- function(object, newdata = NULL, type = c("class", "prob"), ...) {
    
    if (is.null(newdata))
        newdata <- object$data$x
    type <- match.arg(type)
    if (nrow(newdata) > 1 & anyNA(newdata)){
        results <- lapply(seq_len(nrow(newdata)), function(i){
            test_data <- rbind(newdata[i, ])
            test_data <- rbind(test_data[, !is.na(test_data)])
            predict.binomial_naive_bayes(object, newdata = test_data, type = type, ...)})
        if (type == "class")
            return(unlist(results))
        else
            return(do.call(rbind, results))
    }
    class_x <- class(newdata)[1]
    use_Matrix <- FALSE
    if (!is.matrix(newdata) & !use_Matrix)
        stop("predict.binomial_naive_bayes(): newdata must be numeric matrix with at least one row and two named columns.", call. = FALSE)
    if (is.matrix(newdata) & mode(newdata) != "numeric")
        stop("predict.binomial_naive_bayes(): newdata must be a numeric matrix.", call. = FALSE)

    
    lev <- object$levels
    n_lev <- length(lev)
    n_obs <- dim(newdata)[1L]
    ### {Changes
    binomial_n <- object$binomial_n
    prior <- object$prior
    fit_prior <- object$fit_prior
    prob1 <- object$prob1
    if (nrow(newdata) == 1)
        newdata <- rbind(newdata[,!is.na(newdata)])
    ### Changes}
    features <- colnames(newdata)[colnames(newdata) %in% colnames(prob1)]
    n_tables <- ncol(prob1)
    prob1 <- prob1[ ,features, drop = FALSE]
    n_features <- length(features)
    n_features_newdata <- ncol(newdata)

    if (n_features == 0) {
        warning(paste0("predict.binomial_naive_bayes(): no feature in newdata corresponds to ",
                       "features defined in the object. Classification is based on prior probabilities."), call. = FALSE)
        if (type == "class") {
            return(factor(rep(lev[which.max(prior)], n_obs), levels = lev))
        } else {
            return(matrix(prior, ncol = n_lev, nrow = n_obs, byrow = TRUE, dimnames = list(NULL, lev)))
        }
    }
    if (n_features < n_tables) {
        warning(paste0("predict.binomial_naive_bayes(): only ", n_features, " feature(s) in newdata could be matched ",
                       "with ", n_tables, " feature(s) defined in the object."), call. = FALSE)
    }
    if (n_features_newdata > n_features) {
        warning(paste0("predict.binomial_naive_bayes(): newdata contains feature(s) that could not be matched ",
                       "with (", n_features, ") feature(s) defined in the object. Only matching features are used for calculation."), call. = FALSE)
        newdata <- newdata[ ,features, drop = FALSE]
    }
    if (object$laplace == 0) {
        threshold <- 0.001
        eps <- 0
        prob1[prob1 <= eps] <- threshold
        prob1[prob1 >= (1 - eps)] <- 1 - threshold
    }

    lprob1 <- prob1
    ### {Changes
    lprob0 <- log(1 - exp(prob1))
   
    post <- tcrossprod(newdata, lprob1) + tcrossprod((binomial_n - newdata), lprob0)

    if (fit_prior){
        post <- post + log(prior)
    }
    ### Changes}
    if (type == "class") {
        return(factor(lev[max.col(post, "first")], levels = lev))
    }
    else {
        return(apply(post, 2, function(x) { 1 / rowSums(exp(post - x)) }))
    }
}

#' \code{print.binomial_naive_bayes}
#' 
#' @description
#' \code{print.binomial_naive_bayes} is an implementation of the print method for the binomial naive bayes algorithm.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param x a binomial_naive_bayes object
#' @param ... additional arguments
#' @return return a binomial_naive_bayes object
print.binomial_naive_bayes <- function (x, ...) {

    model <- "Binomial Naive Bayes"
    n_char <- getOption("width")
    str_left_right <- paste0(rep("=", floor((n_char - nchar(model)) / 2)),
                             collapse = "")
    str_full <- paste0(str_left_right, " ", model, " ",
                       ifelse(n_char %% 2 != 0, "=", ""), str_left_right)
    len <- nchar(str_full)
    l <- paste0(rep("-", len), collapse = "")
    cat("\n")
    cat(str_full, "\n", "\n", "Call:", "\n")
    print(x$call)
    cat("\n")
    cat(l, "\n", "\n")
    cat( "Laplace smoothing:", x$laplace)
    cat("\n")
    cat("\n")
    cat(l, "\n", "\n")
    cat(" A priori probabilities:", "\n")
    print(x$prior)
    cat("\n")
    ### {Changes
    cat(l, "\n", "\n")
    cat("Fit prior: ", x$fit_prior, "\n")
    cat(l, "\n", "\n")
    cat("Binomial N: ", x$binomial_n, "\n")
    ### Changes}
    cat(l, "\n", "\n")
    cat(" Tables:", "\n")
    tabs <- get_binomial_tables(x$prob1)
    n <- length(tabs)
    indices <- seq_len(min(5,n))
    tabs <- tabs[indices]
    print(tabs)
    if (n > 5) {
        cat("\n\n")
        cat("# ... and", n - 5, ifelse(n - 5 == 1, "more table\n\n", "more tables\n\n"))
        cat(l)
    }
    cat("\n\n")
}

#' \code{coef.binomial_naive_bayes}
#' 
#' @description
#' \code{coef.binomial_naive_bayes} is an implementation of the coef method for the binomial naive bayes algorithm.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param object a binomial_naive_bayes object
#' @param ... additional arguments
#' @return return a data frame of coefficients
#' @export
coef.binomial_naive_bayes  <- function(object, ...) {
    prob1 <- object$prob1
    levels <- object$levels
    nlev <- length(levels)
    m <- cbind(1 - exp(prob1), prob1)
    ind <- rep(seq_len(nlev), each = 2)
    m <- m[ ,ifelse(seq_along(ind) %% 2 != 0, ind, ind + nlev)]
    colnames(m) <- (paste0(rep(levels, each = 2), ":", c("0", "1")))
    as.data.frame(m)
}

#' \code{summary.binomial_naive_bayes}
#' 
#' @description
#' \code{summary.binomial_naive_bayes} is an implementation of the summary method for the binomial naive bayes algorithm.
#' modified from bernoulli_naive_bayes function in the naivebayes package
#' @param object a binomial_naive_bayes object
#' @param ... additional arguments
#' @return return a summary of the binomial_naive_bayes object
#' @export
summary.binomial_naive_bayes <- function(object, ...) {
    model <- "Binomial Naive Bayes"
    n_char <- getOption("width")
    str_left_right <- paste0(rep("=", floor((n_char - nchar(model)) / 2)),
                             collapse = "")
    str_full <- paste0(str_left_right, " ", model, " ",
                       ifelse(n_char %% 2 != 0, "=", ""), str_left_right)
    len <- nchar(str_full)
    l <- paste0(rep("-", len), collapse = "")
    cat("\n")
    cat(str_full, "\n", "\n")
    cat("- Call:", deparse(object$call), "\n")
    cat("- Laplace:", object$laplace, "\n")
    cat("- Classes:", nlevels(object$data$y), "\n")
    cat("- Samples:", length(object$data$y), "\n")
    cat("- Features:", nrow(object$prob1), "\n")
    cat("- Prior probabilities: \n")
    cat("    -", paste0(names(object$prior), ": ", round(object$prior, 4), collapse = "\n    - "))
    cat("\n\n")
    cat(l, "\n")
}