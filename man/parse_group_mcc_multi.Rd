% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addn_search.R
\name{parse_group_mcc_multi}
\alias{parse_group_mcc_multi}
\title{\code{parse_group_mcc_multi}}
\usage{
parse_group_mcc_multi(result, as_string = TRUE)
}
\arguments{
\item{result}{result from \code{find_optimised_snps}}

\item{as_string}{whether to return the result as string or data.frame}
}
\value{
Will return the grouped samples.
}
\description{
\code{parse_group_mcc_multi} is used to put samples according to SNP profile, and
put them into a table format.
}
