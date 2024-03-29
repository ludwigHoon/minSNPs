% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addn_search.R
\name{parse_group_mcc}
\alias{parse_group_mcc}
\title{\code{parse_group_mcc}}
\usage{
parse_group_mcc(pattern, goi, MUST_HAVE_TARGET = TRUE)
}
\arguments{
\item{pattern}{the SNP profile for each samples}

\item{goi}{the samples belonging to the group of interest}

\item{MUST_HAVE_TARGET}{whether to force the profile to have at least 1 target profile
(the profile containing the most goi)}
}
\value{
the parsed group views
}
\description{
\code{parse_group_mcc} is used to group the sample according to SNPs profile and present in a table format
}
