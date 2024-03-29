% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{merge_fasta}
\alias{merge_fasta}
\title{\code{merge_fasta}}
\usage{
merge_fasta(
  fasta_1,
  fasta_2,
  meta_1,
  meta_2,
  ref,
  method = "full",
  bp = BiocParallel::SerialParam(),
  ...
)
}
\arguments{
\item{fasta_1}{fasta read into memory to join}

\item{fasta_2}{fasta read into memory to join}

\item{meta_1}{meta file for `fasta_1` denoting all positions of SNPs
and position in reference genome}

\item{meta_2}{meta file for `fasta_2` denoting all positions of SNPs
and position in reference genome}

\item{ref}{name of the reference genome (needs to be in both fasta files)}

\item{method}{how to join the 2 fasta, currently supported methods are:
inner, full}

\item{bp}{the BiocParallel backend}

\item{...}{all other arguments}
}
\value{
Will return a list containing a merged FASTA and a meta.
}
\description{
\code{merge_fasta} is used to combine 2 fasta.
}
