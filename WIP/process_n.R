# seed <- as.numeric(Sys.time())

devtools::load_all("./MinSNPs")
vivax <- read_fasta("Pv_Nature_communications_2018_VQSLOD3.min0.fasta")
seed <- 2021
set.seed(seed)

# //TODO add on other IUPAC.
translation <- list(
    M = c("A", "C"),
    R = c("A", "G"),
    W = c("A", "T"),
    S = c("C", "G"),
    Y = c("C", "T"),
    K = c("G", "T"),
    V = c("A", "C", "G"),
    H = c("A", "C", "T"),
    D = c("A", "G", "T"),
    B = c("C", "G", "T")
)


resolve_IUPAC_missing <- function(seqc, N_is_any_base = FALSE) { # nolint
    accepted_char <- c("A", "C", "T", "G")

    if (N_is_any_base) {
        translation[["N"]] <- accepted_char
    }

    for (p in seq_len(length(seqc[[1]]))) {

        nucleotides <- lapply(seqc, `[[`, p)
        to_modify <- which(!nucleotides %in% accepted_char)
        valid_replacements <- unlist(
            unique(nucleotides[which(nucleotides %in% accepted_char)]))
        cat("Position:", p, "\n", file = "replace.log", append = TRUE)
        for (t in to_modify) {
            print(seqc[[t]][p])
            print(translation[[seqc[[t]][p]]])
            if (seqc[[t]][p] == "N") {
                if (N_is_any_base) {
                    candidates <- valid_replacements
                } else {
                    candidates <- unique(c(translation[[seqc[[t]][p]]],
                    valid_replacements))
                }
            } else {
                candidates <- unique(translation[[seqc[[t]][p]]])
            }
            replacement <-
                sample(candidates, 1, replace = TRUE)
            cat("Replaced for,", names(seqc)[t], ":",
                seqc[[t]][p], "with", replacement, "\n",
                    file = "replace.log", append = TRUE)
            seqc[[t]][p] <- replacement
        }
    }
    return(seqc)
}

r_vivax <- resolve_IUPAC_missing(vivax)
write_fasta(r_vivax, "n_vivax.fasta")


resolve_missing_IUPAC <- function(seqc, bp=BiocParallel::SerialParam(),  #nolint
    dash_ignore=TRUE, accepted_char=c("A", "C", "T", "G"), ignore_case=TRUE) {

    nuc_seq <- bplapply(seq_len(length(seqc[[1]])), function(pos) {
        nucleotides <- lapply(seqc, `[[`, pos)
        which(nucleotides )
        return(seqc, `[[`, pos)
    }, BPPARAM = bp)
}