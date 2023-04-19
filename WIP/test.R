devtools::load_all("../MinSNPs/")
ref_genome <- minSNPs::read_fasta("../Mu50.fasta")
ref_meta <- read.csv("../out_ref_gene4.csv")

read_ref_meta1 <- read_ref_meta(ref_meta, ref_genome, FALSE, FALSE)

result <- get_fragments_isolates(read_ref_meta1, c(1, 15, 49, 50, 20651))
## Assume the SNPs in the following fasta positions are selected: 1 <1-239>, 15 <243-543>, (49, 50)<2866-3208>, 20651 <2878329-2878529>


# read_ref_meta2 <- read_ref_meta("../out_ref_gene4.csv", "../Mu50.fasta")

test_set_1 <- list(
    "fragment-1" = list(start = 1, end = 10),
    "fragment-2" = list(start = 11, end = 20),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_1)

test_set_2 <- list(
    "fragment-1" = list(start = 1, end = 11),
    "fragment-2" = list(start = 11, end = 20),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_2)

test_set_3 <- list(
    "fragment-1" = list(start = 1, end = 12),
    "fragment-2" = list(start = 11, end = 20),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_3)

test_set_4 <- list(
    "fragment-1" = list(start = 1, end = 10),
    "fragment-2" = list(start = 11, end = 21),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_4)

test_set_5 <- list(
    "fragment-1" = list(start = 1, end = 10),
    "fragment-2" = list(start = 11, end = 22),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_5)

test_set_6 <- list(
    "fragment-1" = list(start = 1, end = 11),
    "fragment-2" = list(start = 11, end = 22),
    "fragment-3" = list(start = 21, end = 30)
    )
collapse_fragments(test_set_6)